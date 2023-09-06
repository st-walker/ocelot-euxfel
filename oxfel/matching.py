from __future__ import annotations
from dataclasses import dataclass

from ocelot.cpbd.beam import Twiss, get_envelope


@dataclass
class MismatchSummary:
    bmag_x: float
    bmad_y: float
    l2loss: float


class BacktrackingLinearMatcher:
    def __init__(self, navi: Navigator, parray0: ParticleArray,
                 goal_twiss: Twiss, quad_names: list[str],
                 twiss_function: Callable[[ParticleArray], Twiss] = get_envelope):
        self.navi = deepcopy(navi)
        self.parray0 = parray0
        self.quad_names = quad_names
        self.goal_twiss = goal_twiss
        self.twiss_function = twiss_function
        self.parray1 = None

    def _get_constraint(self) -> dict[str, dict[str, float]]:
        return {self.navi.lat.sequence[-1]: {"beta_x": self.goal_twiss.beta_x,
                                             "beta_y": self.goal_twiss.beta_y,
                                             "alpha_x": self.goal_twiss.alpha_x,
                                             "alpha_y": self.goal_twiss.alpha_y}}

    def l2loss(self):
        assert self.parray1 is not None
        twiss = self.twiss_function(self.parray1)
        gt = self.goal_twiss
        return np.sqrt((twiss.beta_x - gt.beta_x)**2
                       + (twiss.alpha_x - gt.alpha_x)**2
                       + (twiss.beta_y - gt.beta_y)**2
                       + (twiss.alpha_y - gt.alpha_y)**2)

    def bmags(self):
        assert self.parray1 is not None
        twiss = self.twiss_function(self.parray1)
        def calc_bmag(alpha, beta, alpha_design, beta_design):
            return 0.5 * (
            (beta / beta_design + beta_design / beta)
            + (beta * beta_design * ((alpha_design / beta_design) - (alpha / beta)) ** 2)
        )
        bmagx = calc_bmag(twiss.alpha_x, twiss.beta_x,
                          self.goal_twiss.alpha_x,
                          self.goal_twiss.beta_x)
        bmagy = calc_bmag(twiss.alpha_y, twiss.beta_y,
                          self.goal_twiss.alpha_y,
                          self.goal_twiss.beta_y)

        return bmagx, bmagy

    def track_forwards(self, quad_strengths=None) -> None:
        if quad_strengths is None:
            quad_strengths = self.quad_strengths()
        self.navi.go_to_start()
        _, parray1 = track(self.navi.lat, self.parray0.copy(),
                           navi=self.navi,
                           overwrite_progress=True,
                           print_progress=False)
        self.parray1 = parray1

    def _quad_instances_from_navi(self) -> list[Quadrupole]:
        quads = _get_element_instances_from_mlat(self.navi.lat, self.quad_names)
        return quads

    def quad_strengths(self):
        return [q.k1 for q in self._quad_instances_from_navi()]

    def _set_quads(self, strengths: list[float]) -> None:
        """Set the quadrupole strengths in the Navigator instance"""
        quads = self._quad_instances_from_navi()
        for quad, strength in zip(quads, strengths):
            quad.k1 = strength

    def initial_match(self) -> list[float]:
        """Do the first match before we do any particle tracking"""
        twiss0 = self.twiss_function(self.parray0)
        LOG.debug("Running initial matching")
        matched_strengths = match(self.navi.lat,
                                  self._get_constraint(),
                                  self._quad_instances_from_navi(),
                                  twiss0,
                                  verbose=True)
        return matched_strengths

    def match(self, n: int = 1, quad_strengths=None) -> list[float]:
        if quad_strengths is None:
            quad_strengths = self.initial_match()
        self._set_quads(quad_strengths)
        self.track_forwards(quad_strengths)
        for i in range(n):
            quad_strengths = self._backtrack_twiss_and_retrack()
            self._set_quads(quad_strengths)

        return quad_strengths

    def mismatch_summary(self):
        bmagx, bmagy = self.bmags()
        l2loss = self.l2loss()
        return MismatchSummary(bmagx, bmagy, l2loss)

    def rematch(self):
        if self.parray1 is None:
            self.track_forwards(quad_strengths)
        quad_strengths = self._backtrack_twiss_and_retrack()
        self._set_quads(quad_strengths)
        return quad_strengths

    def _backtrack_twiss_and_retrack(self) -> tuple[ParticleArray, list[float]]:
        # Start at end of lattice section and get twiss
        assert self.parray1 is not None
        twiss1 = self.twiss_function(self.parray1)
        LOG.debug(f"Using twiss function: {self.twiss_function}")
        LOG.debug(f"Twiss at end of line before backtracking: {twiss1.beta_x=}, {twiss1.beta_y}, {twiss1.alpha_x=}, {twiss1.alpha_y}")

        # Negate the alpha parameters for backtracking...
        twiss1.alpha_x *= -1
        twiss1.alpha_y *= -1
        # Reverse sequence and backtrack using the reversed Twiss parmaeters
        reversed_sequence = list(self.navi.lat.sequence)[::-1]
        backlat = MagneticLattice(reversed_sequence)
        backtracked_twiss0 = twiss(backlat, tws0=twiss1)[-1]
        bt0 = backtracked_twiss0
        LOG.debug(f"Twiss at start of line after backtracking: {bt0.beta_x=}, {bt0.beta_y}, {bt0.alpha_x=}, {bt0.alpha_y}")


        # Point the twiss forwards again...
        backtracked_twiss0.alpha_x *= -1
        backtracked_twiss0.alpha_y *= -1

        LOG.debug("Matching using backtracked Twiss Parameters")
        # Now we match forwards, using the backtrracked twiss
        matched_strengths = match(self.navi.lat,
                                  self._get_constraint(),
                                  self._quad_instances_from_navi(),
                                  backtracked_twiss0,
                                  verbose=False)
        # And track forwards using these new matched strengths
        self.track_forwards(matched_strengths)

        # Finally we return this new parray1 at the match point and
        # the new matched quad strengths we derived.
        return matched_strengths



def get_unary_twiss_function(twiss_type, **twissfnkwargs):
    if twiss_type.lower() == "projected":
        twiss_function = partial(get_envelope, **twissfnkwargs)
    elif twiss_type.lower() == "imax":
        twiss_function = partial(twiss_parray_slice, slice="Imax", **twissfnkwargs)
    elif twiss_type.lower() == "emax":
        twiss_function = partial(twiss_parray_slice, slice="Emax", **twissfnkwargs)
    elif callable(twiss_type):
        twiss_function = partial(twiss, **twissfnkwargs)
    else:
        raise ValueError(f"Unknown twiss type selected: {twiss_type}")
    return twiss_function


def match_with_backtracking(navi: Navigator,
                            parray0: ParticleArray,
                            twiss_goal: Twiss,
                            quad_names: list[str],
                            maxiter: int = 1,
                            match: str = "projected",
                            **twissfnkwargs):
    """Match a beam using linear optics and tracking to update the
    linear optics.  Only the horizontal Twiss parameters are
    considered (alpha_x, beta_x, alpha_y, beta_y) from the provided
    twiss_goal Twiss instance.

    :param navi:
    :param parray0:
    :param twiss_goal:
    :param quad_names:
    :param nmatch:
    :param match:
    :param twissfnkwargs:
    :return:

    """
    twiss_function = get_unary_twiss_function(match, **twissfnkwargs)
    matcher = BacktrackingLinearMatcher(navi, parray0, twiss_goal, quad_names, twiss_function)
    quad_strengths = matcher.match(maxiter)

    return quad_strengths, matcher.mismatch_summary()
