import os
from pathlib import Path
import re
from typing import Optional, Any, Union, TypeVar
from dataclasses import dataclass, field
from enum import Enum, auto
from collections import deque
import logging
from functools import reduce
import pickle
from importlib_resources import files

import numpy as np
import ocelot.cpbd.elements as elements
from ocelot.cpbd.beam import ParticleArray
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.beam import Twiss, get_envelope
from ocelot.cpbd.latticeIO import LatticeIO
from ocelot.cpbd.elements.optic_element import OpticElement
import pandas as pd
import toml

from .fel_track import EuXFELSimConfig, Linac
from .optics import START_SIM, MATCH_37, MATCH_52, INJECTOR_MATCHING_QUAD_NAMES
from .longlist import make_default_longlist
from .astra import load_reference_0320_100k_distribution


logging.basicConfig()
LOG = logging.getLogger(__name__)
LOG.setLevel(logging.INFO)

AnyTupT = tuple[Any, ...]
ElementT = TypeVar("ElementT", bound=OpticElement)

# Filtering these not just for clutter reasons but also they often
# overlap with other beamline elements, meaning in effect we have to
# choose..  and these elements i think are more for engineers than
# physicists (guess).
SKIP_TYPE = ["BENDMARK", "RF"]
SKIP_CLASS = []
SKIP_GROUP = [
    "CRYO",
    "VACUUM",
    "MOVER",  # "FASTKICK"
]

TRACKING_CONF_NAME = "tracking-matched-conf.toml"
REAL_CONF_NAME = "real-matched-conf.toml"


def get_bunchsizerc_path(beam_destination):
    # e.g. i1d, b2d, tld, etc..
    outdir = files("oxfel.accelerator.lattice")
    return outdir / "{}-real-optics.pcl".format(beam_destination)


class MarkerPositionType(Enum):
    BEFORE = auto()
    AFTER = auto()


@dataclass
class SMarker:
    name: str
    s: float


@dataclass
class RelativeMarker:
    name: str
    reference: str
    position: MarkerPositionType


@dataclass
class LatticeSection:
    start_name1: str
    stop_name1: str
    name: str
    markers: list[Union[SMarker, RelativeMarker]] = field(default_factory=list)


class UnknownLongListElement(ValueError):
    pass


class MarkerPlacer:
    def __init__(self, new_markers: list[Union[SMarker, RelativeMarker]]):
        self.s_markers = [x for x in new_markers if isinstance(x, SMarker)]
        # Don't use a dictionary cos the referenced element can be used multiple times...
        self.rel_markers = [x for x in new_markers if isinstance(x, RelativeMarker)]
        # self.referenced_markers = {x.reference for x in rel_markers}

    def build_element_sequence_with_markers(self, oelement) -> list:
        result = [oelement]
        if not self.rel_markers:
            return result
        seen = []

        result = deque(result)
        for marker in self.rel_markers:
            if marker.reference != oelement.id:
                continue

            omarker = elements.Marker(eid=marker.name)
            if marker.position is MarkerPositionType.BEFORE:
                result.appendleft(omarker)
            elif marker.position is MarkerPositionType.AFTER:
                result.append(omarker)
            else:
                raise ValueError(f"Unknown Marker position type: {marker.position}")

        for seen_marker in seen:
            self.rel_markers.remove(seen_marker)

        LOG.debug(f"Appending: {result}")

        return list(result)

    def s_markers_in_interval(
        self, s_start: float, s_stop: float
    ) -> list[tuple[elements.Marker, float]]:
        result = []
        for s_marker in self.s_markers:
            s = s_marker.s
            if s >= s_start and s < s_stop:
                result.append((elements.Marker(eid=s_marker.name), s))
        return result


class LongListConverter:
    ROUND_NDIGITS = 6  # For rounding to make drifts and other lenghts look nice.

    def __init__(self, fname: os.PathLike, extra_properties=None):
        self.df = pd.read_excel(fname, sheet_name="LONGLIST", skiprows=[1])
        # So that the columns are well formed python names...
        self.df = self.df.rename(columns={"E1/LAG": "E1_LAG", "E2/FREQ": "E2_FREQ"})
        self.extra_properties = extra_properties if extra_properties else {}
        # # to have unique drift names...
        self.drift_counter = 0

    @property
    def twiss0(self) -> Twiss:
        start = self.df.iloc[0]
        t = Twiss()
        t.beta_x = start.BETX
        t.beta_y = start.BETY
        t.alpha_x = start.ALFX
        t.alpha_y = start.ALFY
        t.E = start.ENERGY
        return t

    def _round(self, value: float) -> float:
        return round(value, self.ROUND_NDIGITS)

    def convert_sections(self, sections: list[LatticeSection]):
        self.drift_counter = 0
        result = {}
        for section in sections:
            result[section.name] = self.convert_section(section)
        return result

    def _filter_on_bad_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        t_mask = [df.TYPE == t for t in SKIP_TYPE]
        c_mask = [df.CLASS == c for c in SKIP_CLASS]
        g_mask = [df.GROUP == g for g in SKIP_GROUP]

        full_bad_mask = reduce(np.logical_or, [*t_mask, *c_mask, *g_mask])

        return df[~full_bad_mask]

    def _filter_bad_entries(self, df: pd.DataFrame) -> pd.DataFrame:
        neg_s = df[df.S < 0]
        for tup in neg_s.itertuples():
            LOG.warning(f"Dropping element with negative S: {tup.Index=}, {tup.NAME1}")
        df = df[df.S >= 0]
        df = self._filter_on_bad_columns(df)

        # Skip markers that are inside other elements.  causes problems later.
        # Define a "marker" as simply any element that has zero length
        thin_elements = df[df.LENGTH == 0]
        starts = df.S - df.LENGTH * 0.5
        stops = df.S + df.LENGTH * 0.5

        bad_rows = []
        for tup in thin_elements.itertuples():
            is_inside = ((tup.S > starts) & (tup.S < stops)).any()
            if is_inside:
                LOG.warning(
                    f"Dropping bad row inside other element: {tup.Index=}, {tup.NAME1}"
                )
                bad_rows.append(tup)

        bad_row_indices = [x.Index for x in bad_rows]
        return df.drop(index=bad_row_indices)

    def convert_section(self, pysec: LatticeSection) -> list[ElementT]:
        # Do the slicing.
        section_df = self.df.set_index("NAME1").loc[
            pysec.start_name1 : pysec.stop_name1
        ]
        section_df = section_df.reset_index()  # put NAME1 back as a column again.
        # section_df = section_df[section_df.S >= 0] # Skip negative S elements...

        section_df = self._filter_bad_entries(section_df)

        # to be returned:
        sequence = []
        marker_placer = MarkerPlacer(pysec.markers)

        # the end position of the last element in S
        running_s_end = section_df.iloc[0].S - 0.5 * section_df.iloc[0].LENGTH

        for row in section_df.itertuples():
            # Convert to OCELOT element
            oelement = self.dispatch(row)
            # Can't use row.LENGTH because it's wrong for example in UNDU.
            # S Is the same/correct though.
            oelement_start_s = self._round(row.S - oelement.l * 0.5)
            oelement_with_markers = self._attach_markers_to_element(
                running_s_end, oelement_start_s, oelement, marker_placer
            )
            sequence.extend(oelement_with_markers)
            running_s_end += sum([x.l for x in oelement_with_markers])

        return sequence

    def _attach_markers_to_element(
        self, s_start: float, s_stop: float, oelement, marker_placer: MarkerPlacer
    ) -> list:
        """s_start, s_stop define region where were look to try to insert markers"""
        assert s_stop >= s_start

        # Attach any markers either side of the element (based on name)
        oelement_with_markers = marker_placer.build_element_sequence_with_markers(
            oelement
        )

        # See if there are any markers placed at S
        markers_with_s = marker_placer.s_markers_in_interval(s_start, s_stop)
        if not markers_with_s:
            drift_length = s_stop - s_start
            if drift_length > 0:
                drift = self._next_drift(l=drift_length)
                oelement_with_markers = [drift] + oelement_with_markers
            return oelement_with_markers

        # Place markers
        for marker, s in markers_with_s:
            # Drift, Marker
            drift_length = self._round(s - s_start)
            if drift_length > 0:
                oelement_with_markers.append(self._next_drift(l=drift_length))
            oelement_with_markers.append(marker)
            s_start = s

        oelement_with_markers.append(self._next_drift(s_stop - s))

        return oelement_with_markers

    def _next_drift(self, l: float) -> elements.Drift:
        drift = elements.Drift(l=l, eid=f"D_{self.drift_counter}")
        self.drift_counter += 1
        return drift

    def _apply_manual_changes(self, element: ElementT) -> None:
        try:
            properties = self.extra_properties[element.id]
        except KeyError:
            return
        else:
            for property_name, value in properties.items():
                setattr(element, property_name, value)
                if property_name == "nperiods" or property_name == "lperiod":
                    element.l = element.nperiods * element.lperiod

    def dispatch(self, row: AnyTupT) -> ElementT:
        ocelot_class_name = row.GROUP.lower()
        try:
            oelement = getattr(self, f"convert_{ocelot_class_name}")(row)
            self._apply_manual_changes(oelement)
            return oelement
        except AttributeError:
            raise UnknownLongListElement(
                "Unknown element to be converted:"
                f"{row.NAME1=}, {row.GROUP=}, {row.CLASS=}, {row.TYPE=}"
            )

    def convert_magnet(
        self, tup: AnyTupT
    ) -> Union[
        elements.Quadrupole,
        elements.Hcor,
        elements.Vcor,
        elements.SBend,
        elements.Solenoid,
        elements.Sextupole,
    ]:
        common_kw = {
            "eid": tup.NAME1,
            "l": tup.LENGTH,
        }

        if tup.TILT != 0.0:
            common_kw["tilt"] = tup.TILT
        if tup.LENGTH != 0.0:
            common_kw["l"] = tup.LENGTH

        assert tup.GROUP in {"RAMPKICK", "MAGNET", "FASTKICK", "FBKICK"}, tup.GROUP

        if tup.CLASS == "QUAD":
            return elements.Quadrupole(k1=tup.STRENGTH / tup.LENGTH, **common_kw)
        elif tup.CLASS == "HKIC":
            return elements.Hcor(angle=tup.STRENGTH, **common_kw)
        elif tup.CLASS == "VKIC":
            return elements.Vcor(angle=tup.STRENGTH, **common_kw)
        elif tup.CLASS == "SBEN":
            return elements.SBend(
                angle=tup.STRENGTH, e1=tup.E1_LAG, e2=tup.E2_FREQ, **common_kw
            )
        elif tup.CLASS == "SOLE":
            return elements.Solenoid(**common_kw)
        elif tup.CLASS == "SEXT":
            return elements.Sextupole(k2=tup.STRENGTH / tup.LENGTH, **common_kw)
        elif tup.CLASS == "RBEN":
            return elements.RBend(
                angle=tup.STRENGTH, e1=tup.E1_LAG, e2=tup.E2_FREQ, **common_kw
            )

        raise UnknownLongListElement(tup)

    def convert_pmagnet(self, row: AnyTupT) -> elements.RBend:
        assert row.CLASS == "RBEN"
        # Not sure what this is but in the longlist it has zero strength...
        # So just assert that it will always have zero strength and if this changes
        # Then address this in the future.
        assert row.STRENGTH == 0.0
        return elements.RBend(angle=row.STRENGTH, l=row.LENGTH, eid=row.NAME1)

    def convert_undu(self, row: AnyTupT) -> elements.Undulator:
        assert row.GROUP == "UNDU"
        assert row.CLASS == "UNDULATOR"
        # Extract sequences of numbers from the TYPE
        type_numbers = re.findall(r"\d+", row.TYPE)
        lperiod = float(type_numbers[0]) * 1e-3

        undulator = elements.Undulator(lperiod=lperiod, nperiods=1, eid=row.NAME1)
        # have to explicitly set the length here, otherwise it is 0.0
        # cos the class invariant is not properly preserved by the
        # class itself.
        undulator.l = undulator.lperiod * undulator.nperiods
        return undulator

    def convert_cavity(self, row: AnyTupT) -> Union[elements.Cavity, elements.TDCavity]:
        assert row.GROUP == "CAVITY"
        common_kw = {
            "eid": row.NAME1,
            "l": row.LENGTH,
            "v": row.STRENGTH * 1e-3,
            "phi": row.E1_LAG * 360,
            "freq": row.E2_FREQ * 1000000,
        }
        if row.TILT:
            common_kw["tilt"] = row.TILT

        if row.TYPE == "C" or row.TYPE == "C3":
            return elements.Cavity(**common_kw)
        elif row.TYPE == "TDSA" or row.TYPE == "TDSB":
            return elements.TDCavity(**common_kw)

        raise UnknownLongListElement(row)

    def convert_diag(self, row: AnyTupT) -> Union[elements.Marker, elements.Monitor]:
        assert row.GROUP == "DIAG"
        if row.CLASS == "INSTR":
            return elements.Marker(eid=row.NAME1)
        elif row.CLASS == "MONI":
            return elements.Monitor(eid=row.NAME1, l=row.LENGTH)
        elif row.CLASS == "CM":
            return elements.Marker(eid=row.NAME1)
        raise UnknownLongListElement(row)

    def convert_mark(self, row: AnyTupT) -> elements.Marker:
        assert row.GROUP == "MARK"
        return elements.Marker(eid=row.NAME1)

    convert_rampkick = convert_magnet
    convert_fastkick = convert_magnet
    convert_fbkick = convert_magnet


def longlist_to_ocelot(
    fname: os.PathLike,
    ftoml: Optional[os.PathLike] = None,
    outdir: Optional[os.PathLike] = "./",
    match: bool = False,
) -> None:
    if ftoml is None:
        config = get_default_config()
    else:
        config = toml.load(ftoml)

    outdir = Path(outdir)

    # Parse the config dictionary
    sections, extras = _parse_config_dict(config)

    llcv = LongListConverter(fname, extra_properties=extras)

    sequences = llcv.convert_sections(sections)

    for i, (name, sequence) in enumerate(sequences.items()):
        if i == 0:
            twiss0 = llcv.twiss0
        else:
            twiss0 = None

        outf = outdir / f"{name.lower()}.py"
        LatticeIO.save_lattice(
            MagneticLattice(sequence),
            file_name=outf,
            remove_rep_drifts=False,
            tws0=twiss0,
        )
        print("Written", outf)

        real_config = make_felconfig_design_to_real(config["extras"]["real"])
    real_matched_conf = generate_real_i1_matched_config(
        sequences["I1"], llcv.twiss0, real_config
    )

    with (outdir / REAL_CONF_NAME).open("w") as f:
        toml.dump(
            {"components": real_matched_conf.components},
            f,
            encoder=toml.TomlNumpyEncoder(),
        )
        print(f"Wrote {f.name}")


def generate_real_i1_matched_config(i1_sequence, twiss0, realconfig):
    # i1_dummy_section = FELSection(i1_sequence)
    just_injector = Linac(i1_sequence, twiss0, felconfig=realconfig)
    match_52_twiss_constraint = make_default_longlist().get_optics_constraint(MATCH_52)
    real_matched_conf = just_injector.match(
        just_injector.twiss0,
        elements=INJECTOR_MATCHING_QUAD_NAMES,
        constraints=match_52_twiss_constraint,
    )
    return real_matched_conf


def make_felconfig_design_to_real(dconf: dict) -> EuXFELSimConfig:
    result = EuXFELSimConfig()

    try:
        lh_angle = dconf["hlc"]["LH"]["angle"]
    except KeyError:
        lh_angle = None

    result.controls["lh"].angle = lh_angle

    try:
        bc0_angle = dconf["hlc"]["BC0"]["angle"]
    except KeyError:
        bc0_angle = None
    result.controls["bc0"].angle = bc0_angle

    try:
        bc1_angle = dconf["hlc"]["BC1"]["angle"]
    except KeyError:
        bc1_angle = None
    result.controls["bc1"].angle = bc1_angle

    try:
        bc2_angle = dconf["hlc"]["BC2"]["angle"]
    except KeyError:
        bc2_angle = None
    result.controls["bc2"].angle = bc2_angle

    result.components = dconf["components"]
    return result


def _parse_new_markers_dict(dconf: dict) -> list[Union[RelativeMarker, SMarker]]:
    markers = []

    for new_marker_name, position_info in dconf.items():
        try:
            defn = SMarker(new_marker_name, position_info["s"])
        except KeyError:
            pass
        else:
            markers.append(defn)
            continue

        try:
            defn = RelativeMarker(
                new_marker_name,
                position_info["ref"],
                MarkerPositionType[position_info["location"].upper()],
            )
        except KeyError:
            pass
        else:
            markers.append(defn)
            continue

    return markers


def _parse_config_dict(dconf: dict) -> tuple[LatticeSection, dict]:
    sections = []
    for section_name, info in dconf["sections"].items():
        try:
            new_markers = _parse_new_markers_dict(info["new_markers"])
        except KeyError:
            new_markers = []

        section = LatticeSection(
            start_name1=info["start_name1"],
            stop_name1=info["stop_name1"],
            name=section_name,
            markers=new_markers,
        )
        sections.append(section)

    extras = dconf["extras"]["design"]

    return sections, extras


def get_default_config() -> dict:
    return toml.load(files("oxfel.accelerator.lattice") / "conversion-config.toml")


def match_real_injector():
    from oxfel import cat_to_i1d

    fel = cat_to_i1d(model_type="real")

    parray032 = load_reference_0320_100k_distribution()

    parray37 = fel.track(parray032, start="start_ocelot", stop=MATCH_37)
    twiss37 = get_envelope(parray37)

    ll = make_default_longlist()
    match_52_twiss_constraint = ll.get_optics_constraint("MATCH.52.I1")

    linear_matched_conf = fel.match(
        twiss37,
        start=MATCH_37,
        stop=MATCH_52,
        constraints=match_52_twiss_constraint,
        elements=INJECTOR_MATCHING_QUAD_NAMES,
        verbose=True,
    )

    tracking_matched_conf = fel.match_beam(
        parray37,
        start=MATCH_37,
        stop=MATCH_52,
        elements=INJECTOR_MATCHING_QUAD_NAMES,
        constraints=match_52_twiss_constraint,
        felconfig=linear_matched_conf,
    )

    toml.load(files("oxfel.accelerator.lattice") / TRACKING_CONF_NAME)
    outdir = files("oxfel.accelerator.lattice")

    with (outdir / TRACKING_CONF_NAME).open("w") as f:
        toml.dump(
            {"components": tracking_matched_conf.components},
            f,
            encoder=toml.TomlNumpyEncoder(),
        )
        print(f"Wrote {f.name}")


def update_bunch_size_data(name, fel: Linac, parray032: ParticleArray):
    _, twiss = fel.track_optics(parray032, start=START_SIM)
    outdir = files("oxfel.accelerator.lattice")
    fpath = get_bunchsizerc_path(name)
    twiss.to_pickle(fpath)
    print(f"Wrote {fpath}")
