import numpy as np
import latdraw
from textwrap import dedent
from collections import Counter
import latdraw

from ocelot.cpbd.elements import (Marker, Drift, Solenoid, Hcor, Vcor,
                                  Monitor, Cavity, Undulator,
                                  SBend,
                                  Quadrupole,
                                  TDCavity)

DIPOLE_WIDTH = 0.5 # 50cm

SKIPPABLE_ELEMENTS = (Marker, Drift, Solenoid, Hcor, Vcor, Monitor, Cavity, Undulator, TDCavity)

class OcelotToAstra:
    SECTION_DELIMITER = "\n\n! ------\n\n"

    def __init__(self):
        self.elements_counter = Counter()
        self.current_rotation_theta = 0.0

    def convert(self, sequence, skippable=None, chicane=None):
        self.elements_counter = Counter()

        if skippable is None:
            skippable = []

        if chicane is None:
            chicane = []

        survey = latdraw.lattice_from_ocelot(sequence).survey()

        self.elements_counter["dipole"] += 1
        self.elements_counter["quadrupole"] += 1
        dipoles = []
        quadrupoles = []
        for esurvey, element in zip(survey.itertuples(), sequence):
            if element.id in skippable:
                print("Explicitly skipped element", element, esurvey)
            elif isinstance(element, SBend) and element.id in chicane:
                dipoles.append(self.make_chicane_dipole(element, esurvey))
                self.elements_counter["dipole"] += 1
            elif isinstance(element, SBend):
                dipoles.append(self.make_rbend(element, esurvey))
                self.elements_counter["dipole"] += 1
            elif isinstance(element, Quadrupole):
                # from IPython import embed; embed()
                quadrupoles.append(self.make_quadrupole(element, esurvey))
                self.elements_counter["quadrupole"] += 1
            elif isinstance(element, SKIPPABLE_ELEMENTS):
                print("Skipping", element)
            else:
                raise RuntimeError(element)

        output_section = self.make_output(esurvey.z,
                                          survey[survey["name"] == "OTRC.64.I1D"],
                                          ["OTRC.64.I1D"])

        return (output_section + 
                self.SECTION_DELIMITER
                + self.quadrupole_section(quadrupoles)
                + self.SECTION_DELIMITER
                + self.dipole_section(dipoles))

    def make_output(self, z_stop, screen_survey, screen_comments):
        output = dedent(f"""\
        &OUTPUT
        ZSTART=0
        ! I1D DUMP
        ZSTOP={z_stop}
        Zphase=40
        EmitS=T
        PhaseS=T
        TrackS=F
        RefS=F
        TcheckS=F
        CathodeS=F
        high_res=T
        """)

        # from IPython import embed; embed()

        for i, (ssurvey, comment) in enumerate(zip(screen_survey.itertuples(), screen_comments)):
            output += dedent(f"""\
            ! {comment}
            SCREEN({i})={ssurvey.z}
            Scr_xrot({i})=-{ssurvey.theta}""")

        return output + "\n/"

    def section_heading(self, name):
        return f"{self.SECTION_DELIMITER} {name}"

    def make_chicane_dipole(self, dipole, dsurvey):
        name = dipole.id
        l = dsurvey.length
        end_z = dsurvey.z
        hw = 0.5 * DIPOLE_WIDTH
        start_z = end_z - l
        d1 = (+hw, start_z)
        d2 = (-hw, start_z)
        d3 = (+hw, start_z + l)
        d4 = (-hw, start_z + l)

        bending_radius = l / dipole.angle

        dipole_index = self.elements_counter["dipole"]

        definition = dedent(f"""\
        ! Dipole: {name}
        D_Type({dipole_index})="horizontal"
        D1({dipole_index})={d1}
        D2({dipole_index})={d2}
        D3({dipole_index})={d3}
        D4({dipole_index})={d4}
        D_radius({dipole_index})={bending_radius}
        D_Gap(1,{dipole_index})=0.0000001
        D_Gap(2,{dipole_index})=0.0000001
        """)
        return definition

    def make_rbend(self, dipole, dsurvey):
        name = dipole.id
        l = dsurvey.length
        end_z = dsurvey.z
        hw = 0.5 * DIPOLE_WIDTH
        start_z = end_z - l
        d1 = (+hw, start_z)
        d2 = (-hw, start_z)
        d3 = (+hw, start_z + l)
        d4 = (-hw, start_z + l)

        bending_radius = l / dipole.angle

        dipole_index = self.elements_counter["dipole"]

        definition = dedent(f"""\
        ! Dipole: {name}
        D_Type({dipole_index})="horizontal"
        D1({dipole_index})={d1}
        D2({dipole_index})={d2}
        D3({dipole_index})={d3}
        D4({dipole_index})={d4}
        D_radius({dipole_index})={bending_radius}
        D_Gap(1,{dipole_index})=0.0000001
        D_Gap(2,{dipole_index})=0.0000001
        """)

        placement_string = f"D_xrot({dipole_index})=-{dipole.angle / 2.0}\n"

        return definition + placement_string

    def make_quadrupole(self, quadrupole, qsurvey):
        name = quadrupole.id
        l = quadrupole.l
        bore = 0.0000001
        k1 = quadrupole.k1
        iquad = self.elements_counter["quadrupole"]

        end = np.array([qsurvey.x, qsurvey.y, qsurvey.z])
        tangent = qsurvey.zlocal
        xmid, ymid, zmid = end - 0.5 * l * tangent

        placement_string = f"Q_pos({iquad})={zmid}\n"
        if abs(xmid) >  1e-6:
            placement_string += f"Q_xoff({iquad})={xmid}\n"
        if abs(ymid) > 1e-6:
            placement_string += f"Q_yoff({iquad})={ymid}\n"
        if qsurvey.theta:
            placement_string += f"Q_xrot({iquad})=-{qsurvey.theta}\n"

        definition = dedent(f"""\
        ! Quadrupole: {name}
        Q_length({iquad})={l}
        Q_K({iquad})={k1}
        Q_Bore({iquad})={bore}
        """)

        return definition + placement_string

    def quadrupole_section(self, quadrupole_definitions):
        preamble = "&QUADRUPOLE\nLoop=.F\nLQUAD=.T\n"
        postamble = "/"
        return preamble + "\n".join(quadrupole_definitions) + postamble

    def dipole_section(self, dipole_definitions):
        preamble = "&DIPOLE\nLoop=.F\nLDipole=.T\n"
        postamble = "/"
        return preamble + "\n".join(dipole_definitions) + postamble

def injector_to_astra(injector_sequence, outfile):
    model_string = OcelotToAstra().convert(injector_sequence, skippable=["QLN.23.I1", "QLS.23.I1"],
                                           chicane=["BL.48I.I1", "BL.48II.I1", "BL.50I.I1", "BL.50II.I1"])

    from IPython import embed; embed()

# def print_all_of_lh_chicane_human_readable():
#     lat = MagneticLattice(i1.cell, start=i1.qi_47_i1, stop=i1.qi_50_i1)

#     brho = 3.3356 * 130e-3
#     print(f"{brho=}")

#     # from IPython import embed; embed()

#     previous_element = None
#     summed_drift = 0
#     for element in flatten(lat.sequence):
#         if element.l == 0:
#             continue
#         if isinstance(element, Drift):
#             summed_drift += element.l
#             continue
#         elif summed_drift != 0:
#             print("drift, l =",summed_drift)
#             summed_drift = 0

#         if isinstance(element, SBend):
#             rho = element.l / element.angle
#             b = brho / rho
#             print(f"{element.id},{element.angle}, {element.l}, {b=}T")
#         elif isinstance(element, Undulator):
#             from scipy.constants import e, m_e, c
#             from math import pi
#             k = element.Kx
#             b = k * 2 * pi * m_e * c  / (e * element.lperiod)
#             print(element.id, f"{b=}")

#         elif isinstance(element, (Hcor, Vcor)):
#             print(f"{element.id=}, {element.l=}")
#         else:
#             print(element, "!!!!!!!!!!")

