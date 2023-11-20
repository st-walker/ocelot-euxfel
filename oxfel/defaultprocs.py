from pathlib import Path
from importlib_resources import files

# from ocelot.utils.section_track import *
from oxfel.lattice import i1
from oxfel.lattice import i1d
from oxfel.lattice import l1
from oxfel.lattice import b1d
from oxfel.lattice import l2
from oxfel.lattice import b2d
from oxfel.lattice import l3
from oxfel.lattice import cl
from oxfel.lattice import tld


from ocelot.cpbd.physics_proc import *

# from oxfel.fel_track import FELSection, MachineSequence

from .processes import PlacedPhysicsProcess, PhysicsList


WAKES_PATH = files("oxfel.accelerator.wakes")

SmoothPar = 1000
LHE = 0 * 5000e-9  # GeV
WakeSampling = 1000
WakeFilterOrder = 10
CSRBin = 400
CSRSigmaFactor = 0.1
SCmesh = [63, 63, 63]
bISR = True
bRandomMesh = True


BEAM_SMOOTHER = SmoothBeam()
BEAM_SMOOTHER.mslice = 1000

LHE = 11000e-9 * 0.74 / 0.8  # GeV
WAKE_SAMPLING = 1000
WAKE_FILTER_ORDER = 10
CSR_SIGMA_FACTOR = 0.1
SC_MESH = [63, 63, 63]
# bISR = True
SC_RANDOM_MESH = True
CSR_N_BIN = 400

class DontTrackHereError(NotImplementedError):
    pass




class FELSection:
    pass


class G1(FELSection):
    def __init__(self):
        sequence = MachineSequence(i1.cell)[: i1.start_ocelot.id]
        super().__init__(sequence)

    def build_physics_processes(self):
        message = "Tracking in G1 is unsupported and should be left to ASTRA."
        return [(ExceptionalProc(message), self.start, self.stop)]

class A1(FELSection):
    def __init__(self):
        start_sim = i1.start_ocelot
        acc1_stop = i1.a1_sim_stop
        self.start = start_sim.id
        self.stop = acc1_stop.id

    def build_physics_processes(self):
        # init physics processes
        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = SCmesh
        sc.random_mesh = bRandomMesh

        sc2 = SpaceCharge()
        sc2.step = 1
        sc2.nmesh_xyz = SCmesh
        sc2.random_mesh = bRandomMesh

        wake = Wake()
        wake.wake_table = WakeTable(WAKES_PATH / "RF/wake_table_A1.dat")
        wake.factor = 1
        wake.step = 10
        wake.w_sampling = WakeSampling
        wake.filter_order = WakeFilterOrder

        smooth = SmoothBeam()
        smooth.mslice = SmoothPar
        # adding physics processes
        acc1_1_stop = i1.a1_1_stop
        physics_processes = [
            PlacedPhysicsProcess(smooth, self.start, self.start),
            PlacedPhysicsProcess(sc, self.start, acc1_1_stop.id),
            PlacedPhysicsProcess(sc2, acc1_1_stop.id, self.stop),
            PlacedPhysicsProcess(wake, i1.c_a1_1_1_i1.id, self.stop),
        ]

        return physics_processes


class AH1(FELSection):
    def __init__(self):
        ah1_start = i1.ah1_sim_start
        acc39_stop = i1.stlat_47_i1
        self.start = ah1_start.id
        self.stop = acc39_stop.id

    def build_physics_processes(self):
        # init physics processes
        sc = make_space_charge(step=1, nmesh_xyz=SCmesh, random_mesh=bRandomMesh)
        wake = make_wake(
            WAKES_PATH / "RF/wake_table_AH1.dat",
            factor=1,
            step=10,
            w_sampling=WakeSampling,
            filter_order=WakeFilterOrder,
        )
        physics_processes = [
            PlacedPhysicsProcess(sc, self.start, self.stop),
            PlacedPhysicsProcess(wake, i1.c3_ah1_1_1_i1.id, self.stop),
        ]

        return physics_processes


class LH(FELSection):
    def __init__(self):
        acc39_stop = i1.stlat_47_i1
        lhm_stop = i1.dump_csr_start  # for going in I1D
        self.start = acc39_stop.id
        self.stop = lhm_stop.id


    def build_physics_processes(self):
        csr = CSR()
        csr.sigma_min = physics_config.get_csr_sigma_min(self.start, self.stop)

        csr.traj_step = 0.0005
        csr.apply_step = 0.005
        sc = SpaceCharge()
        sc.step = 5
        sc.nmesh_xyz = SCmesh
        sc.random_mesh = bRandomMesh
        wake = Wake()
        wake.wake_table = WakeTable(WAKES_PATH / "RF/wake_table_TDS1.dat")
        wake.factor = 1
        wake.step = 10
        wake.w_sampling = WakeSampling
        wake.filter_order = WakeFilterOrder

        lh = LaserModulator()
        lh.dE = LHE
        lh.sigma_l = 300
        lh.sigma_x = 300e-6
        lh.sigma_y = 300e-6
        lh.z_waist = None

        physics_processes = [
            PlacedPhysicsProcess(sc, self.start, self.stop),
            PlacedPhysicsProcess(csr, self.start, self.stop),
            PlacedPhysicsProcess(wake, self.start, self.stop),
            PlacedPhysicsProcess(lh, "lh_start", "lh_stop", active=False)
        ]
        return physics_processes


class I1D(FELSection):
    def __init__(self):
        # init tracking lattice
        i1d_start = i1.dump_csr_start
        i1d_stop = i1d.otrc_64_i1d
        self.start = i1.dump_csr_start.id
        self.stop = i1d.ensec_66_i1d.id

    def build_physics_processes(self):
        # init tracking lattice
        # init physics processes
        sigma_min = physics_config.get_csr_sigma_min(self.start, self.stop)
        csr = make_csr(sigma_min=sigma_min, traj_step=0.0005, apply_step=0.005)
        sc = make_space_charge(step=5, nmesh_xyz=SCmesh, random_mesh=bRandomMesh)

        physics_processes = [
            PlacedPhysicsProcess(sc, self.start, self.stop),
            PlacedPhysicsProcess(csr, self.start, "BPMA.63.I1D"),
        ]
        return physics_processes


# def make_processes(

# class DL(FELSection):
#     def __init__(self):
#         cell = MachineSequence(i1.cell + l1.cell)
#         lh_stop_dl_start = i1.dump_csr_start.id  # "dogleg_start"
#         dogleg_stop = "dogleg_stop_bc0_start"
#         sequence = cell[lh_stop_dl_start:dogleg_stop]
#         super().__init__(sequence)

#     def build_physics_processes(self, physics_config):
#         csr = make_csr(
#             sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
#             traj_step=0.0005,
#             apply_step=0.005,
#             n_bin=CSR_N_BIN,
#         )
#         wake = make_wake(
#             WAKES_PATH / "mod_wake_0070.030_0073.450_MONO.dat",
#             factor=1,
#             w_sampling=WAKE_SAMPLING,
#             filter_order=WAKE_FILTER_ORDER,
#         )
#         sc = make_space_charge(step=25, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)

#         physics_processes = [
#             (csr, self.start, self.stop),
#             (sc, self.start, self.stop),
#             (wake, self.stop, self.stop),
#         ]
#         return physics_processes


# class BC0(FELSection):
#     def __init__(self):
#         l1_cell = MachineSequence(l1.cell)
#         dogleg_stop_bc0_start = "dogleg_stop_bc0_start"
#         bc0_stop_l1_start = "bc0_stop_l1_start"

#         sequence = l1_cell[dogleg_stop_bc0_start:bc0_stop_l1_start]
#         super().__init__(sequence)

#     def build_physics_processes(self, physics_config):
#         csr = make_csr(
#             step=1,
#             n_bin=CSR_N_BIN,
#             sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
#             traj_step=0.0005,
#             apply_step=0.001,
#         )
#         sc = make_space_charge(step=40, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)
#         physics_processes = [(sc, self.start, self.stop), (csr, self.start, self.stop)]
#         return physics_processes


# class L1(FELSection):
#     def __init__(self):
#         l1_cell = MachineSequence(l1.cell)
#         # Just after end of BC0
#         bc0_stop_l1_start = "bc0_stop_l1_start"
#         # This is just before the start of BC1.
#         a2_stop_bc1_start = "l1_stop_bc1_start"

#         sequence = l1_cell[bc0_stop_l1_start:a2_stop_bc1_start]
#         super().__init__(sequence)

#     def build_physics_processes(self, physics_config):
#         sc = make_space_charge(step=50, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)
#         wake = make_wake(
#             "wakes/RF/mod_TESLA_MODULE_WAKE_TAYLOR.dat",
#             factor=4,
#             step=100,
#             w_sampling=WAKE_SAMPLING,
#             filter_order=WAKE_FILTER_ORDER,
#         )
#         wake2 = make_wake(
#             "wakes/mod_wake_0078.970_0159.280_MONO.dat",
#             factor=1,
#             w_sampling=WAKE_SAMPLING,
#             filter_order=WAKE_FILTER_ORDER,
#         )

#         l1_first_cavity = "C.A2.1.1.L1"
#         l1_last_cavity = "C.A2.4.8.L1"

#         physics_processes = [
#             (BEAM_SMOOTHER, self.start, self.start),
#             (sc, self.start, self.stop),
#             (wake, l1_first_cavity, l1_last_cavity),
#             (wake2, self.stop, self.stop),
#         ]

#         return physics_processes


# class BC1(FELSection):
#     def __init__(self):
#         l1_cell = MachineSequence(l1.cell)
#         a2_stop_bc1_start = "l1_stop_bc1_start"
#         bc1_stop_l2_start = "bc1_stop_l2_start"

#         sequence = l1_cell[a2_stop_bc1_start:bc1_stop_l2_start]

#         super().__init__(sequence)

#     def build_physics_processes(self, physics_config):
#         # init physics processes
#         csr = make_csr(
#             step=1,
#             n_bin=CSR_N_BIN,
#             sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
#             traj_step=0.0005,
#             apply_step=0.001,
#         )
#         sc = make_space_charge(step=40, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)

#         return [(csr, self.start, self.stop), (sc, self.start, self.stop)]


# class B1D(FELSection):
#     def __init__(self):
#         cell = MachineSequence(l1.cell + b1d.cell)
#         bc1_stop_b1d_start = "STSEC.229.B1D"
#         sequence = cell[bc1_stop_b1d_start:]
#         super().__init__(sequence)

#     def build_physics_processes(self, physics_config):
#         csr = make_csr(
#             step=1,
#             n_bin=CSR_N_BIN,
#             sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
#             traj_step=0.0005,
#             apply_step=0.001,
#         )
#         sc = make_space_charge(step=50, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)

#         return [(csr, self.start, self.stop), (sc, self.start, self.stop)]


# class L2(FELSection):
#     def __init__(self):
#         cell = MachineSequence(l1.cell + l2.cell)
#         bc1_stop_l2_start = "bc1_stop_l2_start"
#         l2_stop_bc2_start = "l2_stop_bc2_start"
#         sequence = cell[bc1_stop_l2_start:l2_stop_bc2_start]
#         super().__init__(sequence)

#     def build_physics_processes(self, physics_config):
#         sc = make_space_charge(step=100, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)
#         wake = make_wake(
#             WAKES_PATH / "RF/mod_TESLA_MODULE_WAKE_TAYLOR.dat", factor=4 * 3, step=200
#         )
#         wake2 = make_wake(WAKES_PATH / "mod_wake_0179.810_0370.840_MONO.dat", factor=1)

#         first_cavity_l2 = "C.A3.1.1.L2"
#         last_cavity_l2 = "C.A5.4.8.L2"

#         return [
#             (BEAM_SMOOTHER, self.start, self.start),
#             (sc, self.start, self.stop),
#             (wake, first_cavity_l2, last_cavity_l2),
#             (wake2, self.stop, self.stop),
#         ]


# class BC2(FELSection):
#     def __init__(self):
#         l2_cell = MachineSequence(l2.cell)

#         l2_stop_bc2_start = "l2_stop_bc2_start"
#         bc2_stop_b2d_start = "bc2_stop_b2d_start"
#         # init tracking lattice
#         sequence = l2_cell[l2_stop_bc2_start:bc2_stop_b2d_start]
#         super().__init__(sequence)

#     def build_physics_processes(self, physics_config):
#         csr = make_csr(
#             step=1,
#             n_bin=CSR_N_BIN,
#             sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
#             traj_step=0.0005,
#             apply_step=0.001,
#         )
#         sc = make_space_charge(step=50, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)

#         return [(csr, self.start, self.stop), (sc, self.start, self.stop)]


# class B2D(FELSection):
#     def __init__(self):
#         cell = MachineSequence(l2.cell + b2d.cell)

#         bc2_stop_b2d_start = "bc2_stop_b2d_start"
#         sequence = cell[bc2_stop_b2d_start:]
#         super().__init__(sequence)

#     def build_physics_processes(self, physics_config):
#         csr = make_csr(
#             step=1,
#             n_bin=CSR_N_BIN,
#             sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
#             traj_step=0.0005,
#             apply_step=0.001,
#         )
#         sc = make_space_charge(step=50, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)

#         return [(csr, self.start, self.stop), (sc, self.start, self.stop)]


# class L3(FELSection):
#     def __init__(self):
#         bc2_stop = l2.tora_415_b2
#         acc6t26_stop = l3.stop_l3.id

#         # init tracking lattice
#         cell = MachineSequence(l2.cell + l3.cell + cl.cell)
#         sequence = cell[bc2_stop:acc6t26_stop]
#         super().__init__(sequence)

#     def build_physics_processes(self, physics_config):
#         # init physics processes
#         smooth = SmoothBeam()
#         smooth.mslice = SmoothPar

#         sc = SpaceCharge()
#         sc.step = 5
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh

#         wake = Wake()
#         wake.wake_table = WakeTable(WAKES_PATH / "RF/mod_TESLA_MODULE_WAKE_TAYLOR.dat")
#         wake.factor = 4 * 21
#         wake.step = 10
#         wake.w_sampling = 1000
#         wake_add = Wake()
#         wake_add.wake_table = WakeTable(
#             WAKES_PATH / "mod_wake_0391.350_1629.700_MONO.dat"
#         )
#         wake_add.factor = 1

#         app = PhaseSpaceAperture()
#         app.taumin = -5
#         app.taumax = 3

#         return [
#             (app, self.start, self.start),
#             (smooth, self.start, self.start),
#             (wake, l3.c_a6_1_1_l3.id, l3.c_a25_4_8_l3.id),
#             (sc, self.start, self.stop),
#             (wake_add, self.stop, self.stop),
#         ]


# class CL1(FELSection, SpontaneousRadiationSection):
#     def __init__(self):
#         acc6t26_stop = l3.stop_l3.id
#         collimator1_stop = cl.bpma_1746_cl.id

#         cell = MachineSequence(l3.cell + cl.cell)
#         sequence = cell[acc6t26_stop:collimator1_stop]
#         super().__init__(sequence)

#     def build_physics_processes(self, physics_config):
#         sc = make_space_charge(step=10, nmesh_xyz=SCmesh, random_mesh=bRandomMesh)
#         csr = make_csr(
#             traj_step=0.0005,
#             apply_step=0.001,
#             n_bin=CSRBin,
#             sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
#         )

#         physics_processes = [(sc, self.start, self.stop), (csr, self.start, self.stop)]

#         dipole_names = [
#             "BE.1678.CL",
#             "BL.1688.CL",
#             "BL.1695.CL",
#             "BE.1705.CL",
#             "BE.1714.CL",
#             "BL.1724.CL",
#             "BL.1731.CL",
#             "BE.1741.CL",
#         ]
#         if bISR:
#             physics_processes.extend(
#                 self.make_dipole_spontaneous_radiation(dipole_names)
#             )

#         return physics_processes


# class CL2(FELSection):
#     def __init__(self):
#         collimator1_stop = cl.bpma_1746_cl.id
#         collimator2_stop = cl.bpma_1783_cl.id
#         cell = MachineSequence(cl.cell)
#         sequence = cell[collimator1_stop:collimator2_stop]
#         super().__init__(sequence)

#     def build_physics_processes(self, physics_config):
#         sc = make_space_charge(step=1, nmesh_xyz=SCmesh, random_mesh=bRandomMesh)
#         return [(sc, self.start, self.stop)]


# class CL3(FELSection, SpontaneousRadiationSection):
#     def __init__(self):
#         collimator2_stop = cl.bpma_1783_cl.id
#         collimator3_stop = cl.ensec_1854_cl.id
#         # init tracking lattice
#         cell = MachineSequence(cl.cell)
#         sequence = cell[collimator2_stop:collimator3_stop]
#         super().__init__(sequence)

#     def build_physics_processes(self, physics_config):
#         sc = make_space_charge(step=10, nmesh_xyz=SCmesh, random_mesh=bRandomMesh)
#         csr = make_csr(
#             traj_step=0.0005,
#             apply_step=0.001,
#             n_bin=CSRBin,
#             sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
#         )
#         wake_add = make_wake(
#             WAKES_PATH / "mod_wake_1629.700_1831.200_MONO.dat", factor=1
#         )

#         processes = [
#             (csr, self.start, self.stop),
#             (sc, self.start, self.stop),
#             (wake_add, self.stop, self.stop),
#         ]
#         if bISR:
#             dipoles = [
#                 "BE.1786.CL",
#                 "BL.1796.CL",
#                 "BL.1803.CL",
#                 "BE.1813.CL",
#                 "BE.1822.CL",
#                 "BL.1832.CL",
#                 "BL.1839.CL",
#                 "BE.1849.CL",
#             ]
#             processes.extend(self.make_dipole_spontaneous_radiation(dipoles))

#         return processes


# class TLD(FELSection):
#     def __init__(self):
#         collimator3_stop = cl.ensec_1854_cl.id
#         cell = MachineSequence(cl.cell + tld.cell)
#         sequence = cell[collimator3_stop:]
#         super().__init__(sequence)

#     def build_physics_processes(self, physics_config):
#         return []


# class STN10(FELSection):
#     def __init__(self):
#         collimator3_stop = cl.ensec_1854_cl
#         stN10_stop = sase1.ensec_2235_t2
#         # stN10_stop = cl.ensub_1980_tl
#         # init tracking lattice
#         self.lattice = MagneticLattice(
#             cl.cell + tl34_sase1.cell + sase1.cell,
#             start=collimator3_stop,
#             stop=stN10_stop,
#             method=self.method,
#         )
#         # init physics processes
#         sc = SpaceCharge()
#         sc.step = 1
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh
#         wake_add = Wake()
#         wake_add.wake_table = WakeTable(
#             "accelerator/wakes/mod_wake_1831.200_2035.190_MONO.dat"
#         )
#         wake_add.factor = 1
#         wake_add1 = Wake()
#         wake_add1.wake_table = WakeTable(
#             "accelerator/wakes/mod_wake_2035.190_2213.000_MONO.dat"
#         )
#         wake_add1.factor = 1

#         self.add_physics_process(sc, start=collimator3_stop, stop=stN10_stop)
#         self.add_physics_process(
#             wake_add, start=collimator3_stop, stop=collimator3_stop
#         )
#         self.add_physics_process(wake_add1, start=stN10_stop, stop=stN10_stop)
