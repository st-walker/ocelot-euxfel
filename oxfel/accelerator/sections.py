from pathlib import Path
from importlib_resources import files

from ocelot.utils.section_track import *
from .lattice import i1
from .lattice import i1d
from .lattice import l1
from .lattice import b1d
from .lattice import l2
from .lattice import b2d
from .lattice import l3
from .lattice import cl
from .lattice import tld


from ocelot.cpbd.physics_proc import *

# from oxfel.fel_track import FELSection, MachineSequence


class SpontaneousRadiationSection:
    def make_dipole_spontaneous_radiation(self, dipole_names):
        for dipole_name in dipole_names:
            dipole = self.sequence[dipole_name]
            length = dipole.l
            angle = dipole.angle
            radius = length / np.sin(angle)
            proc = SpontanRadEffects(radius=radius, type="dipole")
            self.add_physics_process(proc, dipole_name, dipole_name)

class G1(FELSection):
    def __init__(self):
        sequence = MachineSequence(i1.cell)[: i1.start_ocelot.id]
        super().__init__(sequence)

    def build_physics_processes(self, _):
        message = "Tracking in G1 is unsupported and should be left to ASTRA."
        return [(ExceptionalProc(message), self.start, self.stop)]


class A1(FELSection):
    def __init__(self):
        start_sim = i1.start_ocelot
        acc1_stop = i1.a1_sim_stop
        lattice = MagneticLattice(i1.cell, start=start_sim, stop=acc1_stop)
        sequence = MachineSequence(lattice.sequence)
        super().__init__(sequence)

    def setup_physics(self):
        start_sim = i1.start_ocelot
        acc1_stop = i1.a1_sim_stop

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
        self.add_physics_process(smooth, start=start_sim.id, stop=start_sim.id)
        self.add_physics_process(sc, start=start_sim.id, stop=acc1_1_stop.id)
        self.add_physics_process(sc2, start=acc1_1_stop.id, stop=acc1_stop.id)
        self.add_physics_process(wake, start=i1.c_a1_1_1_i1.id, stop=acc1_stop.id)

        

class AH1(FELSection):
    def __init__(self):
        ah1_start = i1.ah1_sim_start
        acc39_stop = i1.stlat_47_i1
        lattice = MagneticLattice(i1.cell, start=ah1_start, stop=acc39_stop)
        super().__init__(lattice.sequence)

    def build_physics_processes(self, physics_config):
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
            (sc, self.start, self.stop),
            (wake, i1.c3_ah1_1_1_i1.id, self.stop),
        ]

        return physics_processes


class LH(FELSection):
    def __init__(self):
        acc39_stop = i1.stlat_47_i1
        lhm_stop = i1.dump_csr_start  # for going in I1D
        lattice = MagneticLattice(i1.cell + l1.cell, start=acc39_stop, stop=lhm_stop)
        super().__init__(lattice.sequence)

    def build_physics_processes(self, physics_config):
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
            (sc, self.start, self.stop),
            (csr, self.start, self.stop),
            (wake, self.start, self.stop),
        ]
        return physics_processes


class I1D(FELSection):
    def __init__(self):
        # init tracking lattice
        i1d_start = i1.dump_csr_start
        i1d_stop = i1d.otrc_64_i1d
        lattice = MagneticLattice(i1.cell + i1d.cell, start=i1d_start.id)
        super().__init__(lattice.sequence)

    def build_physics_processes(self, physics_config):
        # init tracking lattice
        # init physics processes
        sigma_min = physics_config.get_csr_sigma_min(self.start, self.stop)
        csr = make_csr(sigma_min=sigma_min, traj_step=0.0005, apply_step=0.005)
        sc = make_space_charge(step=5, nmesh_xyz=SCmesh, random_mesh=bRandomMesh)

        physics_processes = [
            (sc, self.start, self.stop),
            (csr, self.start, "BPMA.63.I1D"),
        ]
        return physics_processes



class DL(FELSection):
    def __init__(self):
        cell = MachineSequence(i1.cell + l1.cell)
        lh_stop_dl_start = i1.dump_csr_start.id  # "dogleg_start"
        dogleg_stop = "dogleg_stop_bc0_start"
        sequence = cell[lh_stop_dl_start:dogleg_stop]
        super().__init__(sequence)

    def build_physics_processes(self, physics_config):
        csr = make_csr(
            sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
            traj_step=0.0005,
            apply_step=0.005,
            n_bin=CSR_N_BIN,
        )
        wake = make_wake(
            WAKES_PATH / "mod_wake_0070.030_0073.450_MONO.dat",
            factor=1,
            w_sampling=WAKE_SAMPLING,
            filter_order=WAKE_FILTER_ORDER,
        )
        sc = make_space_charge(step=25, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)

        physics_processes = [
            (csr, self.start, self.stop),
            (sc, self.start, self.stop),
            (wake, self.stop, self.stop),
        ]
        return physics_processes

class BC0(FELSection):
    def __init__(self):
        l1_cell = MachineSequence(l1.cell)
        dogleg_stop_bc0_start = "dogleg_stop_bc0_start"
        bc0_stop_l1_start = "bc0_stop_l1_start"

        sequence = l1_cell[dogleg_stop_bc0_start:bc0_stop_l1_start]
        super().__init__(sequence)

    def build_physics_processes(self, physics_config):
        csr = make_csr(
            step=1,
            n_bin=CSR_N_BIN,
            sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
            traj_step=0.0005,
            apply_step=0.001,
        )
        sc = make_space_charge(step=40, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)
        physics_processes = [(sc, self.start, self.stop), (csr, self.start, self.stop)]
        return physics_processes



class L1(FELSection):
    def __init__(self):
        l1_cell = MachineSequence(l1.cell)
        # Just after end of BC0
        bc0_stop_l1_start = "bc0_stop_l1_start"
        # This is just before the start of BC1.
        a2_stop_bc1_start = "l1_stop_bc1_start"

        sequence = l1_cell[bc0_stop_l1_start:a2_stop_bc1_start]
        super().__init__(sequence)

    def build_physics_processes(self, physics_config):
        sc = make_space_charge(step=50, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)
        wake = make_wake(
            "wakes/RF/mod_TESLA_MODULE_WAKE_TAYLOR.dat",
            factor=4,
            step=100,
            w_sampling=WAKE_SAMPLING,
            filter_order=WAKE_FILTER_ORDER,
        )
        wake2 = make_wake(
            "wakes/mod_wake_0078.970_0159.280_MONO.dat",
            factor=1,
            w_sampling=WAKE_SAMPLING,
            filter_order=WAKE_FILTER_ORDER,
        )

        l1_first_cavity = "C.A2.1.1.L1"
        l1_last_cavity = "C.A2.4.8.L1"

        physics_processes = [
            (BEAM_SMOOTHER, self.start, self.start),
            (sc, self.start, self.stop),
            (wake, l1_first_cavity, l1_last_cavity),
            (wake2, self.stop, self.stop),
        ]

        return physics_processes


class BC1(FELSection):
    def __init__(self):
        l1_cell = MachineSequence(l1.cell)
        a2_stop_bc1_start = "l1_stop_bc1_start"
        bc1_stop_l2_start = "bc1_stop_l2_start"

        sequence = l1_cell[a2_stop_bc1_start:bc1_stop_l2_start]

        super().__init__(sequence)

    def build_physics_processes(self, physics_config):
        # init physics processes
        csr = make_csr(
            step=1,
            n_bin=CSR_N_BIN,
            sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
            traj_step=0.0005,
            apply_step=0.001,
        )
        sc = make_space_charge(step=40, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)

        return [(csr, self.start, self.stop), (sc, self.start, self.stop)]


class B1D(FELSection):
    def __init__(self):
        cell = MachineSequence(l1.cell + b1d.cell)
        bc1_stop_b1d_start = "STSEC.229.B1D"
        sequence = cell[bc1_stop_b1d_start:]
        super().__init__(sequence)

    def build_physics_processes(self, physics_config):
        csr = make_csr(
            step=1,
            n_bin=CSR_N_BIN,
            sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
            traj_step=0.0005,
            apply_step=0.001,
        )
        sc = make_space_charge(step=50, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)

        return [(csr, self.start, self.stop), (sc, self.start, self.stop)]



class L2(FELSection):
    def __init__(self):
        cell = MachineSequence(l1.cell + l2.cell)
        bc1_stop_l2_start = "bc1_stop_l2_start"
        l2_stop_bc2_start = "l2_stop_bc2_start"
        sequence = cell[bc1_stop_l2_start:l2_stop_bc2_start]
        super().__init__(sequence)

    def build_physics_processes(self, physics_config):
        sc = make_space_charge(step=100, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)
        wake = make_wake(
            WAKES_PATH / "RF/mod_TESLA_MODULE_WAKE_TAYLOR.dat", factor=4 * 3, step=200
        )
        wake2 = make_wake(WAKES_PATH / "mod_wake_0179.810_0370.840_MONO.dat", factor=1)

        first_cavity_l2 = "C.A3.1.1.L2"
        last_cavity_l2 = "C.A5.4.8.L2"

        return [
            (BEAM_SMOOTHER, self.start, self.start),
            (sc, self.start, self.stop),
            (wake, first_cavity_l2, last_cavity_l2),
            (wake2, self.stop, self.stop),
        ]


class BC2(FELSection):
    def __init__(self):
        l2_cell = MachineSequence(l2.cell)

        l2_stop_bc2_start = "l2_stop_bc2_start"
        bc2_stop_b2d_start = "bc2_stop_b2d_start"
        # init tracking lattice
        sequence = l2_cell[l2_stop_bc2_start:bc2_stop_b2d_start]
        super().__init__(sequence)

    def build_physics_processes(self, physics_config):
        csr = make_csr(
            step=1,
            n_bin=CSR_N_BIN,
            sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
            traj_step=0.0005,
            apply_step=0.001,
        )
        sc = make_space_charge(step=50, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)

        return [(csr, self.start, self.stop), (sc, self.start, self.stop)]


class B2D(FELSection):
    def __init__(self):
        cell = MachineSequence(l2.cell + b2d.cell)

        bc2_stop_b2d_start = "bc2_stop_b2d_start"
        sequence = cell[bc2_stop_b2d_start:]
        super().__init__(sequence)

    def build_physics_processes(self, physics_config):
        csr = make_csr(
            step=1,
            n_bin=CSR_N_BIN,
            sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
            traj_step=0.0005,
            apply_step=0.001,
        )
        sc = make_space_charge(step=50, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)

        return [(csr, self.start, self.stop), (sc, self.start, self.stop)]


class L3(FELSection):
    def __init__(self):
        bc2_stop = l2.tora_415_b2
        acc6t26_stop = l3.stop_l3.id

        # init tracking lattice
        cell = MachineSequence(l2.cell + l3.cell + cl.cell)
        sequence = cell[bc2_stop:acc6t26_stop]
        super().__init__(sequence)

    def build_physics_processes(self, physics_config):
        # init physics processes
        smooth = SmoothBeam()
        smooth.mslice = SmoothPar

        sc = SpaceCharge()
        sc.step = 5
        sc.nmesh_xyz = SCmesh
        sc.random_mesh = bRandomMesh

        wake = Wake()
        wake.wake_table = WakeTable(WAKES_PATH / "RF/mod_TESLA_MODULE_WAKE_TAYLOR.dat")
        wake.factor = 4 * 21
        wake.step = 10
        wake.w_sampling = 1000
        wake_add = Wake()
        wake_add.wake_table = WakeTable(
            WAKES_PATH / "mod_wake_0391.350_1629.700_MONO.dat"
        )
        wake_add.factor = 1

        app = PhaseSpaceAperture()
        app.taumin = -5
        app.taumax = 3

        return [
            (app, self.start, self.start),
            (smooth, self.start, self.start),
            (wake, l3.c_a6_1_1_l3.id, l3.c_a25_4_8_l3.id),
            (sc, self.start, self.stop),
            (wake_add, self.stop, self.stop),
        ]


class CL1(FELSection, SpontaneousRadiationSection):
    def __init__(self):
        acc6t26_stop = l3.stop_l3.id
        collimator1_stop = cl.bpma_1746_cl.id

        cell = MachineSequence(l3.cell + cl.cell)
        sequence = cell[acc6t26_stop:collimator1_stop]
        super().__init__(sequence)

    def build_physics_processes(self, physics_config):
        sc = make_space_charge(step=10, nmesh_xyz=SCmesh, random_mesh=bRandomMesh)
        csr = make_csr(
            traj_step=0.0005,
            apply_step=0.001,
            n_bin=CSRBin,
            sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
        )

        physics_processes = [(sc, self.start, self.stop), (csr, self.start, self.stop)]

        dipole_names = [
            "BE.1678.CL",
            "BL.1688.CL",
            "BL.1695.CL",
            "BE.1705.CL",
            "BE.1714.CL",
            "BL.1724.CL",
            "BL.1731.CL",
            "BE.1741.CL",
        ]
        if bISR:
            physics_processes.extend(
                self.make_dipole_spontaneous_radiation(dipole_names)
            )

        return physics_processes


class CL2(FELSection):
    def __init__(self):
        collimator1_stop = cl.bpma_1746_cl.id
        collimator2_stop = cl.bpma_1783_cl.id
        cell = MachineSequence(cl.cell)
        sequence = cell[collimator1_stop:collimator2_stop]
        super().__init__(sequence)

    def build_physics_processes(self, physics_config):
        sc = make_space_charge(step=1, nmesh_xyz=SCmesh, random_mesh=bRandomMesh)
        return [(sc, self.start, self.stop)]


class CL3(FELSection, SpontaneousRadiationSection):
    def __init__(self):
        collimator2_stop = cl.bpma_1783_cl.id
        collimator3_stop = cl.ensec_1854_cl.id
        # init tracking lattice
        cell = MachineSequence(cl.cell)
        sequence = cell[collimator2_stop:collimator3_stop]
        super().__init__(sequence)

    def build_physics_processes(self, physics_config):
        sc = make_space_charge(step=10, nmesh_xyz=SCmesh, random_mesh=bRandomMesh)
        csr = make_csr(
            traj_step=0.0005,
            apply_step=0.001,
            n_bin=CSRBin,
            sigma_min=physics_config.get_csr_sigma_min(self.start, self.stop),
        )
        wake_add = make_wake(
            WAKES_PATH / "mod_wake_1629.700_1831.200_MONO.dat", factor=1
        )

        processes = [
            (csr, self.start, self.stop),
            (sc, self.start, self.stop),
            (wake_add, self.stop, self.stop),
        ]
        if bISR:
            dipoles = [
                "BE.1786.CL",
                "BL.1796.CL",
                "BL.1803.CL",
                "BE.1813.CL",
                "BE.1822.CL",
                "BL.1832.CL",
                "BL.1839.CL",
                "BE.1849.CL",
            ]
            processes.extend(self.make_dipole_spontaneous_radiation(dipoles))

        return processes


class TLD(FELSection):
    def __init__(self):
        collimator3_stop = cl.ensec_1854_cl.id
        cell = MachineSequence(cl.cell + tld.cell)
        sequence = cell[collimator3_stop:]
        super().__init__(sequence)

    def build_physics_processes(self, physics_config):
        return []


class STN10(FELSection):
    def __init__(self):
        collimator3_stop = cl.ensec_1854_cl
        stN10_stop = sase1.ensec_2235_t2
        # stN10_stop = cl.ensub_1980_tl
        # init tracking lattice
        self.lattice = MagneticLattice(
            cl.cell + tl34_sase1.cell + sase1.cell,
            start=collimator3_stop,
            stop=stN10_stop,
            method=self.method,
        )
        # init physics processes
        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = SCmesh
        sc.random_mesh = bRandomMesh
        wake_add = Wake()
        wake_add.wake_table = WakeTable(
            "accelerator/wakes/mod_wake_1831.200_2035.190_MONO.dat"
        )
        wake_add.factor = 1
        wake_add1 = Wake()
        wake_add1.wake_table = WakeTable(
            "accelerator/wakes/mod_wake_2035.190_2213.000_MONO.dat"
        )
        wake_add1.factor = 1

        self.add_physics_process(sc, start=collimator3_stop, stop=stN10_stop)
        self.add_physics_process(
            wake_add, start=collimator3_stop, stop=collimator3_stop
        )
        self.add_physics_process(wake_add1, start=stN10_stop, stop=stN10_stop)

# class SASE1(FELSection):
#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)
#         # setting parameters
#         self.lattice_name = 'SASE1'
#         self.unit_step = 5

#         self.input_beam_file = self.particle_dir + 'section_STN10.npz'
#         self.output_beam_file = self.particle_dir + 'section_SASE1.npz'
#         self.tws_file = self.tws_dir + "tws_section_SASE1.npz"
#         # last element sase1 - stsec_2461_t4
#         stN10_stop = sase1.ensec_2235_t2
#         sase1_stop = sase1.stsec_2461_t4
#         # init tracking lattice
#         self.lattice = MagneticLattice(sase1.cell, start=stN10_stop, stop=sase1_stop, method=self.method)

#         # init physics processes
#         wake = Wake()
#         wake.wake_table = WakeTable('accelerator/wakes/Undulator/wake_undulator_OCELOT.txt')
#         wake.step = 10
#         wake.w_sampling = WakeSampling
#         wake.factor = 35 * 6.1
#         sc = SpaceCharge()
#         sc.step = 1
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh

#         sre = SpontanRadEffects()
#         sre.K = 3.9
#         sre.lperiod = 0.04
#         sre.filling_coeff = 5 / 6.1

#         self.add_physics_process(wake, start=sase1.match_2247_sa1,stop=sase1_stop)
#         self.add_physics_process(sc, start=stN10_stop, stop=sase1_stop)
#         self.add_physics_process(sre, start=stN10_stop, stop=sase1_stop)


# class T4(FELSection):
#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)
#         # setting parameters
#         self.lattice_name = 'T4'
#         self.unit_step = 0.2

#         self.input_beam_file = self.particle_dir + 'section_SASE1.npz'
#         self.output_beam_file = self.particle_dir + 'section_T4.npz'
#         self.tws_file = self.tws_dir + "tws_section_T4.npz"
#         # last element sase1 - stsec_2461_t4
#         sase1_stop = sase1.stsec_2461_t4
#         t4_stop = t4.ensub_2800_t4
#         csr_start = t4.t4_start_csr
#         csr_stop = t4.bpma_2606_t4
#         # init tracking lattice
#         self.lattice = MagneticLattice(sase1.cell + t4.cell, start=sase1_stop, stop=t4_stop, method=self.method)

#         # init physics processes

#         sc = SpaceCharge()
#         sc.step = 25
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh

#         csr = CSR()
#         csr.traj_step = 0.0005
#         csr.apply_step = 0.005
#         csr.n_bin = CSRBin
#         csr.sigma_min = Sig_Z[3]*CSRSigmaFactor

#         sc2 = SpaceCharge()
#         sc2.step = 25
#         sc2.nmesh_xyz = [31, 31, 31]
#         sc2.random_mesh = bRandomMesh

#         sc_in_bend = SpaceCharge()
#         sc_in_bend.step = 25
#         sc_in_bend.nmesh_xyz = SCmesh
#         sc_in_bend.random_mesh = bRandomMesh

#         # creation of wake object with parameters
#         wake = Wake()
#         wake.wake_table = WakeTable('accelerator/wakes/Dechirper/wake_hor_axis_500um.txt')
#         wake.w_sampling = 500
#         wake.factor = 1
#         wake.step = 1  # step in Navigator.unit_step, dz = Navigator.unit_step * wake.step [m]


#         # creation of wake object with parameters
#         wake_vert = Wake()
#         wake_vert.factor = 1
#         #wake_vert.wake_table = WakeTable('accelerator/wakes/wake_vert_1m_500mkm.txt')
#         wake_vert.wake_table = WakeTable('accelerator/wakes/Dechirper/wake_vert_axis_500um.txt')
#         wake_vert.w_sampling = 500
#         wake_vert.step = 1  # step in Navigator.unit_step, dz = Navigator.unit_step * wake.step [m]

#         svb4 = SaveBeam(filename=self.particle_dir + "before_structure.npz")
#         svb3 = SaveBeam(filename=self.particle_dir + "after_structure.npz")
#         svb1 = SaveBeam(filename=self.particle_dir + "screen1.npz")
#         svb2 = SaveBeam(filename=self.particle_dir + "screen2.npz")

#         #self.add_physics_process(svb4, start=t4.wake_start, stop=t4.wake_start)
#         #self.add_physics_process(svb3, start=t4.wake_stop, stop=t4.wake_stop)
#         #self.add_physics_process(svb1, start=t4.m_img1, stop=t4.m_img1)
#         #self.add_physics_process(svb2, start=t4.m_img2, stop=t4.m_img2)
#         #self.add_physics_process(sc, start=sase1_stop, stop=csr_start)
#         self.add_physics_process(csr, start=csr_start, stop=csr_stop)
#         #self.add_physics_process(sc2, start=csr_stop, stop=t4.ensub_2800_t4)
#         #self.add_physics_process(wake, start=t4.wake_start, stop=t4.m_tds)
#         #self.add_physics_process(wake_vert, start=t4.m_tds, stop=t4.wake_stop)
#         #self.add_physics_process(sc_in_bend, start=csr_start, stop=csr_stop)


# class SASE3(FELSection):
#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)
#         # setting parameters
#         self.lattice_name = 'SASE3'
#         self.unit_step = 1

#         self.input_beam_file = self.particle_dir + 'section_T4.npz'
#         self.output_beam_file = self.particle_dir + 'section_SASE3.npz'
#         self.tws_file = self.tws_dir + "tws_section_SASE3.npz"

#         start = sase3.ensec_2800_t4
#         stop = sase3.ensec_2940_sa3
#         # init tracking lattice
#         self.lattice = MagneticLattice(sase3.cell, start=start, stop=stop, method=self.method)

#         # init physics processes

#         sc = SpaceCharge()
#         sc.step = 10
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh

#         self.add_physics_process(sc, start=start, stop=stop)


# class T4D(FELSection):
#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)
#         # setting parameters
#         self.lattice_name = 'SASE3'
#         self.unit_step = 1

#         self.input_beam_file = self.particle_dir + 'section_SASE3.npz'
#         self.output_beam_file = self.particle_dir + 'section_T4D.npz'
#         self.tws_file = self.tws_dir + "tws_section_tT4D.npz"

#         start = t4d.stsec_2940_t4d
#         stop = t4d.ensec_3106_t4d
#         # init tracking lattice
#         self.lattice = MagneticLattice(t4d.cell, start=start, stop=stop, method=self.method)

#         # init physics processes

#         sc = SpaceCharge()
#         sc.step = 10
#         sc.nmesh_xyz = SCmesh
#         self.add_physics_process(sc, start=start, stop=stop)
#         sc.random_mesh = bRandomMesh

#         csr = CSR()
#         csr.traj_step = 0.0005
#         csr.apply_step = 0.005
#         csr.n_bin = CSRBin
#         csr.sigma_min = Sig_Z[3]*CSRSigmaFactor

#         self.add_physics_process(csr, start=t4d.tora_3065_t4d, stop=t4d.qk_3090_t4d)


# class T4_short(FELSection):
#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)
#         # setting parameters
#         self.lattice_name = 'T4'
#         self.unit_step = 1
#         self.calc_tws = False

#         self.input_beam_file = self.particle_dir + 'before_structure.npz'
#         self.output_beam_file = self.particle_dir + 'section_T4.npz'
#         self.tws_file = self.tws_dir + "tws_section_T4.npz"
#         # last element sase1 - stsec_2461_t4
#         sase1_stop = t4.wake_start# sase1.stsec_2461_t4
#         t4_stop = t4.m_img1 #t4.ensub_2800_t4
#         csr_start = t4.t4_start_csr
#         csr_stop = t4.bpma_2606_t4
#         # init tracking lattice
#         self.lattice = MagneticLattice(sase1.cell + t4.cell, start=sase1_stop, stop=t4_stop, method=self.method)

#         # init physics processes

#         sc = SpaceCharge()
#         sc.step = 25
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh

#         csr = CSR()
#         csr.traj_step = 0.0005
#         csr.apply_step = 0.005
#         csr.n_bin = CSRBin
#         csr.sigma_min = Sig_Z[3]*CSRSigmaFactor

#         sc2 = SpaceCharge()
#         sc2.step = 25
#         sc2.nmesh_xyz = SCmesh


#         # creation of wake object with parameters
#         wake = Wake()
#         wake.wake_table = WakeTable('accelerator/wakes/tt.txt')

#         # w_sampling - defines the number of the equidistant sampling points for the one-dimensional
#         # wake coefficients in the Taylor expansion of the 3D wake function.
#         wake.w_sampling = 500
#         wake.factor = 1
#         wake.step = 1  # step in Navigator.unit_step, dz = Navigator.unit_step * wake.step [m]


#         # creation of wake object with parameters
#         wake_vert = Wake()
#         wake_vert.factor = 1
#         wake_vert.wake_table = WakeTable('accelerator/wakes/tt.txt')
#         wake_vert.w_sampling = 500
#         wake_vert.step = 1  # step in Navigator.unit_step, dz = Navigator.unit_step * wake.step [m]


#         #svb1 = SaveBeam(filename=self.particle_dir + "screen1.npz")

#         #self.add_physics_process(svb1, start=t4.m_img1, stop=t4.m_img1)

#         #svb2 = SaveBeam(filename=self.particle_dir + "screen2.npz")
#         svb3 = SaveBeam(filename=self.particle_dir + "after_structure.npz")
#         #svb4 = SaveBeam(filename=self.particle_dir + "before_structure.npz")
#         #self.add_physics_process(svb2, start=t4.m_img2, stop=t4.m_img2)

#         self.add_physics_process(wake_vert, start=t4.wake_start, stop=t4.m_tds)
#         #self.add_physics_process(wake, start=t4.m_tds, stop=t4.wake_stop)

#         #self.add_physics_process(svb3, start=t4.wake_stop, stop=t4.wake_stop)
#         #self.add_physics_process(svb4, start=t4.wake_start, stop=t4.wake_start)

#         #self.add_physics_process(sc, start=sase1_stop, stop=csr_start)
#         #self.add_physics_process(sc, start=sase1_stop, stop=csr_start)
#         #self.add_physics_process(csr, start=csr_start, stop=csr_stop)
#         #self.add_physics_process(sc2, start=csr_stop, stop=t4.ensub_2800_t4)

#         sc_in_bend = SpaceCharge()
#         sc_in_bend.step = 25
#         sc_in_bend.nmesh_xyz = SCmesh
#         #self.add_physics_process(sc_in_bend, start=csr_start, stop=csr_stop)

# ################  SASE2 branch ####################################################

# class T1(FELSection):
#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)

#         # setting parameters
#         self.lattice_name = 'T1'
#         self.unit_step = 0.2

#         self.input_beam_file = self.particle_dir + 'section_CL3.npz'
#         self.output_beam_file = self.particle_dir + 'section_T1.npz'
#         self.tws_file = self.tws_dir + "tws_section_T1.npz"

#         collimator3_stop = cl.ensec_1854_cl
#         t1_stop = t1.stsec_2197_sa2
#         # init tracking lattice
#         self.lattice = MagneticLattice(cl.cell + tl34.cell + t1.cell, start=collimator3_stop, stop=t1_stop, method=self.method)

#         # init physics processes
#         csr = CSR()
#         csr.traj_step = 0.0005
#         csr.apply_step = 0.001
#         csr.n_bin = CSRBin
#         csr.sigma_min = Sig_Z[3] * CSRSigmaFactor
#         #csr.rk_traj=True
#         #csr.energy = 17.5

#         csr2 = CSR()
#         csr2.traj_step = 0.0005
#         csr2.apply_step = 0.001
#         csr2.n_bin = CSRBin
#         csr2.sigma_min = Sig_Z[3] *CSRSigmaFactor

#         csr3 = CSR()
#         csr3.traj_step = 0.0005
#         csr3.apply_step = 0.001
#         csr3.n_bin = CSRBin
#         csr3.sigma_min = Sig_Z[3] *CSRSigmaFactor

#         wake_add = Wake()
#         wake_add.wake_table = WakeTable('accelerator/wakes/mod_wake_1831.200_2035.190_MONO.dat')
#         wake_add.factor = 1
#         wake_add1 = Wake()
#         wake_add1.wake_table = WakeTable('accelerator/wakes/mod_wake_2035.190_2213.000_MONO.dat')
#         wake_add1.factor = 1

#         self.add_physics_process(wake_add, start=collimator3_stop, stop=collimator3_stop)
#         self.add_physics_process(wake_add1, start=t1_stop, stop=t1_stop)
#         #self.add_physics_process(csr, start=tl34.stsub_1997_tl, stop=t1.bend_stop_sa2)
#         self.add_physics_process(csr, start=tl34.stsub_1997_tl, stop=t1.before_XYQuad_sa2)
#         self.add_physics_process(csr2, start=t1.after_XYQuad_sa2, stop=t1.bend1_stop_sa2)
#         self.add_physics_process(csr3, start=t1.bend1_stop_sa2, stop=t1.bend_stop_sa2)

#         if bISR:
#             LD = tl34.kl_1998_tl.l;  teta = tl34.kl_1998_tl.angle; ro = LD / np.sin(teta);
#             sre1 = SpontanRadEffects(); sre1.radius = ro;  sre1.type = 'dipole';
#             self.add_physics_process(sre1, tl34.M1kl_1998_tl, tl34.M2kl_1998_tl)
#             LD = tl34.kl_1999_tl.l;  teta = tl34.kl_1999_tl.angle; ro = LD / np.sin(teta);
#             sre2 = SpontanRadEffects(); sre2.radius = ro;  sre2.type = 'dipole';
#             self.add_physics_process(sre2, tl34.M1kl_1999_tl, tl34.M2kl_1999_tl)
#             LD = tl34.kl_2000_tl.l;  teta = tl34.kl_2000_tl.angle; ro = LD / np.sin(teta);
#             sre3 = SpontanRadEffects(); sre3.radius = ro;  sre3.type = 'dipole';
#             self.add_physics_process(sre3, tl34.M1kl_2000_tl, tl34.M2kl_2000_tl)
#             LD = tl34.kl_2001_tl.l;  teta = tl34.kl_2001_tl.angle; ro = LD / np.sin(teta);
#             sre4 = SpontanRadEffects(); sre4.radius = ro;  sre4.type = 'dipole';
#             self.add_physics_process(sre4, tl34.M1kl_2001_tl, tl34.M2kl_2001_tl)
#             LD = tl34.kl_2002_tl.l;  teta = tl34.kl_2002_tl.angle; ro = LD / np.sin(teta);
#             sre5 = SpontanRadEffects(); sre5.radius = ro;  sre5.type = 'dipole';
#             self.add_physics_process(sre5, tl34.M1kl_1998_tl, tl34.M2kl_1998_tl)
#             LD = tl34.qf_2012_tl.l;  teta = tl34.qf_2012_tl.angle; ro = LD / np.sin(teta);
#             sre6 = SpontanRadEffects(); sre6.radius = ro;  sre6.type = 'dipole';
#             self.add_physics_process(sre6, tl34.M1qf_2012_tl, tl34.M2qf_2012_tl)
#             LD = t1.bz_2025_t1.l;  teta = t1.bz_2025_t1.angle; ro = LD / np.sin(teta);
#             sre7 = SpontanRadEffects(); sre7.radius = ro;  sre7.type = 'dipole';
#             self.add_physics_process(sre7, t1.M1bz_2025_t1, t1.M2bz_2025_t1)
#             LD = t1.bz_2030_t1.l;  teta = t1.bz_2030_t1.angle; ro = LD / np.sin(teta);
#             sre8 = SpontanRadEffects(); sre8.radius = ro;  sre8.type = 'dipole';
#             self.add_physics_process(sre8, t1.M1bz_2030_t1, t1.M2bz_2030_t1)
#             LD = t1.bz_2031_t1.l;  teta = t1.bz_2031_t1.angle; ro = LD / np.sin(teta);
#             sre9 = SpontanRadEffects(); sre9.radius = ro;  sre9.type = 'dipole';
#             self.add_physics_process(sre9, t1.M1bz_2031_t1, t1.M2bz_2031_t1)
#             LD = t1.bz_2033_t1.l;  teta = t1.bz_2033_t1.angle; ro = LD / np.sin(teta);
#             sre10 = SpontanRadEffects(); sre10.radius = ro;  sre10.type = 'dipole';
#             self.add_physics_process(sre10, t1.M1bz_2033_t1, t1.M2bz_2033_t1)
#             LD = t1.bd_2050_t1.l;  teta = t1.bd_2050_t1.angle; ro = LD / np.sin(teta);
#             sre11 = SpontanRadEffects(); sre11.radius = ro;  sre11.type = 'dipole';
#             self.add_physics_process(sre11, t1.M1bz_2050_t1, t1.M2bz_2050_t1)
#             LD = t1.bd_2062_t1.l;  teta = t1.bd_2062_t1.angle; ro = LD / np.sin(teta);
#             sre12 = SpontanRadEffects(); sre12.radius = ro;  sre12.type = 'dipole';
#             self.add_physics_process(sre12, t1.M1bz_2062_t1, t1.M2bz_2062_t1)
#             LD = t1.bd_2077_t1.l;  teta = t1.bd_2077_t1.angle; ro = LD / np.sin(teta);
#             sre13 = SpontanRadEffects(); sre13.radius = ro;  sre13.type = 'dipole';
#             self.add_physics_process(sre13, t1.M1bz_2077_t1, t1.M2bz_2077_t1)
#             LD = t1.bd_2079_t1.l;  teta = t1.bd_2079_t1.angle; ro = LD / np.sin(teta);
#             sre14 = SpontanRadEffects(); sre14.radius = ro;  sre14.type = 'dipole';
#             self.add_physics_process(sre14, t1.M1bz_2079_t1, t1.M2bz_2079_t1)
#             LD = t1.bd_2080_t1.l;  teta = t1.bd_2080_t1.angle; ro = LD / np.sin(teta);
#             sre15 = SpontanRadEffects(); sre15.radius = ro;  sre15.type = 'dipole';
#             self.add_physics_process(sre15, t1.M1bz_2080_t1, t1.M2bz_2080_t1)
#             LD = t1.bd_2082_t1.l;  teta = t1.bd_2082_t1.angle; ro = LD / np.sin(teta);
#             sre16 = SpontanRadEffects(); sre16.radius = ro;  sre16.type = 'dipole';
#             self.add_physics_process(sre16, t1.M1bz_2082_t1, t1.M2bz_2082_t1)
#             LD = t1.bd_2084_t1.l;  teta = t1.bd_2084_t1.angle; ro = LD / np.sin(teta);
#             sre17 = SpontanRadEffects(); sre17.radius = ro;  sre17.type = 'dipole';
#             self.add_physics_process(sre17, t1.M1bz_2084_t1, t1.M2bz_2084_t1)
#             LD = t1.bd_2097_t1.l;  teta = t1.bd_2097_t1.angle; ro = LD / np.sin(teta);
#             sre18 = SpontanRadEffects(); sre18.radius = ro;  sre18.type = 'dipole';
#             self.add_physics_process(sre18, t1.M1bz_2097_t1, t1.M2bz_2097_t1)

# class SASE2(FELSection):
#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)
#         # setting parameters
#         self.lattice_name = 'SASE2'
#         self.unit_step = 2
#         self.input_beam_file = self.particle_dir + 'section_T1.npz'
#         self.output_beam_file = self.particle_dir + 'section_SASE2.npz'
#         self.tws_file = self.tws_dir + "tws_section_SASE2.npz"

#         sase2_start = sase2.stsec_2197_sa2
#         sase2_stop=sase2.stsec_2423_t3

#         self.lattice = MagneticLattice(sase2.cell, start=sase2_start, stop=sase2_stop,
#                                            method=self.method)
#         # init physics processes
#         wake = Wake()
#         wake.wake_table = WakeTable('accelerator/wakes/Undulator/wake_undulator_OCELOT.txt')
#         wake.step = 10
#         wake.w_sampling = 1000
#         wake.factor = 35*6.1
#         sc = SpaceCharge()
#         sc.step = 1
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh

#         sre=SpontanRadEffects()
#         sre.K = 3.9
#         sre.lperiod = 0.04
#         sre.filling_coeff=5/6.1

#         self.add_physics_process(wake, start=sase2_stop, stop=sase2_stop)
#         self.add_physics_process(sc, start=sase2_start, stop=sase2_stop)
#         #self.add_physics_process(sre, start=sase2_start, stop=sase2_stop)

# class T3(FELSection):
#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)
#         # setting parameters
#         self.lattice_name = 'T3'
#         self.unit_step = 10
#         self.input_beam_file = self.particle_dir + 'section_SASE2.npz'
#         self.output_beam_file = self.particle_dir + 'section_T3.npz'
#         self.tws_file = self.tws_dir + "tws_section_T3.npz"

#         self.lattice = MagneticLattice(t3.cell, start=t3.stsec_2423_t3, stop=t3.stsec_2743_t5,
#                                            method=self.method)

#         csr = CSR()
#         csr.traj_step = 0.0005
#         csr.apply_step = 0.001
#         csr.n_bin = CSRBin
#         csr.sigma_min = Sig_Z[3] * CSRSigmaFactor
#         self.add_physics_process(csr, t3.M1be_2546_t3, t3.endCSR_t3)
#         if bISR:
#             LD = t3.be_2546_t3.l;  teta = t3.be_2546_t3.angle; ro = LD / np.sin(teta);
#             sre1 = SpontanRadEffects(); sre1.radius = ro;  sre1.type = 'dipole';
#             self.add_physics_process(sre1, t3.M1be_2546_t3, t3.M2be_2546_t3)
#             LD = t3.be_2566_t3.l;  teta = t3.be_2566_t3.angle; ro = LD / np.sin(teta);
#             sre1 = SpontanRadEffects(); sre1.radius = ro;  sre1.type = 'dipole';
#             self.add_physics_process(sre1, t3.M1be_2566_t3, t3.M2be_2566_t3)

# class T5(FELSection):
#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)
#         # setting parameters
#         self.lattice_name = 'T5'
#         self.unit_step = 10
#         self.input_beam_file = self.particle_dir + 'section_T3.npz'
#         self.output_beam_file = self.particle_dir + 'section_T5.npz'
#         self.tws_file = self.tws_dir + "tws_section_T5.npz"

#         self.lattice = MagneticLattice(t5.cell, start=t5.stsec_2743_t5, stop=t5.stsec_3039_t5d,
#                                            method=self.method)
