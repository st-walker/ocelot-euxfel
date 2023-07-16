from pathlib import Path
from importlib_resources import files

from ocelot.utils.section_track import *
from .lattice import i1
from .lattice import i1d
from .lattice import l1
from .lattice import l2
from .lattice import b2d


from ocelot.cpbd.physics_proc import *

from oxfel.fel_track import FELSection, MachineSequence

# Sig_Z=(0.0019996320155001497, 0.0006893836215002082, 0.0001020391309281775, 1.25044082708419e-05) #500pC 5kA
# Sig_Z=(0.0019996320155001497, 0.0006817907866411071, 9.947650872824487e-05, 7.13045869665955e-06)  #500pC 10kA
# Sig_Z=(0.0018761888067590127, 0.0006359220169656093, 9.204477386791353e-05, 7.032551498646372e-06) #250pC 5kA
# Sig_Z=(0.0018856911379360524, 0.0005463919476045524, 6.826162032352288e-05, 1.0806534547678727e-05) #100pC 1kA
Sig_Z = (
    0.0018732376720197858,
    0.000545866016784069,
    7.09234589639138e-05,
    2.440742745010469e-06,
)  # 100 pC 5kA
# Sig_Z=(0.0013314283765668853, 0.0004502566926198658, 4.64037216210807e-05, 2.346018397815618e-06) #100 pC 5kA SC
# Sig_Z=(0.0013314187263949542, 0.00045069372029991764, 4.537451914820527e-05, 4.0554988027793585e-06)#100 pC 2.5kA SC


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

# 250pC, no compression.
Sig_Z = 0.0013
# Sig_Z = Sig_Z


def make_space_charge(*, step, nmesh_xyz=None, random_mesh=None):
    sc = SpaceCharge()
    sc.step = step
    if nmesh_xyz is not None:
        sc.nmesh_xyz = nmesh_xyz
    if random_mesh is not None:
        sc.random_mesh = random_mesh
    return sc


def make_wake(rel_path, *, factor, step=None, w_sampling=None, filter_order=None):
    wake = Wake()
    wake.wake_table = WakeTable(Path(__file__).parent / rel_path)
    wake.factor = factor
    if step is not None:
        wake.step = step

    if w_sampling is not None:
        wake.w_sampling = w_sampling
    if filter_order is not None:
        wake.filter_order = filter_order

    return wake


def make_beam_transform(*, betx, bety, alfx, alfy, tr_slice):
    tws = Twiss()
    tws.beta_x = betx
    tws.beta_y = bety
    tws.alpha_x = alfx
    tws.alpha_y = alfy
    tws.gamma_x = (1 + alfx**2) / betx
    tr = BeamTransform(tws=tws)
    tr.slice = tr_slice
    return tr


def make_laser_modulator():
    lh = LaserModulator()
    lh.dE = 300
    lh.sigma_l = 300
    lh.sigma_x = 300e-6
    lh.sigma_y = 300e-6
    lh.z_waist = None
    return lh


def make_csr(*, sigma_min, traj_step, apply_step, n_bin=None, step=None):
    csr = CSR()
    csr.sigma_min = sigma_min
    csr.traj_step = traj_step
    csr.apply_step = apply_step
    if n_bin is not None:
        csr.n_bin = n_bin
    if step is not None:
        csr.step = step
    return csr


class G1(FELSection):
    def __init__(self):
        lattice = MagneticLattice(i1.cell, stop=i1.stop_astra)
        sequence = MachineSequence(lattice.sequence)
        super().__init__(sequence)


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

    def setup_physics(self):
        ah1_start = i1.ah1_sim_start
        acc39_stop = i1.stlat_47_i1
        # init physics processes
        sc = SpaceCharge()
        sc.step = 5
        sc.nmesh_xyz = SCmesh
        sc.random_mesh = bRandomMesh
        wake = Wake()
        wake.wake_table = WakeTable(WAKES_PATH / "RF/wake_table_AH1.dat")
        wake.factor = 1
        wake.step = 10
        wake.w_sampling = WakeSampling
        wake.filter_order = WakeFilterOrder
        # adding physics processes
        self.add_physics_process(sc, start=ah1_start.id, stop=acc39_stop.id)
        self.add_physics_process(wake, start=i1.c3_ah1_1_1_i1.id, stop=acc39_stop.id)


class LH(FELSection):
    def __init__(self):
        acc39_stop = i1.stlat_47_i1
        lhm_stop = i1.dump_csr_start  # for going in I1D
        lattice = MagneticLattice(i1.cell + l1.cell, start=acc39_stop, stop=lhm_stop)
        super().__init__(lattice.sequence)

    def setup_physics(self):
        acc39_stop = i1.stlat_47_i1
        lhm_stop = i1.dump_csr_start  # for going in I1D
        csr = CSR()
        # csr.sigma_min = Sig_Z[0] * CSRSigmaFactor
        csr.sigma_min = Sig_Z * CSRSigmaFactor
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

        self.add_physics_process(sc, start=acc39_stop.id, stop=lhm_stop.id)
        self.add_physics_process(csr, start=acc39_stop.id, stop=lhm_stop.id)
        self.add_physics_process(wake, start=acc39_stop.id, stop=lhm_stop.id)


class I1D_Screen(FELSection):
    def __init__(self):
        # init tracking lattice
        i1d_start = i1.dump_csr_start
        i1d_stop = i1d.otrc_64_i1d
        lattice = MagneticLattice(i1.cell + i1d.cell, start=i1d_start.id)
        super().__init__(lattice.sequence)

    def setup_physics(self):
        i1d_start = i1.dump_csr_start
        i1d_stop = i1d.otrc_64_i1d
        # init tracking lattice
        # init physics processes
        sigma = Sig_Z
        csr = CSR()
        csr.sigma_min = sigma * 0.1
        csr.traj_step = 0.0005
        csr.apply_step = 0.005

        sc = SpaceCharge()
        sc.step = 5
        sc.nmesh_xyz = SCmesh
        sc.random_mesh = bRandomMesh

        self.add_physics_process(sc, start=i1d_start.id, stop=i1d_stop.id)
        self.add_physics_process(csr, start=i1d_start.id, stop=i1d.bpma_63_i1d.id)
        # self.add_physics_process(sv_dump, start=i1d_start.id, stop=i1d_start.id)


class DL(FELSection):
    def __init__(self):
        cell = MachineSequence(i1.cell + l1.cell)
        lh_stop_dl_start = i1.dump_csr_start.id  # "dogleg_start"
        dogleg_stop = "dogleg_stop_bc0_start"
        sequence = cell[lh_stop_dl_start:dogleg_stop]
        super().__init__(sequence)

    def setup_physics(self):
        lh_stop_dl_start = i1.dump_csr_start.id  # "dogleg_start"
        dogleg_stop = "dogleg_stop_bc0_start"

        csr = make_csr(
            sigma_min=Sig_Z * CSR_SIGMA_FACTOR,
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

        self.add_physics_process(csr, start=lh_stop_dl_start, stop=dogleg_stop)
        self.add_physics_process(sc, start=lh_stop_dl_start, stop=dogleg_stop)
        self.add_physics_process(wake, start=dogleg_stop, stop=dogleg_stop)


class BC0(FELSection):
    def __init__(self):
        l1_cell = MachineSequence(l1.cell)
        dogleg_stop_bc0_start = "dogleg_stop_bc0_start"
        bc0_stop_l1_start = "bc0_stop_l1_start"

        sequence = l1_cell[dogleg_stop_bc0_start:bc0_stop_l1_start]
        super().__init__(sequence)

    def setup_physics(self):
        dogleg_stop_bc0_start = "dogleg_stop_bc0_start"
        bc0_stop_l1_start = "bc0_stop_l1_start"

        csr = make_csr(
            step=1,
            n_bin=CSR_N_BIN,
            sigma_min=Sig_Z * CSR_SIGMA_FACTOR,
            traj_step=0.0005,
            apply_step=0.001,
        )
        sc = make_space_charge(step=40, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)

        self.add_physics_process(
            sc, start=dogleg_stop_bc0_start, stop=bc0_stop_l1_start
        )
        self.add_physics_process(
            csr, start=dogleg_stop_bc0_start, stop=bc0_stop_l1_start
        )


# class BC0(FELSection):

#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)

#         # setting parameters
#         self.lattice_name = 'BC0'
#         self.unit_step = 0.02


#         self.input_beam_file = self.particle_dir / 'section_DL.npz'
#         self.output_beam_file = self.particle_dir / 'section_BC0.npz'
#         self.tws_file = self.tws_dir / "tws_section_BC0.npz"
#         # init tracking lattice
#         st4_stop = l1.stlat_96_i1
#         bc0_stop = l1.enlat_101_i1
#         self.lattice = MagneticLattice(l1.cell, start=st4_stop, stop=bc0_stop, method=self.method)

#         # init physics processes
#         csr = CSR()
#         csr.step = 1
#         csr.n_bin = CSRBin
#         csr.sigma_min = Sig_Z[1]*CSRSigmaFactor
#         csr.traj_step = 0.0005
#         csr.apply_step = 0.001

#         sc = SpaceCharge()
#         sc.step = 40
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh
#         match_bc0 = st4_stop
#         self.add_physics_process(sc, start=match_bc0, stop=bc0_stop)
#         self.add_physics_process(csr, start=match_bc0, stop=bc0_stop)
#         self.dipoles = [l1.bb_96_i1, l1.bb_98_i1, l1.bb_100_i1, l1.bb_101_i1]
#         self.dipole_len = 0.5
#         self.bc_gap=1.0


class L1(FELSection):
    def __init__(self):
        l1_cell = MachineSequence(l1.cell)
        # Just after end of BC0
        bc0_stop_l1_start = "bc0_stop_l1_start"
        # This is just before the start of BC1.
        a2_stop_bc1_start = "l1_stop_bc1_start"

        sequence = l1_cell[bc0_stop_l1_start:a2_stop_bc1_start]
        super().__init__(sequence)

    def setup_physics(self):
        # Just after end of BC0
        bc0_stop_l1_start = "bc0_stop_l1_start"
        # This is just before the start of BC1.
        a2_stop_bc1_start = "l1_stop_bc1_start"

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
        self.add_physics_process(
            BEAM_SMOOTHER, start=bc0_stop_l1_start, stop=bc0_stop_l1_start
        )
        self.add_physics_process(sc, start=bc0_stop_l1_start, stop=a2_stop_bc1_start)
        self.add_physics_process(wake, start=l1_first_cavity, stop=l1_last_cavity)
        self.add_physics_process(wake2, start=a2_stop_bc1_start, stop=a2_stop_bc1_start)


# class L1(FELSection):

#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)

#         # setting parameters
#         self.lattice_name = 'L1'
#         self.unit_step = 0.02


#         self.input_beam_file = self.particle_dir / 'section_BC0.npz'
#         self.output_beam_file = self.particle_dir / 'section_L1.npz'
#         self.tws_file = self.tws_dir / "tws_section_L1.npz"
#         bc0_stop = l1.enlat_101_i1
#         acc2_stop = l1.stlat_182_b1
#         # init tracking lattice
#         self.lattice = MagneticLattice(l1.cell, start=bc0_stop, stop=acc2_stop, method=self.method)

#         # init physics processes
#         smooth = SmoothBeam()
#         smooth.mslice = SmoothPar

#         sc = SpaceCharge()
#         sc.step = 50
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh
#         wake = Wake()
#         wake.wake_table = WakeTable(WAKES_PATH / 'mod_TESLA_MODULE_WAKE_TAYLOR.dat')
#         wake.factor = 4
#         wake.step = 100
#         wake.w_sampling = WakeSampling
#         wake.filter_order = WakeFilterOrder
#         wake_add = Wake()
#         wake_add.wake_table = WakeTable(WAKES_PATH / 'mod_wake_0078.970_0159.280_MONO.dat')
#         wake_add.factor = 1
#         wake_add.w_sampling = WakeSampling
#         wake_add.filter_order = WakeFilterOrder
#         match_acc2 = bc0_stop
#         L1_wake_kick = acc2_stop
#         self.add_physics_process(smooth, start=match_acc2.id, stop=match_acc2.id)
#         self.add_physics_process(sc, start=match_acc2.id, stop=L1_wake_kick.id)
#         self.add_physics_process(wake, start=l1.c_a2_1_1_l1.id, stop=l1.c_a2_4_8_l1.id)
#         self.add_physics_process(wake_add, start=L1_wake_kick.id, stop=L1_wake_kick.id)


class BC1(FELSection):
    def __init__(self):
        l1_cell = MachineSequence(l1.cell)
        a2_stop_bc1_start = "l1_stop_bc1_start"
        bc1_stop_l2_start = "bc1_stop_l2_start"

        sequence = l1_cell[a2_stop_bc1_start:bc1_stop_l2_start]

        super().__init__(sequence)

    def setup_physics(self):
        a2_stop_bc1_start = "l1_stop_bc1_start"
        bc1_stop_l2_start = "bc1_stop_l2_start"

        # init physics processes
        csr = make_csr(
            step=1,
            n_bin=CSR_N_BIN,
            sigma_min=Sig_Z * CSR_SIGMA_FACTOR,
            traj_step=0.0005,
            apply_step=0.001,
        )
        sc = make_space_charge(step=40, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)

        self.add_physics_process(csr, start=a2_stop_bc1_start, stop=bc1_stop_l2_start)
        self.add_physics_process(sc, start=a2_stop_bc1_start, stop=bc1_stop_l2_start)

        # self.dipoles = [l1.bb_182_b1, l1.bb_191_b1, l1.bb_193_b1, l1.bb_202_b1]
        # self.dipole_len = 0.5
        # self.bc_gap=8.5


# class BC1(FELSection):

#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)

#         # setting parameters
#         self.lattice_name = 'BC1'
#         self.unit_step = 0.02

#         self.input_beam_file = self.particle_dir / 'section_L1.npz'
#         self.output_beam_file = self.particle_dir / 'section_BC1.npz'
#         self.tws_file = self.tws_dir / "tws_section_BC1.npz"

#         acc2_stop = l1.stlat_182_b1
#         bc1_stop = l1.tora_203_b1
#         # init tracking lattice
#         self.lattice = MagneticLattice(l1.cell, start=acc2_stop, stop=bc1_stop, method=self.method)

#         # init physics processes
#         csr = CSR()
#         csr.step = 1
#         csr.n_bin = CSRBin
#         csr.sigma_min = Sig_Z[2]*CSRSigmaFactor
#         csr.traj_step = 0.0005
#         csr.apply_step = 0.001

#         sc = SpaceCharge()
#         sc.step = 40
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh
#         match_bc1 = acc2_stop
#         self.add_physics_process(csr, start=match_bc1, stop=bc1_stop)
#         self.add_physics_process(sc, start=match_bc1, stop=bc1_stop)
#         self.dipoles = [l1.bb_182_b1, l1.bb_191_b1, l1.bb_193_b1, l1.bb_202_b1]
#         self.dipole_len = 0.5
#         self.bc_gap=8.5


class L2(FELSection):
    def __init__(self):
        cell = MachineSequence(l1.cell + l2.cell)
        bc1_stop_l2_start = "bc1_stop_l2_start"
        l2_stop_bc2_start = "l2_stop_bc2_start"
        sequence = cell[bc1_stop_l2_start:l2_stop_bc2_start]
        super().__init__(sequence)

    def setup_physics(self):
        bc1_stop_l2_start = "bc1_stop_l2_start"
        l2_stop_bc2_start = "l2_stop_bc2_start"
        sc = make_space_charge(step=100, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)
        wake = make_wake(
            WAKES_PATH / "RF/mod_TESLA_MODULE_WAKE_TAYLOR.dat", factor=4 * 3, step=200
        )
        wake2 = make_wake(WAKES_PATH / "mod_wake_0179.810_0370.840_MONO.dat", factor=1)

        first_cavity_l2 = "C.A3.1.1.L2"
        last_cavity_l2 = "C.A5.4.8.L2"

        self.add_physics_process(
            BEAM_SMOOTHER, start=bc1_stop_l2_start, stop=bc1_stop_l2_start
        )
        self.add_physics_process(sc, start=bc1_stop_l2_start, stop=l2_stop_bc2_start)
        self.add_physics_process(wake, start=first_cavity_l2, stop=last_cavity_l2)
        self.add_physics_process(wake2, start=l2_stop_bc2_start, stop=l2_stop_bc2_start)


# class L2(FELSection):

#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)

#         # setting parameters
#         self.lattice_name = 'L2'
#         self.unit_step = 0.02

#         self.input_beam_file = self.particle_dir / 'section_BC1.npz'
#         self.output_beam_file = self.particle_dir / 'section_L2.npz'
#         self.tws_file = self.tws_dir / "tws_section_L2.npz"

#         bc1_stop = l1.tora_203_b1
#         acc3t5_stop = l2.stlat_393_b2
#         # init tracking lattice
#         self.lattice = MagneticLattice(l1.cell + l2.cell, start=bc1_stop, stop=acc3t5_stop, method=self.method)

#         # init physics processes
#         smooth = SmoothBeam()
#         smooth.mslice = SmoothPar

#         sc = SpaceCharge()
#         sc.step = 100
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh
#         wake = Wake()
#         wake.wake_table = WakeTable(WAKES_PATH / 'mod_TESLA_MODULE_WAKE_TAYLOR.dat')
#         wake.factor = 4 * 3
#         wake.step = 200
#         wake_add = Wake()
#         wake_add.wake_table = WakeTable(WAKES_PATH / 'mod_wake_0179.810_0370.840_MONO.dat')
#         wake_add.factor = 1
#         self.add_physics_process(smooth, start=bc1_stop, stop=bc1_stop)
#         self.add_physics_process(sc, start=bc1_stop, stop=acc3t5_stop)
#         self.add_physics_process(wake, start=l2.c_a3_1_1_l2, stop=l2.c_a5_4_8_l2)
#         self.add_physics_process(wake_add, start=acc3t5_stop, stop=acc3t5_stop)


class BC2(FELSection):
    def __init__(self):
        l2_cell = MachineSequence(l2.cell)

        l2_stop_bc2_start = "l2_stop_bc2_start"
        bc2_stop_b2d_start = "bc2_stop_b2d_start"
        # init tracking lattice
        sequence = l2_cell[l2_stop_bc2_start:bc2_stop_b2d_start]
        super().__init__(sequence)

    def setup_physics(self):
        l2_stop_bc2_start = "l2_stop_bc2_start"
        bc2_stop_b2d_start = "bc2_stop_b2d_start"

        csr = make_csr(
            step=1,
            n_bin=CSR_N_BIN,
            sigma_min=Sig_Z * CSR_SIGMA_FACTOR,
            traj_step=0.0005,
            apply_step=0.001,
        )
        sc = make_space_charge(step=50, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)

        self.add_physics_process(csr, start=l2_stop_bc2_start, stop=bc2_stop_b2d_start)
        self.add_physics_process(sc, start=l2_stop_bc2_start, stop=bc2_stop_b2d_start)


# class BC2(FELSection):

#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)

#         # setting parameters
#         self.lattice_name = 'BC2'
#         self.dipoles = [l2.bb_393_b2, l2.bb_402_b2, l2.bb_404_b2, l2.bb_413_b2]
#         self.dipole_len = 0.5
#         self.bc_gap = 8.5

#         self.unit_step = 0.02

#         self.input_beam_file = self.particle_dir / 'section_L2.npz'
#         self.output_beam_file = self.particle_dir / 'section_BC2.npz'
#         self.tws_file = self.tws_dir / "tws_section_BC2.npz"

#         acc3t5_stop = l2.stlat_393_b2
#         bc2_stop = l2.tora_415_b2
#         # init tracking lattice
#         self.lattice = MagneticLattice(l2.cell, start=acc3t5_stop, stop=bc2_stop, method=self.method)

#         # init physics processes

#         csr = CSR()
#         csr.step=1
#         csr.n_bin = CSRBin
#         csr.sigma_min = Sig_Z[3]*CSRSigmaFactor
#         csr.traj_step = 0.0005
#         csr.apply_step = 0.001
#         #csr.rk_traj = True
#         #csr.energy = 2.4

#         sc = SpaceCharge()
#         sc.step = 50 # 50
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh

#         self.add_physics_process(csr, start=acc3t5_stop, stop=bc2_stop)
#         self.add_physics_process(sc, start=acc3t5_stop, stop=bc2_stop)


class B2D(FELSection):
    def __init__(self):
        cell = MachineSequence(l2.cell + b2d.cell)

        bc2_stop_b2d_start = "bc2_stop_b2d_start"
        sequence = cell[bc2_stop_b2d_start:]
        super().__init__(sequence)

    def setup_physics(self):
        bc2_stop_b2d_start = "bc2_stop_b2d_start"
        csr = make_csr(
            step=1,
            n_bin=CSR_N_BIN,
            sigma_min=Sig_Z * CSR_SIGMA_FACTOR,
            traj_step=0.0005,
            apply_step=0.001,
        )
        sc = make_space_charge(step=50, nmesh_xyz=SC_MESH, random_mesh=SC_RANDOM_MESH)

        self.add_physics_process(csr, start=bc2_stop_b2d_start, stop="BPMA.471.B2D")
        self.add_physics_process(sc, start=bc2_stop_b2d_start, stop=self.sequence[-1].id)


# class L3(FELSection):

#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)

#         # setting parameters
#         self.lattice_name = 'L3'
#         self.unit_step = 0.2

#         self.input_beam_file = self.particle_dir + 'section_BC2.npz'
#         self.output_beam_file = self.particle_dir + 'section_L3.npz'
#         self.tws_file = self.tws_dir + "tws_section_L3.npz"

#         bc2_stop = l2.tora_415_b2
#         acc6t26_stop =l3.stop_l3
#         #acc6t26_stop=l3.qd_470_b2

#         # init tracking lattice
#         self.lattice = MagneticLattice(l2.cell + l3.cell + cl.cell, start=bc2_stop, stop=acc6t26_stop, method=self.method)

#         # init physics processes
#         smooth = SmoothBeam()
#         smooth.mslice = SmoothPar


#         sc = SpaceCharge()
#         sc.step = 5
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh

#         wake = Wake()
#         wake.wake_table = WakeTable('accelerator/wakes/RF/mod_TESLA_MODULE_WAKE_TAYLOR.dat')
#         wake.factor = 4 * 21
#         wake.step = 10
#         wake.w_sampling = 1000
#         wake_add = Wake()
#         wake_add.wake_table = WakeTable('accelerator/wakes/mod_wake_0391.350_1629.700_MONO.dat')
#         wake_add.factor = 1

#         app = PhaseSpaceAperture()
#         app.taumin = -5
#         app.taumax = 3

#         self.add_physics_process(app, start=bc2_stop, stop=bc2_stop)
#         self.add_physics_process(smooth, start=bc2_stop, stop=bc2_stop)
#         self.add_physics_process(wake, start=l3.c_a6_1_1_l3,stop=l3.c_a25_4_8_l3)
#         self.add_physics_process(sc, start=bc2_stop, stop=acc6t26_stop)
#         self.add_physics_process(wake_add, start=acc6t26_stop, stop=acc6t26_stop)


# class CL1(FELSection):
#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)

#         # setting parameters
#         self.lattice_name = 'CL1'
#         self.unit_step = 0.2


#         self.input_beam_file = self.particle_dir + 'section_L3.npz'
#         self.output_beam_file = self.particle_dir + 'section_CL1.npz'
#         self.tws_file = self.tws_dir + "tws_section_CL1.npz"

#         #acc6t26_stop = cl.match_1673_cl
#         acc6t26_stop = l3.stop_l3
#         collimator1_stop = cl.bpma_1746_cl
#         # init tracking lattice
#         self.lattice = MagneticLattice(l3.cell + cl.cell, start=acc6t26_stop,stop=collimator1_stop, method=self.method)

#         # init physics processes

#         sc = SpaceCharge()
#         sc.step = 10
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh
#         self.add_physics_process(sc, start=acc6t26_stop, stop=collimator1_stop)

#         csr = CSR()
#         csr.traj_step = 0.0005
#         csr.apply_step = 0.001
#         csr.n_bin = CSRBin
#         csr.sigma_min = Sig_Z[3]*CSRSigmaFactor
#         self.add_physics_process(csr, start=acc6t26_stop, stop=collimator1_stop)

#         if bISR:
#             LD = cl.be_1678_cl.l;  teta = cl.be_1678_cl.angle; ro = LD / np.sin(teta);
#             sre1 = SpontanRadEffects(); sre1.radius = ro;  sre1.type = 'dipole';
#             self.add_physics_process(sre1, cl.M1be_1678_cl, cl.M2be_1678_cl)
#             LD = cl.bl_1688_cl.l;        teta = cl.bl_1688_cl.angle;        ro = LD / np.sin(teta);
#             sre2 = SpontanRadEffects();        sre2.radius = ro;        sre2.type = 'dipole';
#             self.add_physics_process(sre2, cl.M1be_1688_cl, cl.M2be_1688_cl)
#             LD = cl.bl_1695_cl.l;        teta = cl.bl_1695_cl.angle;        ro = LD / np.sin(teta);
#             sre3 = SpontanRadEffects();        sre3.radius = ro;        sre3.type = 'dipole';
#             self.add_physics_process(sre3, cl.M1be_1695_cl, cl.M2be_1695_cl)
#             LD = cl.be_1705_cl.l;        teta = cl.be_1705_cl.angle;        ro = LD / np.sin(teta);
#             sre4 = SpontanRadEffects();        sre4.radius = ro;        sre4.type = 'dipole';
#             self.add_physics_process(sre4, cl.M1be_1705_cl, cl.M2be_1705_cl)
#             LD = cl.be_1714_cl.l;  teta = cl.be_1714_cl.angle; ro = LD / np.sin(teta);
#             sre5 = SpontanRadEffects(); sre5.radius = ro;  sre5.type = 'dipole';
#             self.add_physics_process(sre5, cl.M1be_1714_cl, cl.M2be_1714_cl)
#             LD = cl.bl_1724_cl.l;        teta = cl.bl_1724_cl.angle;        ro = LD / np.sin(teta);
#             sre6 = SpontanRadEffects();        sre6.radius = ro;        sre6.type = 'dipole';
#             self.add_physics_process(sre6, cl.M1be_1724_cl, cl.M2be_1724_cl)
#             LD = cl.bl_1731_cl.l;        teta = cl.bl_1731_cl.angle;        ro = LD / np.sin(teta);
#             sre7 = SpontanRadEffects();        sre7.radius = ro;        sre7.type = 'dipole';
#             self.add_physics_process(sre7, cl.M1be_1731_cl, cl.M2be_1731_cl)
#             LD = cl.be_1741_cl.l;  teta = cl.be_1741_cl.angle; ro = LD / np.sin(teta);
#             sre8 = SpontanRadEffects(); sre8.radius = ro;  sre8.type = 'dipole';
#             self.add_physics_process(sre8, cl.M1be_1741_cl, cl.M2be_1741_cl)

# class CL2(FELSection):
#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)

#         # setting parameters
#         self.lattice_name = 'CL2'
#         self.unit_step = 1


#         self.input_beam_file = self.particle_dir + 'section_CL1.npz'
#         self.output_beam_file = self.particle_dir + 'section_CL2.npz'
#         self.tws_file = self.tws_dir + "tws_section_CL2.npz"

#         collimator1_stop = cl.bpma_1746_cl
#         collimator2_stop = cl.bpma_1783_cl
#         # init tracking lattice
#         self.lattice = MagneticLattice(cl.cell, start=collimator1_stop, stop=collimator2_stop, method=self.method)

#         # init physics processes

#         sc = SpaceCharge()
#         sc.step = 1
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh
#         self.add_physics_process(sc, start=collimator1_stop, stop=collimator2_stop)


# class CL3(FELSection):
#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)

#         # setting parameters
#         self.lattice_name = 'CL3'
#         self.unit_step = 0.2


#         self.input_beam_file = self.particle_dir + 'section_CL2.npz'
#         self.output_beam_file = self.particle_dir + 'section_CL3.npz'
#         self.tws_file = self.tws_dir + "tws_section_CL3.npz"

#         collimator2_stop = cl.bpma_1783_cl
#         collimator3_stop = cl.ensec_1854_cl
#         # init tracking lattice
#         self.lattice = MagneticLattice(cl.cell, start=collimator2_stop, stop=collimator3_stop, method=self.method)

#         # init physics processes

#         sc = SpaceCharge()
#         sc.step = 10
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh

#         csr = CSR()
#         csr.traj_step = 0.0005
#         csr.apply_step = 0.001
#         csr.n_bin = CSRBin
#         csr.sigma_min = Sig_Z[3]*CSRSigmaFactor

#         wake_add = Wake()
#         wake_add.wake_table = WakeTable('accelerator/wakes/mod_wake_1629.700_1831.200_MONO.dat')
#         wake_add.factor = 1

#         self.add_physics_process(csr, start=collimator2_stop, stop=collimator3_stop)
#         self.add_physics_process(sc, start=collimator2_stop, stop=collimator3_stop)
#         self.add_physics_process(wake_add, start=collimator3_stop, stop=collimator3_stop)

#         if bISR:
#             LD = cl.be_1786_cl.l;  teta = cl.be_1786_cl.angle; ro = LD / np.sin(teta);
#             sre1 = SpontanRadEffects(); sre1.radius = ro;  sre1.type = 'dipole';
#             self.add_physics_process(sre1, cl.M1be_1786_cl, cl.M2be_1786_cl)
#             LD = cl.bl_1796_cl.l;        teta = cl.bl_1796_cl.angle;        ro = LD / np.sin(teta);
#             sre2 = SpontanRadEffects();        sre2.radius = ro;        sre2.type = 'dipole';
#             self.add_physics_process(sre2, cl.M1be_1796_cl, cl.M2be_1796_cl)
#             LD = cl.bl_1803_cl.l;        teta = cl.bl_1803_cl.angle;        ro = LD / np.sin(teta);
#             sre3 = SpontanRadEffects();        sre3.radius = ro;        sre3.type = 'dipole';
#             self.add_physics_process(sre3, cl.M1be_1803_cl, cl.M2be_1803_cl)
#             LD = cl.be_1813_cl.l;        teta = cl.be_1813_cl.angle;        ro = LD / np.sin(teta);
#             sre4 = SpontanRadEffects();        sre4.radius = ro;        sre4.type = 'dipole';
#             self.add_physics_process(sre4, cl.M1be_1813_cl, cl.M2be_1813_cl)
#             LD = cl.be_1822_cl.l;  teta = cl.be_1822_cl.angle; ro = LD / np.sin(teta);
#             sre5 = SpontanRadEffects(); sre5.radius = ro;  sre5.type = 'dipole';
#             self.add_physics_process(sre5, cl.M1be_1822_cl, cl.M2be_1822_cl)
#             LD = cl.bl_1832_cl.l;        teta = cl.bl_1832_cl.angle;        ro = LD / np.sin(teta);
#             sre6 = SpontanRadEffects();        sre6.radius = ro;        sre6.type = 'dipole';
#             self.add_physics_process(sre6, cl.M1be_1832_cl, cl.M2be_1832_cl)
#             LD = cl.bl_1839_cl.l;        teta = cl.bl_1839_cl.angle;        ro = LD / np.sin(teta);
#             sre7 = SpontanRadEffects();        sre7.radius = ro;        sre7.type = 'dipole';
#             self.add_physics_process(sre7, cl.M1be_1839_cl, cl.M2be_1839_cl)
#             LD = cl.be_1849_cl.l;  teta = cl.be_1849_cl.angle; ro = LD / np.sin(teta);
#             sre8 = SpontanRadEffects(); sre8.radius = ro;  sre8.type = 'dipole';
#             self.add_physics_process(sre8, cl.M1be_1849_cl, cl.M2be_1849_cl)

# class STN10(FELSection):
#     def __init__(self, data_dir, *args, **kwargs):
#         super().__init__(data_dir)
#         # setting parameters
#         self.lattice_name = 'ST10'
#         self.unit_step = 10
#         self.input_beam_file = self.particle_dir + 'section_CL3.npz'
#         self.output_beam_file = self.particle_dir + 'section_STN10.npz'
#         self.tws_file = self.tws_dir + "tws_section_STN10.npz"
#         collimator3_stop = cl.ensec_1854_cl
#         stN10_stop = sase1.ensec_2235_t2
#         #stN10_stop = cl.ensub_1980_tl
#         # init tracking lattice
#         self.lattice = MagneticLattice(cl.cell + tl34_sase1.cell + sase1.cell, start=collimator3_stop, stop=stN10_stop, method=self.method)
#         # init physics processes
#         sc = SpaceCharge()
#         sc.step = 1
#         sc.nmesh_xyz = SCmesh
#         sc.random_mesh = bRandomMesh
#         wake_add = Wake()
#         wake_add.wake_table = WakeTable('accelerator/wakes/mod_wake_1831.200_2035.190_MONO.dat')
#         wake_add.factor = 1
#         wake_add1 = Wake()
#         wake_add1.wake_table = WakeTable('accelerator/wakes/mod_wake_2035.190_2213.000_MONO.dat')
#         wake_add1.factor = 1

#         self.add_physics_process(sc, start=collimator3_stop, stop=stN10_stop)
#         self.add_physics_process(wake_add, start=collimator3_stop, stop=collimator3_stop)
#         self.add_physics_process(wake_add1, start=stN10_stop, stop=stN10_stop)


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
