from ocelot import * 

# Drifts
d_671 = Drift(l=0.0014499999999770807, eid='D_671')
d_672 = Drift(l=0.5814499999999612, eid='D_672')
d_673 = Drift(l=0.17394999999999072, eid='D_673')
d_674 = Drift(l=0.18394999999998163, eid='D_674')
d_675 = Drift(l=1.0437200000000075, eid='D_675')
d_676 = Drift(l=0.1362799999999993, eid='D_676')
d_677 = Drift(l=0.17999999999994998, eid='D_677')
d_678 = Drift(l=0.20000000000004547, eid='D_678')
d_679 = Drift(l=0.9839499999999362, eid='D_679')
d_680 = Drift(l=0.43628000000001066, eid='D_680')
d_681 = Drift(l=0.4191200000000208, eid='D_681')
d_682 = Drift(l=0.5163769999999772, eid='D_682')
d_683 = Drift(l=0.16901999999998907, eid='D_683')
d_684 = Drift(l=0.25897999999995136, eid='D_684')
d_685 = Drift(l=0.224899999999991, eid='D_685')
d_686 = Drift(l=0.2340199999999868, eid='D_686')
d_687 = Drift(l=0.6636799999999994, eid='D_687')
d_688 = Drift(l=0.22500000000002274, eid='D_688')
d_689 = Drift(l=0.125, eid='D_689')
d_690 = Drift(l=0.2852699999999686, eid='D_690')
d_691 = Drift(l=0.7280000000000086, eid='D_691')

# Quadrupoles
qf_469_b2d = Quadrupole(l=0.5321, k1=-2.19942871499906, eid='QF.469.B2D')
qe_471_b2d = Quadrupole(l=0.24, k1=1.335539651, eid='QE.471.B2D')
qf_472_b2d = Quadrupole(l=0.5321, k1=-2.19942871499906, eid='QF.472.B2D')
qf_476_b2d = Quadrupole(l=0.5321, k1=3.13097893, eid='QF.476.B2D')
qf_477_b2d = Quadrupole(l=0.5321, k1=0.7703645572993798, eid='QF.477.B2D')

# SBends
bg_467_b2d = SBend(l=1.5971, angle=0.2094395102, tilt=1.570796327, eid='BG.467.B2D')
bg_474_b2d = SBend(l=1.5971, angle=-0.2094395102, tilt=1.570796327, eid='BG.474.B2D')

# Hcors
cfx_470_b2d = Hcor(l=0.1, eid='CFX.470.B2D')
cfx_477_b2d = Hcor(l=0.1, eid='CFX.477.B2D')

# Vcors
cfy_468_b2d = Vcor(l=0.1, eid='CFY.468.B2D')
cfy_471_b2d = Vcor(l=0.1, eid='CFY.471.B2D')
cfy_476_b2d = Vcor(l=0.1, eid='CFY.476.B2D')

# Monitors
bpma_469_b2d = Monitor(eid='BPMA.469.B2D')
bpma_471_b2d = Monitor(eid='BPMA.471.B2D')
bpma_477_b2d = Monitor(eid='BPMA.477.B2D')
bpmd_479_b2d = Monitor(eid='BPMD.479.B2D')

# Markers
stsec_466_b2d = Marker(eid='STSEC.466.B2D')
otra_473_b2d = Marker(eid='OTRA.473.B2D')
otrd_478_b2d = Marker(eid='OTRD.478.B2D')
torc_479_b2d = Marker(eid='TORC.479.B2D')
duflange_479_b2d = Marker(eid='DUFLANGE.479.B2D')
duabsorb_480_b2d = Marker(eid='DUABSORB.480.B2D')
ensec_480_b2d = Marker(eid='ENSEC.480.B2D')

# Lattice 
cell = [stsec_466_b2d, d_671, bg_467_b2d, d_672, cfy_468_b2d, d_673, qf_469_b2d, d_674, bpma_469_b2d, 
d_675, cfx_470_b2d, d_676, qe_471_b2d, d_677, bpma_471_b2d, d_678, cfy_471_b2d, d_679, qf_472_b2d, 
d_680, otra_473_b2d, d_681, bg_474_b2d, d_682, cfy_476_b2d, d_683, qf_476_b2d, d_684, bpma_477_b2d, 
d_685, cfx_477_b2d, d_686, qf_477_b2d, d_687, otrd_478_b2d, d_688, torc_479_b2d, d_689, bpmd_479_b2d, 
d_690, duflange_479_b2d, d_691, duabsorb_480_b2d, ensec_480_b2d]
