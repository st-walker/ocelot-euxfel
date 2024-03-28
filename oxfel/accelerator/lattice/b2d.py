from ocelot import * 

# Drifts
d_687 = Drift(l=0.0014499999999770807, eid='D_687')
d_688 = Drift(l=0.581450000000018, eid='D_688')
d_689 = Drift(l=0.17394999999999072, eid='D_689')
d_690 = Drift(l=0.18394999999998163, eid='D_690')
d_691 = Drift(l=1.0437200000000075, eid='D_691')
d_692 = Drift(l=0.1362799999999993, eid='D_692')
d_693 = Drift(l=0.17999999999994998, eid='D_693')
d_694 = Drift(l=0.19999999999998863, eid='D_694')
d_695 = Drift(l=0.983949999999993, eid='D_695')
d_696 = Drift(l=0.43628000000001066, eid='D_696')
d_697 = Drift(l=0.4191200000000208, eid='D_697')
d_698 = Drift(l=0.5163770000000341, eid='D_698')
d_699 = Drift(l=0.16901999999993222, eid='D_699')
d_700 = Drift(l=0.2589800000000082, eid='D_700')
d_701 = Drift(l=0.224899999999991, eid='D_701')
d_702 = Drift(l=0.2340199999999868, eid='D_702')
d_703 = Drift(l=0.6636799999999994, eid='D_703')
d_704 = Drift(l=0.22500000000002274, eid='D_704')
d_705 = Drift(l=0.125, eid='D_705')
d_706 = Drift(l=0.2852699999999686, eid='D_706')
d_707 = Drift(l=0.7280000000000086, eid='D_707')

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
cell = (stsec_466_b2d, d_687, bg_467_b2d, d_688, cfy_468_b2d, d_689, qf_469_b2d, d_690, bpma_469_b2d, 
d_691, cfx_470_b2d, d_692, qe_471_b2d, d_693, bpma_471_b2d, d_694, cfy_471_b2d, d_695, qf_472_b2d, 
d_696, otra_473_b2d, d_697, bg_474_b2d, d_698, cfy_476_b2d, d_699, qf_476_b2d, d_700, bpma_477_b2d, 
d_701, cfx_477_b2d, d_702, qf_477_b2d, d_703, otrd_478_b2d, d_704, torc_479_b2d, d_705, bpmd_479_b2d, 
d_706, duflange_479_b2d, d_707, duabsorb_480_b2d, ensec_480_b2d)
