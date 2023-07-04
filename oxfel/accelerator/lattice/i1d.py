from ocelot import * 

# Drifts
d_93 = Drift(l=0.002879000000000076, eid='D_93')
d_94 = Drift(l=0.969028999999999, eid='D_94')
d_95 = Drift(l=0.15615000000000379, eid='D_95')
d_96 = Drift(l=0.4060000000000059, eid='D_96')
d_97 = Drift(l=0.15014999999999645, eid='D_97')
d_98 = Drift(l=0.24921000000000504, eid='D_98')
d_99 = Drift(l=0.22500000000000142, eid='D_99')
d_100 = Drift(l=0.125, eid='D_100')
d_101 = Drift(l=0.28522999999999854, eid='D_101')
d_102 = Drift(l=0.25200000000000244, eid='D_102')
d_103 = Drift(l=1.695999999999998, eid='D_103')
d_104 = Drift(l=0.09799999999999898, eid='D_104')

# Quadrupoles
qi_63_i1d = Quadrupole(l=0.2377, k1=4.401795, eid='QI.63.I1D')
qi_64_i1d = Quadrupole(l=0.2377, eid='QI.64.I1D')

# SBends
bb_62_i1d = SBend(l=0.5, angle=0.5235987756, e1=0.261799388, e2=0.261799388, eid='BB.62.I1D')

# Monitors
bpma_63_i1d = Monitor(eid='BPMA.63.I1D')
bpmd_64_i1d = Monitor(eid='BPMD.64.I1D')

# Markers
stsec_62_i1d = Marker(eid='STSEC.62.I1D')
otrc_64_i1d = Marker(eid='OTRC.64.I1D')
otrd_64_i1d = Marker(eid='OTRD.64.I1D')
torc_64_i1d = Marker(eid='TORC.64.I1D')
duflange_65_i1d = Marker(eid='DUFLANGE.65.I1D')
duconcrete_65_i1d = Marker(eid='DUCONCRETE.65.I1D')
bhm_66_i1d = Marker(eid='BHM.66.I1D')
duabsorb_66_i1d = Marker(eid='DUABSORB.66.I1D')
ensec_66_i1d = Marker(eid='ENSEC.66.I1D')

# Lattice 
cell = [stsec_62_i1d, d_93, bb_62_i1d, d_94, qi_63_i1d, d_95, bpma_63_i1d, d_96, otrc_64_i1d, 
d_97, qi_64_i1d, d_98, otrd_64_i1d, d_99, torc_64_i1d, d_100, bpmd_64_i1d, d_101, duflange_65_i1d, 
d_102, duconcrete_65_i1d, d_103, bhm_66_i1d, d_104, duabsorb_66_i1d, ensec_66_i1d]
