from ocelot import *

# Drifts
d_395 = Drift(l=0.254587000000015, eid="D_395")
d_396 = Drift(l=1.3195570000000032, eid="D_396")
d_397 = Drift(l=0.15064999999998463, eid="D_397")
d_398 = Drift(l=0.13164999999997917, eid="D_398")
d_399 = Drift(l=1.231650000000002, eid="D_399")
d_400 = Drift(l=0.6896500000000003, eid="D_400")
d_401 = Drift(l=0.1810000000000116, eid="D_401")
d_402 = Drift(l=2.584000000000003, eid="D_402")
d_403 = Drift(l=0.14829999999997767, eid="D_403")
d_404 = Drift(l=0.15430000000000632, eid="D_404")
d_405 = Drift(l=0.4252999999999929, eid="D_405")
d_406 = Drift(l=0.724000000000018, eid="D_406")

# Quadrupoles
qd_231_b1d = Quadrupole(l=0.2367, k1=-3.0, eid="QD.231.B1D")
qd_232_b1d = Quadrupole(l=0.2367, eid="QD.232.B1D")

# SBends
bb_229_b1d = SBend(
    l=0.5, angle=0.2094395102, e2=0.20943951, tilt=1.570796327, eid="BB.229.B1D"
)

# Hcors
ccx_233_b1d = Hcor(l=0.1, eid="CCX.233.B1D")

# Vcors
ccy_231_b1d = Vcor(l=0.1, eid="CCY.231.B1D")

# Monitors
bpma_231_b1d = Monitor(eid="BPMA.231.B1D")
bpma_233_b1d = Monitor(eid="BPMA.233.B1D")
bpma_236_b1d = Monitor(eid="BPMA.236.B1D")

# Markers
stsec_229_b1d = Marker(eid="STSEC.229.B1D")
otrc_236_b1d = Marker(eid="OTRC.236.B1D")
tora_236_b1d = Marker(eid="TORA.236.B1D")
duflange_237_b1d = Marker(eid="DUFLANGE.237.B1D")
duabsorb_237_b1d = Marker(eid="DUABSORB.237.B1D")
ensec_237_b1d = Marker(eid="ENSEC.237.B1D")

# Lattice
cell = [
    stsec_229_b1d,
    d_395,
    bb_229_b1d,
    d_396,
    bpma_231_b1d,
    d_397,
    qd_231_b1d,
    d_398,
    ccy_231_b1d,
    d_399,
    qd_232_b1d,
    d_400,
    bpma_233_b1d,
    d_401,
    ccx_233_b1d,
    d_402,
    otrc_236_b1d,
    d_403,
    tora_236_b1d,
    d_404,
    bpma_236_b1d,
    d_405,
    duflange_237_b1d,
    d_406,
    duabsorb_237_b1d,
    ensec_237_b1d,
]
