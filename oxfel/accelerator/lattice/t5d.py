from ocelot import * 

# Drifts
d_2585 = Drift(l=1.0, eid='D_2585')
d_2586 = Drift(l=12.471400000000358, eid='D_2586')
d_2587 = Drift(l=15.514999999999873, eid='D_2587')
d_2588 = Drift(l=0.20895000000018626, eid='D_2588')
d_2589 = Drift(l=0.15394999999944048, eid='D_2589')
d_2590 = Drift(l=11.204999999999927, eid='D_2590')
d_2591 = Drift(l=0.20894999999973152, eid='D_2591')
d_2592 = Drift(l=0.15394999999944048, eid='D_2592')
d_2593 = Drift(l=60.20499999999993, eid='D_2593')
d_2594 = Drift(l=0.20894999999973152, eid='D_2594')
d_2595 = Drift(l=0.15394999999989523, eid='D_2595')
d_2596 = Drift(l=5.359999999999673, eid='D_2596')
d_2597 = Drift(l=0.15394999999944048, eid='D_2597')
d_2598 = Drift(l=0.26395000000002256, eid='D_2598')
d_2599 = Drift(l=0.09600000000000364, eid='D_2599')
d_2600 = Drift(l=0.11079999999992651, eid='D_2600')
d_2601 = Drift(l=0.5158219999998437, eid='D_2601')
d_2602 = Drift(l=9.900000031848322e-05, eid='D_2602')
d_2603 = Drift(l=0.5001979999997275, eid='D_2603')
d_2604 = Drift(l=1.9300990000001548, eid='D_2604')
d_2605 = Drift(l=0.14239999999972497, eid='D_2605')
d_2606 = Drift(l=0.197400000000016, eid='D_2606')
d_2607 = Drift(l=0.5535000000004402, eid='D_2607')
d_2608 = Drift(l=0.9034999999998945, eid='D_2608')
d_2609 = Drift(l=0.8535000000001673, eid='D_2609')
d_2610 = Drift(l=0.5534999999999854, eid='D_2610')
d_2611 = Drift(l=0.197400000000016, eid='D_2611')
d_2612 = Drift(l=0.14240000000017972, eid='D_2612')
d_2613 = Drift(l=1.9300989999997, eid='D_2613')
d_2614 = Drift(l=0.500199000000066, eid='D_2614')
d_2615 = Drift(l=0.6117990000002465, eid='D_2615')
d_2616 = Drift(l=0.16069999999990614, eid='D_2616')
d_2617 = Drift(l=0.34480000000030486, eid='D_2617')
d_2618 = Drift(l=0.34480000000030486, eid='D_2618')
d_2619 = Drift(l=0.34480000000030486, eid='D_2619')
d_2620 = Drift(l=0.21439999999984138, eid='D_2620')
d_2621 = Drift(l=0.5579999999999927, eid='D_2621')
d_2622 = Drift(l=0.49999999999954525, eid='D_2622')
d_2623 = Drift(l=0.8219999999996617, eid='D_2623')
d_2624 = Drift(l=0.20000000000027285, eid='D_2624')
d_2625 = Drift(l=0.4499999999998181, eid='D_2625')
d_2626 = Drift(l=0.31300000000010186, eid='D_2626')
d_2627 = Drift(l=3.869999999999891, eid='D_2627')
d_2628 = Drift(l=0.3223000000002685, eid='D_2628')
d_2629 = Drift(l=0.8577999999997701, eid='D_2629')
d_2630 = Drift(l=1.073240000000169, eid='D_2630')
d_2631 = Drift(l=2.0144000000000233, eid='D_2631')

# Quadrupoles
qe_3052_t5d = Quadrupole(l=0.24, k1=0.1781107015, eid='QE.3052.T5D')
qf_3068_t5d = Quadrupole(l=0.5321, k1=-0.17408057190001877, eid='QF.3068.T5D')
qf_3081_t5d = Quadrupole(l=0.5321, k1=0.041949297430934035, eid='QF.3081.T5D')
qf_3142_t5d = Quadrupole(l=0.5321, k1=0.1523450185002819, eid='QF.3142.T5D')
qf_3148_t5d = Quadrupole(l=0.5321, k1=-0.08121677321931968, eid='QF.3148.T5D')
qk_3158_t5d = Quadrupole(l=1.0552, k1=-0.172994430300417, eid='QK.3158.T5D')
qk_3163_t5d = Quadrupole(l=1.0552, k1=-0.172994430300417, eid='QK.3163.T5D')
qk_3172_t5d = Quadrupole(l=1.0552, k1=-0.18514803489954512, eid='QK.3172.T5D')
qk_3174_t5d = Quadrupole(l=1.0552, k1=-0.18514803489954512, eid='QK.3174.T5D')
qk_3175_t5d = Quadrupole(l=1.0552, k1=-0.18514803489954512, eid='QK.3175.T5D')
qk_3177_t5d = Quadrupole(l=1.0552, k1=-0.18514803489954512, eid='QK.3177.T5D')

# SBends
bv_3151_t5d = SBend(l=2.5, angle=0.043633231, e1=0.021816616, e2=0.021816616, tilt=1.570796327, eid='BV.3151.T5D')
bv_3154_t5d = SBend(l=2.5, angle=0.043633231, e1=0.021816616, e2=0.021816616, tilt=1.570796327, eid='BV.3154.T5D')
bv_3167_t5d = SBend(l=2.5, angle=0.043633231, e1=0.021816616, e2=0.021816616, tilt=1.570796327, eid='BV.3167.T5D')
bv_3170_t5d = SBend(l=2.5, angle=0.043633231, e1=0.021816616, e2=0.021816616, tilt=1.570796327, eid='BV.3170.T5D')

# RBends
sweep_3178_t5d = RBend(l=0.64, tilt=-1.570796327, eid='SWEEP.3178.T5D')
sweep_3179_t5d = RBend(l=0.64, eid='SWEEP.3179.T5D')

# Sextupoles
sk_3159_t5d = Sextupole(l=0.343, k2=2.020228531, tilt=1.570796327, eid='SK.3159.T5D')
sk_3161_t5d = Sextupole(l=0.343, k2=2.020228531, tilt=1.570796327, eid='SK.3161.T5D')

# Hcors
cfx_3081_t5d = Hcor(l=0.1, eid='CFX.3081.T5D')
cfx_3142_t5d = Hcor(l=0.1, eid='CFX.3142.T5D')
cnx_3157_t5d = Hcor(l=0.3, eid='CNX.3157.T5D')

# Vcors
cfy_3069_t5d = Vcor(l=0.1, eid='CFY.3069.T5D')
cfy_3148_t5d = Vcor(l=0.1, eid='CFY.3148.T5D')
cny_3164_t5d = Vcor(l=0.3, eid='CNY.3164.T5D')

# Monitors
bpma_3068_t5d = Monitor(eid='BPMA.3068.T5D')
bpma_3080_t5d = Monitor(eid='BPMA.3080.T5D')
bpma_3141_t5d = Monitor(eid='BPMA.3141.T5D')
bpmf_3149_t5d = Monitor(eid='BPMF.3149.T5D')
bpmd_3159_t5d = Monitor(eid='BPMD.3159.T5D')
bpmd_3162_t5d = Monitor(eid='BPMD.3162.T5D')
bpmd_3172_t5d = Monitor(eid='BPMD.3172.T5D')
bpmd_3177_t5d = Monitor(eid='BPMD.3177.T5D')
bpmd_3180_t5d = Monitor(eid='BPMD.3180.T5D')
bpmw_3185_t5d = Monitor(eid='BPMW.3185.T5D')

# Markers
stsec_3039_t5d = Marker(eid='STSEC.3039.T5D')
stsub_3039_t5d = Marker(eid='STSUB.3039.T5D')
tora_3040_t5d = Marker(eid='TORA.3040.T5D')
midbpmf_3149_t5d = Marker(eid='MIDBPMF.3149.T5D')
tora_3149_t5d = Marker(eid='TORA.3149.T5D')
ensub_3149_t5d = Marker(eid='ENSUB.3149.T5D')
stsub_3149_t5d = Marker(eid='STSUB.3149.T5D')
otrd_3160_t5d = Marker(eid='OTRD.3160.T5D')
otrd_3181_t5d = Marker(eid='OTRD.3181.T5D')
torc_3181_t5d = Marker(eid='TORC.3181.T5D')
bhm_3181_t5d = Marker(eid='BHM.3181.T5D')
scrw_3186_t5d = Marker(eid='SCRW.3186.T5D')
duwindow_3186_t5d = Marker(eid='DUWINDOW.3186.T5D')
duflange_3186_t5d = Marker(eid='DUFLANGE.3186.T5D')
duconcrete_3187_t5d = Marker(eid='DUCONCRETE.3187.T5D')
duabsorb_3189_t5d = Marker(eid='DUABSORB.3189.T5D')
ensub_3189_t5d = Marker(eid='ENSUB.3189.T5D')
ensec_3189_t5d = Marker(eid='ENSEC.3189.T5D')

# Lattice 
cell = (stsec_3039_t5d, stsub_3039_t5d, d_2585, tora_3040_t5d, d_2586, qe_3052_t5d, d_2587, bpma_3068_t5d, d_2588, 
qf_3068_t5d, d_2589, cfy_3069_t5d, d_2590, bpma_3080_t5d, d_2591, qf_3081_t5d, d_2592, cfx_3081_t5d, d_2593, 
bpma_3141_t5d, d_2594, qf_3142_t5d, d_2595, cfx_3142_t5d, d_2596, cfy_3148_t5d, d_2597, qf_3148_t5d, d_2598, 
midbpmf_3149_t5d, d_2599, bpmf_3149_t5d, d_2600, tora_3149_t5d, d_2601, ensub_3149_t5d, stsub_3149_t5d, d_2602, bv_3151_t5d, 
d_2603, bv_3154_t5d, d_2604, cnx_3157_t5d, d_2605, qk_3158_t5d, d_2606, bpmd_3159_t5d, d_2607, sk_3159_t5d, 
d_2608, otrd_3160_t5d, d_2609, sk_3161_t5d, d_2610, bpmd_3162_t5d, d_2611, qk_3163_t5d, d_2612, cny_3164_t5d, 
d_2613, bv_3167_t5d, d_2614, bv_3170_t5d, d_2615, bpmd_3172_t5d, d_2616, qk_3172_t5d, d_2617, qk_3174_t5d, 
d_2618, qk_3175_t5d, d_2619, qk_3177_t5d, d_2620, bpmd_3177_t5d, d_2621, sweep_3178_t5d, d_2622, sweep_3179_t5d, 
d_2623, bpmd_3180_t5d, d_2624, otrd_3181_t5d, d_2625, torc_3181_t5d, d_2626, bhm_3181_t5d, d_2627, bpmw_3185_t5d, 
d_2628, scrw_3186_t5d, duwindow_3186_t5d, d_2629, duflange_3186_t5d, d_2630, duconcrete_3187_t5d, d_2631, duabsorb_3189_t5d, ensub_3189_t5d, 
ensec_3189_t5d)
