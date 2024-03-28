from ocelot import * 

#Initial Twiss parameters
tws0 = Twiss()
tws0.beta_x = 55.7981
tws0.beta_y = 55.7981
tws0.alpha_x = 18.1886
tws0.alpha_y = 18.1886
tws0.E = 0.005

# Drifts
d_0 = Drift(l=0.27599999999999997, eid='D_0')
d_1 = Drift(l=0.21600000000000008, eid='D_1')
d_2 = Drift(l=0.1080000000000001, eid='D_2')
d_3 = Drift(l=0.0030000000000000027, eid='D_3')
d_4 = Drift(l=1.1102230246251565e-16, eid='D_4')
d_5 = Drift(l=0.04699999999999993, eid='D_5')
d_6 = Drift(l=0.1499999999999999, eid='D_6')
d_7 = Drift(l=0.19878200000000001, eid='D_7')
d_8 = Drift(l=0.4392180000000001, eid='D_8')
d_9 = Drift(l=0.2250000000000001, eid='D_9')
d_10 = Drift(l=0.08799999999999963, eid='D_10')
d_11 = Drift(l=0.11100000000000021, eid='D_11')
d_12 = Drift(l=0.3099999999999996, eid='D_12')
d_13 = Drift(l=0.12639999999999985, eid='D_13')
ctbi_26_i1 = Drift(l=0.5512, eid='CTBI.26.I1')
d_14 = Drift(l=0.22160000000000002, eid='D_14')
d_15 = Drift(l=0.3459000000000003, eid='D_15')
d_16 = Drift(l=0.34589999999999854, eid='D_16')
d_17 = Drift(l=0.34589999999999854, eid='D_17')
d_18 = Drift(l=0.3459000000000003, eid='D_18')
d_19 = Drift(l=0.34589999999999854, eid='D_19')
d_20 = Drift(l=0.34589999999999854, eid='D_20')
d_21 = Drift(l=0.3459000000000003, eid='D_21')
d_22 = Drift(l=0.20429999999999815, eid='D_22')
d_23 = Drift(l=0.04320000000000057, eid='D_23')
d_24 = Drift(l=0.04320000000000057, eid='D_24')
d_25 = Drift(l=0.08499999999999908, eid='D_25')
d_26 = Drift(l=0.4579000000000004, eid='D_26')
d_27 = Drift(l=0.22109999999999985, eid='D_27')
d_28 = Drift(l=0.12820000000000142, eid='D_28')
d_29 = Drift(l=0.043199999999998795, eid='D_29')
d_30 = Drift(l=0.2020000000000035, eid='D_30')
d_31 = Drift(l=0.262000000000004, eid='D_31')
d_32 = Drift(l=0.26200000000000045, eid='D_32')
d_33 = Drift(l=0.262000000000004, eid='D_33')
d_34 = Drift(l=0.262000000000004, eid='D_34')
d_35 = Drift(l=0.262000000000004, eid='D_35')
d_36 = Drift(l=0.262000000000004, eid='D_36')
d_37 = Drift(l=0.262000000000004, eid='D_37')
d_38 = Drift(l=0.15500000000000114, eid='D_38')
cfbi_44_i1 = Drift(l=1.04, eid='CFBI.44.I1')
d_39 = Drift(l=1.2963999999999984, eid='D_39')
d_40 = Drift(l=0.3330500000000036, eid='D_40')
d_41 = Drift(l=0.08614999999999995, eid='D_41')
d_42 = Drift(l=0.12150000000000105, eid='D_42')
d_43 = Drift(l=0.09600000000000009, eid='D_43')
d_44 = Drift(l=0.21369999999999933, eid='D_44')
d_45 = Drift(l=0.12345000000000184, eid='D_45')
d_46 = Drift(l=0.08115000000000094, eid='D_46')
d_47 = Drift(l=0.051165000000001015, eid='D_47')
d_48 = Drift(l=0.10082699999999534, eid='D_48')
d_49 = Drift(l=0.1884149999999991, eid='D_49')
d_50 = Drift(l=0.09600000000000009, eid='D_50')
d_51 = Drift(l=0.17575000000000074, eid='D_51')
d_52 = Drift(l=0.19399999999999906, eid='D_52')
d_53 = Drift(l=0.19399999999999906, eid='D_53')
d_54 = Drift(l=0.12016499999999652, eid='D_54')
d_55 = Drift(l=0.10082799999999992, eid='D_55')
d_56 = Drift(l=0.05116499999999746, eid='D_56')
d_57 = Drift(l=0.08115000000000094, eid='D_57')
d_58 = Drift(l=0.1896400000000007, eid='D_58')
d_59 = Drift(l=0.1999999999999993, eid='D_59')
d_60 = Drift(l=0.3030000000000008, eid='D_60')
d_61 = Drift(l=0.10149999999999793, eid='D_61')
d_62 = Drift(l=0.09600000000000009, eid='D_62')
d_63 = Drift(l=0.1836500000000001, eid='D_63')
d_64 = Drift(l=3.552713678800501e-15, eid='D_64')
d_65 = Drift(l=0.05114999999999981, eid='D_65')
d_66 = Drift(l=0.05114999999999981, eid='D_66')
d_67 = Drift(l=0.27815000000000367, eid='D_67')
d_68 = Drift(l=0.0861500000000035, eid='D_68')
d_69 = Drift(l=0.0861500000000035, eid='D_69')
d_70 = Drift(l=0.28500000000000014, eid='D_70')
d_71 = Drift(l=0.03499999999999659, eid='D_71')
d_72 = Drift(l=0.17500000000000426, eid='D_72')
d_73 = Drift(l=0.15000000000000568, eid='D_73')
d_74 = Drift(l=0.1311499999999981, eid='D_74')
d_75 = Drift(l=0.0861499999999964, eid='D_75')
d_76 = Drift(l=0.3200000000000074, eid='D_76')
d_77 = Drift(l=0.17499999999999716, eid='D_77')
d_78 = Drift(l=0.15000000000000568, eid='D_78')
d_79 = Drift(l=0.1311499999999981, eid='D_79')
d_80 = Drift(l=0.0861499999999964, eid='D_80')
d_81 = Drift(l=7.105427357601002e-15, eid='D_81')
d_82 = Drift(l=0.3200000000000003, eid='D_82')
d_83 = Drift(l=0.2750000000000057, eid='D_83')
d_84 = Drift(l=0.18114999999999526, eid='D_84')
d_85 = Drift(l=0.20114999999999839, eid='D_85')
d_86 = Drift(l=0.38000000000000256, eid='D_86')
d_87 = Drift(l=0.2811499999999967, eid='D_87')
d_88 = Drift(l=0.26285000000000025, eid='D_88')
d_89 = Drift(l=0.36430000000000007, eid='D_89')
d_90 = Drift(l=0.14999999999999858, eid='D_90')
d_91 = Drift(l=0.2351499999999973, eid='D_91')
d_92 = Drift(l=0.43115, eid='D_92')
d_93 = Drift(l=0.09999999999999432, eid='D_93')

# Quadrupoles
qln_23_i1 = Quadrupole(l=0.05, eid='QLN.23.I1')
qls_23_i1 = Quadrupole(l=0.05, tilt=0.785398163, eid='QLS.23.I1')
q_37_i1 = Quadrupole(l=0.2136, k1=-1.3023162, eid='Q.37.I1')
q_38_i1 = Quadrupole(l=0.2136, k1=1.3324576310018725, eid='Q.38.I1')
qi_46_i1 = Quadrupole(l=0.2377, k1=0.3065153557004628, eid='QI.46.I1')
qi_47_i1 = Quadrupole(l=0.2377, k1=0.16454147949936895, eid='QI.47.I1')
qi_50_i1 = Quadrupole(l=0.2377, k1=-0.8646623293984014, eid='QI.50.I1')
qi_52_i1 = Quadrupole(l=0.2377, k1=-0.3522076140008414, eid='QI.52.I1')
qi_53_i1 = Quadrupole(l=0.2377, k1=2.104794185999159, eid='QI.53.I1')
qi_54_i1 = Quadrupole(l=0.2377, k1=0.7943661063020614, eid='QI.54.I1')
qi_55_i1 = Quadrupole(l=0.2377, k1=-2.6383360780016827, eid='QI.55.I1')
qi_57_i1 = Quadrupole(l=0.2377, k1=3.2780620240008416, eid='QI.57.I1')
qi_59_i1 = Quadrupole(l=0.2377, k1=-2.6383360780016827, eid='QI.59.I1')
qi_60_i1 = Quadrupole(l=0.2377, k1=1.9778194619983174, eid='QI.60.I1')
qi_61_i1 = Quadrupole(l=0.2377, k1=0.8708619870004207, eid='QI.61.I1')

# SBends
bl_48i_i1 = SBend(l=0.2, angle=-0.099484, e2=-0.099484, eid='BL.48I.I1')
bl_48ii_i1 = SBend(l=0.2, angle=0.099484, e1=0.099484, eid='BL.48II.I1')
bl_50i_i1 = SBend(l=0.2, angle=0.099484, e2=0.099484, eid='BL.50I.I1')
bl_50ii_i1 = SBend(l=0.2, angle=-0.099484, e1=-0.099484, eid='BL.50II.I1')

# Hcors
ckx_23_i1 = Hcor(l=0.025, eid='CKX.23.I1')
kix_24_i1 = Hcor(l=0.1, eid='KIX.24.I1')
ckx_24_i1 = Hcor(l=0.025, eid='CKX.24.I1')
ckx_25_i1 = Hcor(l=0.025, eid='CKX.25.I1')
cx_37_i1 = Hcor(eid='CX.37.I1')
cx_39_i1 = Hcor(eid='CX.39.I1')
cix_51_i1 = Hcor(l=0.1, eid='CIX.51.I1')
kjx_54_i1 = Hcor(l=0.175, eid='KJX.54.I1')
kax_55_i1 = Hcor(l=0.35, eid='KAX.55.I1')
kax_56_i1 = Hcor(l=0.35, eid='KAX.56.I1')
cix_57_i1 = Hcor(l=0.1, eid='CIX.57.I1')
kjx_57_i1 = Hcor(l=0.175, eid='KJX.57.I1')

# Vcors
cky_23_i1 = Vcor(l=0.025, eid='CKY.23.I1')
kiy_24_i1 = Vcor(l=0.1, eid='KIY.24.I1')
cky_24_i1 = Vcor(l=0.025, eid='CKY.24.I1')
cky_25_i1 = Vcor(l=0.025, eid='CKY.25.I1')
cy_37_i1 = Vcor(eid='CY.37.I1')
cy_39_i1 = Vcor(eid='CY.39.I1')
ciy_51_i1 = Vcor(l=0.1, eid='CIY.51.I1')
kjy_54_i1 = Vcor(l=0.175, eid='KJY.54.I1')
ciy_55_i1 = Vcor(l=0.1, eid='CIY.55.I1')
kjy_57_i1 = Vcor(l=0.175, eid='KJY.57.I1')
ciy_58_i1 = Vcor(l=0.1, eid='CIY.58.I1')

# Undulators
u74_49_i1 = Undulator(lperiod=0.074, nperiods=10, eid='U74.49.I1')

# Cavitys
c_a1_1_1_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, eid='C.A1.1.1.I1')
c_a1_1_2_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, eid='C.A1.1.2.I1')
c_a1_1_3_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, eid='C.A1.1.3.I1')
c_a1_1_4_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, eid='C.A1.1.4.I1')
c_a1_1_5_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, eid='C.A1.1.5.I1')
c_a1_1_6_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, eid='C.A1.1.6.I1')
c_a1_1_7_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, eid='C.A1.1.7.I1')
c_a1_1_8_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, eid='C.A1.1.8.I1')
c3_ah1_1_1_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, eid='C3.AH1.1.1.I1')
c3_ah1_1_2_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, eid='C3.AH1.1.2.I1')
c3_ah1_1_3_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, eid='C3.AH1.1.3.I1')
c3_ah1_1_4_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, eid='C3.AH1.1.4.I1')
c3_ah1_1_5_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, eid='C3.AH1.1.5.I1')
c3_ah1_1_6_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, eid='C3.AH1.1.6.I1')
c3_ah1_1_7_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, eid='C3.AH1.1.7.I1')
c3_ah1_1_8_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, eid='C3.AH1.1.8.I1')

# TDCavitys
tdsa_52_i1 = TDCavity(l=0.7, freq=2800000000.0, tilt=1.570796327, eid='TDSA.52.I1')

# Monitors
bpmg_24_i1 = Monitor(eid='BPMG.24.I1')
bpmg_25i_i1 = Monitor(eid='BPMG.25I.I1')
bpmc_38i_i1 = Monitor(eid='BPMC.38I.I1')
bpmr_38ii_i1 = Monitor(eid='BPMR.38II.I1')
bpmf_47_i1 = Monitor(eid='BPMF.47.I1')
bpmf_48_i1 = Monitor(eid='BPMF.48.I1')
bpmf_52_i1 = Monitor(eid='BPMF.52.I1')
bpma_55_i1 = Monitor(eid='BPMA.55.I1')
bpma_57_i1 = Monitor(eid='BPMA.57.I1')
bpma_59_i1 = Monitor(eid='BPMA.59.I1')
bpma_60_i1 = Monitor(eid='BPMA.60.I1')

# Markers
stsec_23_i1 = Marker(eid='STSEC.23.I1')
stsub_23_i1 = Marker(eid='STSUB.23.I1')
gun_23_i1 = Marker(eid='GUN.23.I1')
scrn_24_i1 = Marker(eid='SCRN.24.I1')
fcup_24_i1 = Marker(eid='FCUP.24.I1')
ensub_24_i1 = Marker(eid='ENSUB.24.I1')
stsub_24i_i1 = Marker(eid='STSUB.24I.I1')
tora_25_i1 = Marker(eid='TORA.25.I1')
scrn_25i_i1 = Marker(eid='SCRN.25I.I1')
fcup_25i_i1 = Marker(eid='FCUP.25I.I1')
dcm_25_i1 = Marker(eid='DCM.25.I1')
stac_26_i1 = Marker(eid='STAC.26.I1')
a1_1_stop = Marker(eid='a1_1_stop')
a1_sim_stop = Marker(eid='a1_sim_stop')
ah1_sim_start = Marker(eid='ah1_sim_start')
match_37_i1 = Marker(eid='MATCH.37.I1')
enac_38_i1 = Marker(eid='ENAC.38.I1')
stac_38_i1 = Marker(eid='STAC.38.I1')
enac_44_i1 = Marker(eid='ENAC.44.I1')
tora_46_i1 = Marker(eid='TORA.46.I1')
bam_47_i1 = Marker(eid='BAM.47.I1')
midbpmf_47_i1 = Marker(eid='MIDBPMF.47.I1')
dcm_47_i1 = Marker(eid='DCM.47.I1')
stlat_47_i1 = Marker(eid='STLAT.47.I1')
midbpmf_48_i1 = Marker(eid='MIDBPMF.48.I1')
otrl_48_i1 = Marker(eid='OTRL.48.I1')
lh_start = Marker(eid='lh_start')
lh_stop = Marker(eid='lh_stop')
otrl_50_i1 = Marker(eid='OTRL.50.I1')
enlat_50_i1 = Marker(eid='ENLAT.50.I1')
eod_51_i1 = Marker(eid='EOD.51.I1')
midbpmf_52_i1 = Marker(eid='MIDBPMF.52.I1')
match_52_i1 = Marker(eid='MATCH.52.I1')
stlat_55_i1 = Marker(eid='STLAT.55.I1')
match_55_i1 = Marker(eid='MATCH.55.I1')
otrc_55_i1 = Marker(eid='OTRC.55.I1')
otrc_56_i1 = Marker(eid='OTRC.56.I1')
otrc_58_i1 = Marker(eid='OTRC.58.I1')
otrc_59_i1 = Marker(eid='OTRC.59.I1')
tora_60_i1 = Marker(eid='TORA.60.I1')
bpmatest_61_i1 = Marker(eid='BPMATEST.61.I1')
ensub_62_i1 = Marker(eid='ENSUB.62.I1')
dump_csr_start = Marker(eid='DUMP.CSR.START')

# Lattice 
cell = (stsec_23_i1, stsub_23_i1, gun_23_i1, d_0, qln_23_i1, qls_23_i1, d_1, ckx_23_i1, cky_23_i1, 
d_2, kix_24_i1, kiy_24_i1, d_3, ckx_24_i1, d_4, cky_24_i1, d_5, bpmg_24_i1, d_6, 
scrn_24_i1, fcup_24_i1, d_7, ensub_24_i1, stsub_24i_i1, d_8, tora_25_i1, d_9, scrn_25i_i1, fcup_25i_i1, 
d_10, bpmg_25i_i1, d_11, dcm_25_i1, d_12, ckx_25_i1, cky_25_i1, d_13, ctbi_26_i1, stac_26_i1, 
d_14, c_a1_1_1_i1, a1_1_stop, d_15, c_a1_1_2_i1, d_16, c_a1_1_3_i1, d_17, c_a1_1_4_i1, d_18, 
c_a1_1_5_i1, d_19, c_a1_1_6_i1, d_20, c_a1_1_7_i1, d_21, c_a1_1_8_i1, a1_sim_stop, ah1_sim_start, d_22, 
match_37_i1, d_23, q_37_i1, d_24, cx_37_i1, cy_37_i1, d_25, bpmc_38i_i1, d_26, enac_38_i1, 
stac_38_i1, d_27, bpmr_38ii_i1, d_28, q_38_i1, d_29, cx_39_i1, cy_39_i1, d_30, c3_ah1_1_1_i1, 
d_31, c3_ah1_1_2_i1, d_32, c3_ah1_1_3_i1, d_33, c3_ah1_1_4_i1, d_34, c3_ah1_1_5_i1, d_35, c3_ah1_1_6_i1, 
d_36, c3_ah1_1_7_i1, d_37, c3_ah1_1_8_i1, d_38, enac_44_i1, cfbi_44_i1, d_39, tora_46_i1, d_40, 
qi_46_i1, d_41, bam_47_i1, d_42, bpmf_47_i1, d_43, midbpmf_47_i1, d_44, dcm_47_i1, d_45, 
qi_47_i1, d_46, stlat_47_i1, d_47, bl_48i_i1, d_48, bl_48ii_i1, d_49, midbpmf_48_i1, d_50, 
bpmf_48_i1, d_51, otrl_48_i1, d_52, lh_start, u74_49_i1, lh_stop, d_53, otrl_50_i1, d_54, 
bl_50i_i1, d_55, bl_50ii_i1, d_56, enlat_50_i1, d_57, qi_50_i1, d_58, eod_51_i1, d_59, 
ciy_51_i1, d_60, cix_51_i1, d_61, bpmf_52_i1, d_62, midbpmf_52_i1, d_63, match_52_i1, d_64, 
qi_52_i1, d_65, tdsa_52_i1, d_66, qi_53_i1, d_67, kjx_54_i1, kjy_54_i1, d_68, qi_54_i1, 
d_69, kax_55_i1, d_70, stlat_55_i1, d_71, match_55_i1, otrc_55_i1, d_72, ciy_55_i1, d_73, 
bpma_55_i1, d_74, qi_55_i1, d_75, kax_56_i1, d_76, otrc_56_i1, d_77, cix_57_i1, d_78, 
bpma_57_i1, d_79, qi_57_i1, d_80, kjx_57_i1, d_81, kjy_57_i1, d_82, otrc_58_i1, d_83, 
ciy_58_i1, d_84, qi_59_i1, d_85, bpma_59_i1, d_86, otrc_59_i1, d_87, qi_60_i1, d_88, 
tora_60_i1, d_89, bpma_60_i1, d_90, bpmatest_61_i1, d_91, qi_61_i1, ensub_62_i1, d_92, dump_csr_start, 
d_93)
