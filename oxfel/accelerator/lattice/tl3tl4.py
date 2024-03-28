from ocelot import * 

# Drifts
d_1939 = Drift(l=1.4724009999999907, eid='D_1939')
d_1940 = Drift(l=13.047401000000036, eid='D_1940')
d_1941 = Drift(l=1.1589500000002317, eid='D_1941')
d_1942 = Drift(l=0.1539500000001226, eid='D_1942')
d_1943 = Drift(l=0.16000000000008185, eid='D_1943')
d_1944 = Drift(l=0.07000000000016371, eid='D_1944')
d_1945 = Drift(l=0.07000000000016371, eid='D_1945')
d_1946 = Drift(l=0.2800000000002001, eid='D_1946')
d_1947 = Drift(l=0.07000000000016371, eid='D_1947')
d_1948 = Drift(l=0.07000000000016371, eid='D_1948')
d_1949 = Drift(l=0.20500000000015461, eid='D_1949')
d_1950 = Drift(l=0.19000000000005457, eid='D_1950')
d_1951 = Drift(l=5.3800000000003365, eid='D_1951')
d_1952 = Drift(l=0.445699999999988, eid='D_1952')
d_1953 = Drift(l=0.15429999999992106, eid='D_1953')
d_1954 = Drift(l=0.20895100000007005, eid='D_1954')
d_1955 = Drift(l=0.15395100000023376, eid='D_1955')
d_1956 = Drift(l=0.20000000000004547, eid='D_1956')
d_1957 = Drift(l=7.784700000000157, eid='D_1957')
d_1958 = Drift(l=0.09999999999968168, eid='D_1958')
d_1959 = Drift(l=1.4347000000000207, eid='D_1959')
d_1960 = Drift(l=0.1999999999998181, eid='D_1960')
d_1961 = Drift(l=0.15000000000009095, eid='D_1961')
vcv100_2023_tl = Drift(l=0.085, eid='VCV100.2023.TL')
d_1962 = Drift(l=0.585700000000088, eid='D_1962')
vcb100_2024_tl = Drift(l=0.18, eid='VCB100.2024.TL')
vcabsa_2024_tl = Drift(l=0.61, eid='VCABSA.2024.TL')
d_1963 = Drift(l=0.2500199999999495, eid='D_1963')

# Quadrupoles
qk_1982_tl = Quadrupole(l=1.0552, k1=0.09035960007960576, eid='QK.1982.TL')
qf_1997_tl = Quadrupole(l=0.5321, k1=-0.17919084760007517, eid='QF.1997.TL')
qf_2012_tl = Quadrupole(l=0.5321, k1=0.17919084760007517, eid='QF.2012.TL')

# RBends
kl_1998_tl = RBend(l=0.93, tilt=1.570796327, eid='KL.1998.TL')
kl_1999_tl = RBend(l=0.93, tilt=1.570796327, eid='KL.1999.TL')
kl_2000_tl = RBend(l=0.93, tilt=1.570796327, eid='KL.2000.TL')
kl_2001_tl = RBend(l=0.93, tilt=1.570796327, eid='KL.2001.TL')
kl_2002_tl = RBend(l=0.93, tilt=1.570796327, eid='KL.2002.TL')
kl_2003_tl = RBend(l=0.93, tilt=1.570796327, eid='KL.2003.TL')
kl_2005_tl = RBend(l=0.93, eid='KL.2005.TL')

# Hcors
chx_2012_tl = Hcor(l=0.2, eid='CHX.2012.TL')
cnx_2021_tl = Hcor(l=0.3, eid='CNX.2021.TL')

# Vcors
chy_1997_tl = Vcor(l=0.2, eid='CHY.1997.TL')
chy_2004_tl = Vcor(l=0.2, eid='CHY.2004.TL')
cfy_2010_tl = Vcor(l=0.1, eid='CFY.2010.TL')
chy_2012_tl = Vcor(l=0.2, eid='CHY.2012.TL')
cny_2021_tl = Vcor(l=0.3, eid='CNY.2021.TL')

# Monitors
bpma_1995_tl = Monitor(eid='BPMA.1995.TL')
bpma_2011_tl = Monitor(eid='BPMA.2011.TL')
bpmd_2022_tl = Monitor(eid='BPMD.2022.TL')

# Markers
stsub_1980_tl = Marker(eid='STSUB.1980.TL')
ensub_1997_tl = Marker(eid='ENSUB.1997.TL')
stsub_1997_tl = Marker(eid='STSUB.1997.TL')
tora_2011_tl = Marker(eid='TORA.2011.TL')
otre_2023_tl = Marker(eid='OTRE.2023.TL')
ensub_2025_tl = Marker(eid='ENSUB.2025.TL')

# Lattice 
cell = (stsub_1980_tl, d_1939, qk_1982_tl, d_1940, bpma_1995_tl, d_1941, qf_1997_tl, d_1942, chy_1997_tl, 
ensub_1997_tl, stsub_1997_tl, d_1943, kl_1998_tl, d_1944, kl_1999_tl, d_1945, kl_2000_tl, d_1946, kl_2001_tl, 
d_1947, kl_2002_tl, d_1948, kl_2003_tl, d_1949, chy_2004_tl, d_1950, kl_2005_tl, d_1951, cfy_2010_tl, 
d_1952, tora_2011_tl, d_1953, bpma_2011_tl, d_1954, qf_2012_tl, d_1955, chx_2012_tl, d_1956, chy_2012_tl, 
d_1957, cnx_2021_tl, d_1958, cny_2021_tl, d_1959, bpmd_2022_tl, d_1960, otre_2023_tl, d_1961, vcv100_2023_tl, 
d_1962, vcb100_2024_tl, vcabsa_2024_tl, d_1963, ensub_2025_tl)
