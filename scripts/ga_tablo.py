#!/usr/bin/env python3
# ga_tablo.py : SilenitB - makale sayilarinin TEK kaynagi (0 c-saat). Tum enerjiler run.out/XML'den okunur; elle sayi yok.
# T1 Bi25GaO39 siteleri (E, E_form O-rich/O-poor, geometri) | T2 referanslar + hull | T3 Bi24Ga2O39 | T4 gap (PBE/PAW/HSE)
# Kullanim (dizin /arf/scratch/nparan/SilenitB): python3 99_analiz/ga_tablo.py | tee 99_analiz/ga_tablolar.txt
import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
import ga_f2 as g2
import ga_elek as ge
B, F, RY = g.B, g.F, 13.605693
A = 19.58225862 * 0.529177; SC = 10.183 / 10.3625

def E(path, need_bf=False):
    r = g.oku1(path)
    if r is None or not r['jd'] or (need_bf and not r['bf']): g.dur('eksik/yakinsamamis: ' + path)
    return r

def ref(name):                                  # Fe projesi 06_referans (SG15, 80/320; teyit 12.11 ve 17.11)
    t = open('%s/06_referans/%s.out' % (F, name), errors='ignore').read()
    return float(re.findall(r'^!\s+total energy\s+=\s+(-?[\d.]+)', t, re.M)[-1])

S40 = E(B + '/02_O40/gaS_rel/run.out', True)
eO2, eBi, eBi2O3 = ref('O2'), ref('Bi_bulk') / 2, ref('Bi2O3') / 4
eGa = E(B + '/06_ref/aGa_scf/run.out')['E'] / 8; eGa2O3 = E(B + '/06_ref/bGa2O3_scf/run.out')['E'] / 4
muR = eO2 / 2; muP = (S40['E'] - 25 * eBi - eGa) / 40
print('# kaynak: %s (+ Fe 06_referans: O2, Bi_bulk, Bi2O3)' % B)
print('mu_O O-rich %.8f Ry | O-poor (Bi+Ga metal sinir) %.8f Ry | aralik %.3f eV' % (muR, muP, (muR - muP) * RY))
print('\nT1  Bi25GaO39 O bosluk siteleri (80/320, k2, SG15 relaks)')
print('%-5s %10s %9s %9s %9s %9s %18s' % ('site', 'dE (eV)', 'Ef rich', 'Ef poor', 'Bt||111', 'Ga||111', 'Ga-O<2.3 (olcekli)'))
R = {s: E('%s/03_site/ga%s_rel/run.out' % (B, s), True) for s in 'ABC'}
for s in 'BAC':
    q = g2.geo(g2.konum(R[s]['t']), A); b = [d for d in q['gaO'] if d < 2.3]; ef = (R[s]['E'] - S40['E'] + muR) * RY
    print('%-5s %+10.4f %+9.3f %+9.3f %+9.3f %+9.3f %6d x %.3f' % (s, (R[s]['E'] - R['B']['E']) * RY, ef, ef - (muR - muP) * RY, q['kay'], q['gak'], len(b), SC * sum(b) / len(b)))
print('\nT2  referanslar (PBE, 80/320)')
dHBi = (eBi2O3 - 2 * eBi - 1.5 * eO2) * RY; dHGa = (eGa2O3 - 2 * eGa - 1.5 * eO2) * RY
h40 = (S40['E'] - 12.5 * eBi2O3 - 0.5 * eGa2O3 - 0.5 * eO2) * RY; h39 = (R['B']['E'] - 12.5 * eBi2O3 - 0.5 * eGa2O3) * RY
print('dH_f Bi2O3 %.3f eV/f.u. | dH_f beta-Ga2O3 %.3f eV/f.u.' % (dHBi, dHGa))
print('hull: Bi25GaO40 %+.3f eV/hucre (%+.4f eV/atom) | Bi25GaO39 %+.3f eV/hucre (%+.4f eV/atom)' % (h40, h40 / 66, h39, h39 / 65))
print('\nT3  Bi24Ga2O39 (Bi_tet -> Ga)')
R2 = {s: E('%s/03_site/ga2%s_rel/run.out' % (B, s), True) for s in 'ABC'}; k0 = min(R2, key=lambda s: R2[s]['E'])
for s in 'ABC': print('ga2%s  dE %+.4f eV  P %.2f kbar' % (s, (R2[s]['E'] - R2[k0]['E']) * RY, R2[s]['P']))
dEd = (R2[k0]['E'] - 24 / 25 * R['B']['E'] - 13 / 25 * eGa2O3) * RY; dP = (R2[k0]['P'] - R['B']['P']) / 10
print('ayrisma -> 24/25 Bi25GaO39 + 13/25 Ga2O3: %+.3f eV/f.u. (%+.2f meV/atom) | da (birinci mertebe, B0 45.6 GPa) %+.4f A olcekli' % (dEd, dEd / 65 * 1000, A * dP / (3 * 45.6) * SC))
print('\nT4  bant araliklari (eV)')
for y in 'SB':
    n = ge.xml_bant('%s/ga%s_nscf/out/%s.save/data-file-schema.xml' % (ge.E5, y, ge.YAP[y][1]))['gap']
    p = ge.xml_bant('%s/ga%s_paw/out/ga%s_paw.save/data-file-schema.xml' % (ge.E5, y, y))['gap']
    h = []
    for d in (['ga%s_hse' % y] + (['gaS_f320_hse'] if y == 'S' else [])):
        r = E('%s/07_hse/%s/run.out' % (B, d)); L = re.findall(r'highest occupied, lowest unoccupied level \(ev\):\s*(-?[\d.]+)\s+(-?[\d.]+)', r['t'])
        h.append('%s %.3f' % (d, float(L[-1][1]) - float(L[-1][0])))
    print('%s (%s): PBE SG15 k4 %.3f | PAW k2 %.3f | HSE06 k2: %s' % (y, 'Bi25GaO40' if y == 'S' else 'Bi25GaO39', n, p, ' | '.join(h)))
print('deney (durum 2b): 2.68 (Wang 2012), 2.95 (Lin 2008), 2.99 (Kargin derlemesi)')
