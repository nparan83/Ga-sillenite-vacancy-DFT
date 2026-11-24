#!/usr/bin/env python3
# ga_sekil_veri.py : SilenitB - sekil verilerini tek pakette topla (0 c-saat). Colab icin 99_analiz/sekil_veri_SilenitB.tar.gz
# Icerik: EOS .dat (gaB, gaS, ozet, si, + ga2A uretilir), PDOS (S, B), ga_tablolar.txt, *_oku.log, son yapilar XSF (Kartezyen A).
# Kullanim (dizin /arf/scratch/nparan/SilenitB): python3 99_analiz/ga_sekil_veri.py
import sys, os, re, shutil, tarfile, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
B, RY, B2A = g.B, 13.605693, 0.529177
O = B + '/99_analiz/sekil_veri'
ZORUNLU = ['99_analiz/eos_gaB.dat', '99_analiz/eos_gaS.dat', '99_analiz/eos_si.dat', '99_analiz/ga_tablolar.txt', '05_elek/pdos_gaS.dat', '05_elek/pdos_gaB.dat']
ISTEGE = ['99_analiz/eos_ozet.dat'] + ['99_analiz/%s_oku.log' % k for k in ('ga_f2', 'ga_f2a', 'ga_f2c', 'ga_f2cb', 'ga_f4si', 'ga_ref', 'ga_elek', 'ga_hse', 'ga_f4_ozet', 'ga_f4c')]
YAP = [('02_O40/gaS_rel', 'Bi25GaO40_S'), ('03_site/gaA_rel', 'Bi25GaO39_A'), ('03_site/gaB_rel', 'Bi25GaO39_B'), ('03_site/gaC_rel', 'Bi25GaO39_C'),
       ('03_site/ga2A_rel', 'Bi24Ga2O39_2A'), ('03_site/ga2B_rel', 'Bi24Ga2O39_2B'), ('03_site/ga2C_rel', 'Bi24Ga2O39_2C'), ('04_eos/si_s1.000', 'Bi12SiO20_s1.000')]

def konum(t):
    k = t.rfind('ATOMIC_POSITIONS'); X = []
    if k < 0 or 'crystal' not in t[k:k + 40].lower(): return X
    for l in t[k:].splitlines()[1:]:
        w = l.split()
        if len(w) < 4 or not re.fullmatch(r'[A-Z][a-z]?', w[0]): break
        X.append((w[0], [float(q) for q in w[1:4]]))
    return X

def xsf(d, ad):
    r = g.oku1('%s/%s/run.out' % (B, d))
    if r is None or not (r['jd'] and r['bf']): return '%s: JD/bfgs yok, atlandi' % d
    m = re.search(r'celldm\(1\)\s*=\s*([0-9.dDeE+-]+)', open('%s/%s/run.in' % (B, d)).read())
    if not re.search(r'ibrav\s*=\s*1\b', open('%s/%s/run.in' % (B, d)).read()) or not m: return '%s: ibrav=1/celldm yok, atlandi' % d
    a = float(m.group(1).lower().replace('d', 'e')) * B2A; X = konum(r['t'])
    with open('%s/%s.xsf' % (O, ad), 'w') as f:
        f.write('# %s/%s son konumlar; E = %.8f Ry\nCRYSTAL\nPRIMVEC\n' % (d, 'run.out', r['E']))
        for i in range(3): f.write(' '.join('%.8f' % (a if j == i else 0.) for j in range(3)) + '\n')
        f.write('PRIMCOORD\n%d 1\n' % len(X))
        for e, p in X: f.write('%-3s %12.6f %12.6f %12.6f\n' % (e, *(a * q for q in p)))
    return '%s -> %s.xsf (%d atom, a %.4f A)' % (d, ad, len(X), a)

def main():
    eks = [f for f in ZORUNLU if not os.path.exists('%s/%s' % (B, f))]
    if eks: g.dur('zorunlu dosya yok: ' + ', '.join(eks))
    if os.path.exists(O): shutil.rmtree(O)
    os.makedirs(O)
    for f in ZORUNLU + ISTEGE:
        if os.path.exists('%s/%s' % (B, f)): shutil.copy('%s/%s' % (B, f), O)
        else: print('istege bagli yok:', f)
    with open(O + '/eos_ga2A.dat', 'w') as o:
        o.write('# s V_A3 E_eV P_QE_GPa JDBF  (Bi24Ga2O39 2A, 90/360; 04_eos/ga2A_s*)\n')
        for s in (0.97, 0.98, 0.99, 1.00):
            r = g.oku1('%s/04_eos/ga2A_s%.3f/run.out' % (B, s))
            if r: o.write('%.3f %.4f %.6f %.4f %d\n' % (s, (19.58225862 * s * B2A) ** 3, r['E'] * RY, r['P'] / 10, r['jd'] and r['bf']))
    for d, ad in YAP: print(xsf(d, ad))
    tgz = B + '/99_analiz/sekil_veri_SilenitB.tar.gz'
    with tarfile.open(tgz, 'w:gz') as t: t.add(O, arcname='sekil_veri')
    print('\n'.join(sorted(os.listdir(O))))
    print('yazildi: %s (%.1f kB, md5 %s)' % (tgz, os.path.getsize(tgz) / 1024, hashlib.md5(open(tgz, 'rb').read()).hexdigest()))

if __name__ == '__main__':
    main()
