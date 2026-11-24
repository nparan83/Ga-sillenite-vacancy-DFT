#!/usr/bin/env python3
# ga_f4si.py : SilenitB F4-Si - yontem dogrulamasi: Bi12SiO20 (Bi24Si2O40, 66 atom) EOS, Ga silenit F4 ile AYNI yontem
# (SG15 ONCV, 90/360 Ry, k 2x2x2, nspin=1, sabit hacim relax, BM3 E(V) fiti; sablon = 02_O40/gaS_rel/run.in, yalniz Ga -> Si ve konumlar).
# Yapi: COD 1533225 (Neov, Reehuis, Marinova, Sonntag, Appl. Phys. A 74 (2002) S1016; notron, a = 10.104 A) -> 06_ref/cif/1533225.cif
# Deney: a = 10.105(1) A; B0 = 63(1) GPa (RUS, Haussuhl ve ark. Mater. Res. Express 7 (2020) 025701, Tablo 4; XRD 63(4) Wiehl ve ark. 2010)
# Kullanim (dizin /arf/scratch/nparan/SilenitB): python3 99_analiz/ga_f4si.py uret | oku     Is: sbatch -J eSi0.98 job_eos.sh si_s0.980
import sys, os, re, glob, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
import ga_f4 as g4
import ga_ref as gr
import ga_f2cb as gcb
B, F, RY, B2A = g.B, g.F, 13.605693, 0.529177
CIF = B + '/06_ref/cif/1533225.cif'
AC = 10.105 * 10.3567 / 10.183                    # izgara merkezi (A): deney a x PBE/deney orani (Bi25GaO39); sonucu etkilemez
OL = (0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02)
AEXP, BEXP = 10.105, 63.0
GA = dict(a_pbe=10.3567, a_exp=10.183, B0=45.6)  # Bi25GaO39: F2c-b P(V) a0; Radaev 1989; F4 BM3 sikistirma kolu
def ad(s): return 'si_s%.3f' % s

def mesafe(X, a, el1, el2, n):                    # her el1 icin en yakin n el2 mesafesi (kubik, en yakin goruntu)
    out = []
    for e, p in X:
        if e != el1: continue
        d = sorted(a * sum(min(abs(u - v), 1 - abs(u - v)) ** 2 for u, v in zip(p, q)) ** .5 for e2, q in X if e2 == el2 and q is not p)
        out.append(d[:n])
    return out

def konum(t):                                     # son ATOMIC_POSITIONS blogu (Bi/Si/O)
    k = t.rfind('ATOMIC_POSITIONS'); X = []
    for l in (t[k:].splitlines()[1:] if k >= 0 else []):
        w = l.split()
        if len(w) < 4 or w[0] not in ('Bi', 'Si', 'O'): break
        X.append((w[0], [float(q) for q in w[1:4]]))
    return X

def yapi():
    if not os.path.exists(CIF): g.dur('CIF yok: %s (COD 1533225 indirilecek, SilenitB_kodlari 18)' % CIF)
    cel, X = gr.cif(CIF)
    if max(abs(c - 90) for c in cel[3:]) > 1e-3 or max(cel[:3]) - min(cel[:3]) > 1e-4 or not 10.0 < cel[0] < 10.2: g.dur('CIF hucresi kubik/10.1 A degil: %s' % cel)
    say = lambda Y: tuple(sum(1 for e, _ in Y if e == k) for k in ('Bi', 'Si', 'O'))
    if say(X) == (12, 1, 20):                     # merkezleme islemleri CIF'te yoksa I ekle
        X = X + [(e, [(u + .5) % 1 for u in p]) for e, p in X]
    if say(X) != (24, 2, 40): g.dur('CIF atom sayimi %s (beklenen Bi24 Si2 O40)' % (say(X),))
    so = mesafe(X, cel[0], 'Si', 'O', 4); bo = mesafe(X, cel[0], 'Bi', 'O', 1)
    if not all(1.55 < d < 1.75 for r in so for d in r): g.dur('Si-O mesafeleri beklenmedik: %s' % so)
    if min(r[0] for r in bo) < 1.9: g.dur('Bi-O < 1.9 A')
    return cel[0], X, so, bo

def psi():
    h = sorted(glob.glob(B + '/00_pseudo/Si_ONCV_PBE-*.upf'))
    if not h:
        k = sorted(glob.glob(F + '/00_pseudo/Si_ONCV_PBE-*.upf') + glob.glob('/arf/scratch/nparan/*/00_pseudo/Si_ONCV_PBE-*.upf'))
        k = [x for x in k if x.endswith('-1.0.upf')] or k
        if not k: g.dur('SG15 Si pseudo (Si_ONCV_PBE-*.upf) bulunamadi -> 00_pseudo/ altina konacak (SilenitB_kodlari 18)')
        shutil.copy(k[-1], B + '/00_pseudo/'); print('kopyalandi:', k[-1]); h = [B + '/00_pseudo/' + os.path.basename(k[-1])]
    v = [x for x in h if x.endswith('-1.0.upf')]                 # Bi/O/Ga ile ayni SG15 surumu (1.0) tercih
    u = os.path.basename((v or h)[-1])
    if abs(g.zval(u) - 4) > 1e-6: g.dur(u + ': z_valence 4 degil')
    return u

def uret():
    if 'temizlik' not in open(B + '/job_eos.sh').read(): g.dur('job_eos.sh temizliksiz (SilenitB_kodlari 9)')
    a, X, so, bo = yapi(); u = psi()
    T = open(B + '/02_O40/gaS_rel/run.in').read()
    if not (re.search(r"calculation\s*=\s*'relax'", T) and re.search(r'ibrav\s*=\s*1\b', T) and re.search(r'nat\s*=\s*66\b', T)): g.dur('sablon relax/ibrav=1/nat=66 degil')
    T, n = re.subn(r'^\s*Ga\s+\S+\s+\S+\s*$', '  Si  28.0855  ' + u, T, flags=re.M)
    if n != 1: g.dur('sablonda tek Ga tur satiri yok')
    sp = dict(re.findall(r'^\s*(Bi|Si|O)\s+[0-9.]+\s+(\S+\.upf)\s*$', T, re.M | re.I))
    ne = 24 * g.zval(sp['Bi']) + 2 * g.zval(sp['Si']) + 40 * g.zval(sp['O'])
    if abs(ne - round(ne)) > 1e-6 or int(round(ne)) % 2: g.dur('elektron sayisi %.1f cift degil' % ne)
    L = T.splitlines(); i = [k for k, l in enumerate(L) if l.strip().upper().startswith('ATOMIC_POSITIONS')]
    if len(i) != 1 or 'crystal' not in L[i[0]].lower(): g.dur('ATOMIC_POSITIONS crystal tek degil')
    i = i[0]; j = i + 1
    while j < len(L) and L[j].split() and L[j].split()[0] in ('Bi', 'Ga', 'O'): j += 1
    if j - i - 1 != 66: g.dur('sablonda %d konum' % (j - i - 1))
    T = '\n'.join(L[:i + 1] + ['%-3s %14.10f %14.10f %14.10f' % (e, *p) for e, p in X] + L[j:]) + '\n'
    print('CIF a %.4f A | Si-O %s A | Bi-O min %.3f A | %s | nelec %d' % (a, ' '.join('%.3f' % d for d in so[0]), min(r[0] for r in bo), u, ne))
    for s in OL:
        t = T
        for k, v in (('celldm(1)', '%.8f' % (AC * s / B2A)), ('ecutwfc', '90.0'), ('ecutrho', '360.0'), ('prefix', "'%s'" % ad(s)), ('outdir', "'./out'"), ('tstress', '.true.')):
            t = g.setk(t, k, v)
        g.yaz('%s/04_eos/%s' % (B, ad(s)), t); print('%s: a %.4f A, V %.2f A3 | yazildi' % (ad(s), AC * s, (AC * s) ** 3))

def oku():
    R = []
    for s in OL:
        r = g.oku1('%s/04_eos/%s/run.out' % (B, ad(s)))
        if r is None: print(ad(s), 'YOK'); continue
        a = AC * s; X = konum(r['t'])
        so = mesafe(X, a, 'Si', 'O', 4) if len(X) == 66 else []
        print('%s: JD %d bfgs %d adim %d P %+.2f GPa | Si-O %s | %s' % (ad(s), r['jd'], r['bf'], r['ns'], r['P'] / 10, ' '.join('%.3f' % d for d in so[0]) if so else '-', r['w']))
        R.append((s, a ** 3, r['E'] * RY, r['P'] / 10, r['jd'] and r['bf']))
    G = [x for x in R if x[4]]
    if len(G) < 5: print('### S1: en az 5 gecerli nokta gerekir'); return
    f = g4.bm3([x[1] for x in G], [x[2] for x in G])
    if f is None: print('### BM3 minimum yok'); return
    dP = max(abs(f['P'](x[1]) - x[3]) for x in G); vp, rp = gcb.p0([(x[1], x[3]) for x in G])
    print('BM3 E(V): V0 %.2f A3, a0 %.4f A, B0 %.1f GPa, B0p %.2f | maks artik %.2f meV | maks |P_fit-P_QE| %.2f GPa | V0 aralikta %s' % (
        f['V0'], f['a0'], f['B0'], f['Bp'], f['res'], dP, min(x[1] for x in G) < f['V0'] < max(x[1] for x in G)))
    print('P(V)=0 (ikinci derece): a0 %.4f A (maks artik %.3f GPa)' % (vp ** (1 / 3), rp))
    da, dg = f['a0'] / AEXP - 1, GA['a_pbe'] / GA['a_exp'] - 1; dB = f['B0'] / BEXP - 1
    print('S4 kafes: a0/a_exp - 1 = %+.2f %% (Si) vs %+.2f %% (Ga) -> fark %.2f %% : %s' % (100 * da, 100 * dg, 100 * abs(da - dg), 'olcek 0.9827 TASINABILIR' if abs(da - dg) <= 0.003 else 'olcek tasinabilirligi ZAYIF'))
    print('S5 B0: PBE %.1f vs deney %.0f GPa -> %+.1f %% : %s | Ga B0 ham %.1f, olceklenmis tahmin %.1f GPa, oran B0(Ga)/B0(Si) = %.3f' % (
        f['B0'], BEXP, 100 * dB, '|dB| <= %10' if abs(dB) <= 0.10 else '|dB| > %10 -> Ga B0 duzeltme notuyla', GA['B0'], GA['B0'] * BEXP / f['B0'], GA['B0'] / f['B0']))
    with open(B + '/99_analiz/eos_si.dat', 'w') as o:
        o.write('# s V_A3 E_eV P_QE_GPa JDBF | BM3 V0 %.4f B0 %.3f B0p %.4f E0 %.6f\n' % (f['V0'], f['B0'], f['Bp'], f['E0']))
        for x in R: o.write('%.3f %.4f %.6f %.4f %d\n' % (x[0], x[1], x[2], x[3], x[4]))
    print('yazildi: 99_analiz/eos_si.dat')

if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] not in ('uret', 'oku'): g.dur('kullanim: ga_f4si.py uret | oku')
    uret() if sys.argv[1] == 'uret' else oku()
