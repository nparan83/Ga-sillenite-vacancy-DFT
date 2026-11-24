#!/usr/bin/env python3
# ga_f4.py : SilenitB F4 - Ga silenit EOS: sabit hacim relax, 90/360 Ry, k2, nspin=1; BM3 E(V) fiti (numpy GEREKMEZ)
# Kullanim (dizin /arf/scratch/nparan/SilenitB):  python3 99_analiz/ga_f4.py uret B [S]  |  python3 99_analiz/ga_f4.py oku B [S]
# B = Bi25GaO39 (gaB_rel son konumlari), S = Bi25GaO40 (gaS_rel son konumlari). Is: sbatch -J eB0.98 job_eos.sh gaB_s0.980
import sys, os, re, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
import ga_f2 as g2
B, RY, B2A, GPA = g.B, 13.605693, 0.529177, 160.21766
A0 = 19.58225862                                   # celldm(1), bohr (Fe hucresi, F1b ile teyitli)
OL = (0.94, 0.96, 0.98, 0.99, 1.00, 1.01, 1.02)    # s = a / a0
YAP = {'B': ('03_site/gaB_rel', 65), 'S': ('02_O40/gaS_rel', 66)}
JOB = g.HEAD.replace('@J@', 'gaB_eos').replace('@T@', '12:00:00') + """cd @B@/04_eos/$1 || exit 1
[ -s run.out ] && { echo "### run.out zaten var -> DUR"; exit 1; }
mpirun pw.x -npool 2 -in run.in > run.out
"""

def ad(y, s): return 'ga%s_s%.3f' % (y, s)

def kaynak(y):
    d, n = YAP[y]; r = g.oku1('%s/%s/run.out' % (B, d))
    if r is None or not (r['jd'] and r['bf']): g.dur('%s: run.out yok ya da JD/bfgs yok (K1)' % d)
    X = g2.konum(r['t'])
    if len(X) != n: g.dur('%s: %d atom okundu (beklenen %d)' % (d, len(X), n))
    return open('%s/%s/run.in' % (B, d)).read(), X

def uret(sec):
    for y in sec:
        T, X = kaynak(y)
        if not re.search(r"calculation\s*=\s*'relax'", T, re.I): g.dur(y + ": sablon 'relax' degil")
        if not re.search(r'ibrav\s*=\s*1\b', T): g.dur(y + ': ibrav=1 degil')
        m = re.search(r'celldm\(1\)\s*=\s*([0-9.dDeE+-]+)', T)
        if not m or abs(float(m.group(1).lower().replace('d', 'e')) - A0) > 1e-6: g.dur(y + ': celldm(1) != %.8f' % A0)
        L = T.splitlines()
        i = [k for k, l in enumerate(L) if l.strip().upper().startswith('ATOMIC_POSITIONS')]
        if len(i) != 1 or 'crystal' not in L[i[0]].lower(): g.dur(y + ': ATOMIC_POSITIONS (crystal) tek degil')
        i = i[0]; j = i + 1
        while j < len(L) and L[j].split() and L[j].split()[0] in ('Bi', 'Ga', 'O'): j += 1
        if j - i - 1 != len(X): g.dur('%s: sablonda %d konum, son konumlar %d' % (y, j - i - 1, len(X)))
        T2 = '\n'.join(L[:i + 1] + ['%-3s %14.10f %14.10f %14.10f' % (s, *x) for s, x in X] + L[j:]) + '\n'
        for s in OL:
            n, t = ad(y, s), T2
            for k, v in (('celldm(1)', '%.8f' % (A0 * s)), ('ecutwfc', '90.0'), ('ecutrho', '360.0'),
                         ('prefix', "'%s'" % n), ('outdir', "'./out'"), ('tstress', '.true.')):
                t = g.setk(t, k, v)
            g.yaz('%s/04_eos/%s' % (B, n), t)
            print('%-12s a %.5f A  V %8.2f A3  nat %d  90/360 relax' % (n, A0 * s * B2A, (A0 * s * B2A) ** 3, len(X)))
    open(B + '/job_eos.sh', 'w').write(JOB.replace('@B@', B))
    print('yazildi: 04_eos/, job_eos.sh')

def coz(M, b):                                    # Gauss eleme, kismi pivot
    n = len(b); A = [M[r][:] + [b[r]] for r in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(A[r][c])); A[c], A[p] = A[p], A[c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            for k in range(c, n + 1): A[r][k] -= f * A[c][k]
    x = [0.] * n
    for r in range(n - 1, -1, -1): x[r] = (A[r][n] - sum(A[r][k] * x[k] for k in range(r + 1, n))) / A[r][r]
    return x

def bm3(V, E):                                    # E = sum c_k t^k, t = (V/Vm)^(-2/3) - 1  (BM3 ile esdeger)
    Vm = sum(V) / len(V); t = [(v / Vm) ** (-2. / 3) - 1 for v in V]
    c = coz([[sum(q ** (a + b) for q in t) for b in range(4)] for a in range(4)], [sum(e * q ** a for e, q in zip(E, t)) for a in range(4)])
    p = lambda q: c[0] + c[1] * q + c[2] * q * q + c[3] * q ** 3
    p1 = lambda q: c[1] + 2 * c[2] * q + 3 * c[3] * q * q
    p2 = lambda q: 2 * c[2] + 6 * c[3] * q
    tt = lambda v: (v / Vm) ** (-2. / 3) - 1
    t1 = lambda v: -2. / 3 / Vm * (v / Vm) ** (-5. / 3)
    t2 = lambda v: 10. / 9 / Vm ** 2 * (v / Vm) ** (-8. / 3)
    P = lambda v: -p1(tt(v)) * t1(v)
    Bm = lambda v: v * (p2(tt(v)) * t1(v) ** 2 + p1(tt(v)) * t2(v))
    D = 4 * c[2] ** 2 - 12 * c[1] * c[3]
    if D < 0: return None
    kok = [(-2 * c[2] + z * math.sqrt(D)) / (6 * c[3]) for z in (1, -1)] if abs(c[3]) > 1e-14 else [-c[1] / (2 * c[2])]
    kok = [q for q in kok if p2(q) > 0 and q > -1]
    if not kok: return None
    q0 = min(kok, key=abs); V0 = Vm * (1 + q0) ** -1.5; h = 1e-3 * V0
    Bp = (Bm(V0 + h) - Bm(V0 - h)) / (P(V0 + h) - P(V0 - h))
    res = max(abs(e - p(q)) for e, q in zip(E, t)) * 1000
    return dict(V0=V0, a0=V0 ** (1. / 3), B0=Bm(V0) * GPA, Bp=Bp, E0=p(q0), res=res, P=lambda v: P(v) * GPA)

def oku(sec):
    for y in sec:
        print('== %s (%s) ==' % (y, 'Bi25GaO39' if y == 'B' else 'Bi25GaO40'))
        print('%-12s %2s %2s %4s %9s %11s %8s %8s %8s %6s %9s' % ('ad', 'JD', 'BF', 'adim', 'V(A3)', 'E-E1(eV)', 'P(GPa)', 'Bt||111', 'GaO_ort', 'n<2.3', 'wall'))
        R, Z = [], {s: g.oku1('%s/04_eos/%s/run.out' % (B, ad(y, s))) for s in OL}
        E1 = Z[1.0]['E'] * RY if Z.get(1.0) else None
        for s in OL:
            n, r = ad(y, s), Z[s]
            if r is None: print('%-12s YOK' % n); continue
            a = A0 * s * B2A; V = a ** 3; E = r['E'] * RY
            X = g2.konum(r['t']); q = g2.geo(X, a) if len(X) == YAP[y][1] else None
            bg = [d for d in q['gaO'] if d < 2.3] if q else []
            ok = r['jd'] and r['bf']
            R.append((s, V, E, r['P'] / 10, ok))
            print('%-12s %2d %2d %4d %9.2f %11.4f %8.2f %8s %8s %6s %9s' % (n, r['jd'], r['bf'], r['ns'], V, E - E1 if E1 else float('nan'), r['P'] / 10,
                  '%+.3f' % q['kay'] if q else '-', '%.3f' % (sum(bg) / len(bg)) if bg else '-', len(bg) if q else '-', r['w']))
        G = [x for x in R if x[4]]
        if len(G) < len(R): print('### %d nokta JD/bfgs yok -> fitten cikarildi (K1)' % (len(R) - len(G)))
        if len(G) < 5: print('### fit icin en az 5 gecerli nokta gerekir'); continue
        f = bm3([x[1] for x in G], [x[2] for x in G])
        if f is None: print('### BM3 fitinde minimum bulunamadi'); continue
        print('BM3 E(V): V0 %.2f A3 (a0 %.4f A) | B0 %.1f GPa | B0p %.2f | maks artik %.2f meV' % (f['V0'], f['a0'], f['B0'], f['Bp'], f['res']))
        dP = [f['P'](x[1]) - x[3] for x in G]
        print('ic tutarlilik P_fit - P_QE (GPa): %s | maks |dP| %.2f' % (' '.join('%+.2f' % d for d in dP), max(abs(d) for d in dP)))
        print('V0 tarama araliginda: %s | maks P %.1f GPa' % (min(x[1] for x in G) < f['V0'] < max(x[1] for x in G), max(x[3] for x in G)))
        with open('%s/99_analiz/eos_ga%s.dat' % (B, y), 'w') as o:
            o.write('# s V_A3 E_eV P_QE_GPa JDBF | BM3 V0 %.4f B0 %.3f B0p %.4f E0 %.6f\n' % (f['V0'], f['B0'], f['Bp'], f['E0']))
            for x in R: o.write('%.3f %.4f %.6f %.4f %d\n' % (x[0], x[1], x[2], x[3], x[4]))
        print('yazildi: 99_analiz/eos_ga%s.dat' % y)

if __name__ == '__main__':
    if len(sys.argv) < 3 or sys.argv[1] not in ('uret', 'oku') or not set(sys.argv[2:]) <= set(YAP): g.dur('kullanim: ga_f4.py uret|oku B [S]')
    (uret if sys.argv[1] == 'uret' else oku)(sys.argv[2:])
