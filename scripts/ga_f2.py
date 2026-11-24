#!/usr/bin/env python3
# ga_f2.py : SilenitB F2 - Ga silenit A/B/C bosluk relakslari (uret) ve analiz (oku)
# Kullanim (dizin /arf/scratch/nparan/SilenitB):  python3 99_analiz/ga_f2.py uret | oku
# ga_f1.py ile ayni sablon/ayarlar (safrel.in, 80/320, k 2x2x2, nspin=1); konumlar 08_paw_U/pos{A,B,C}_relax.pos (Fe -> Ga)
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
g.POS['C'] = 'C'
B = g.B
SITE = ('A', 'B', 'C')
RY = 13.605693
A_EXP = 10.183            # Radaev 1989, Bi25GaO39 notron (Dokl. 306, 624)
EXP = {'kayma': 0.3175, 'GaO': 1.879}   # Radaev 1989: Bi_T kaymasi (0.018*sqrt3*a), Ga-O
E_S40 = -4806.81440564    # 02_O40/gaS_rel (80 Ry, k2) - yalniz bilgi

JOB = g.HEAD.replace('@J@', 'gaB_site').replace('@T@', '12:00:00') + """cd @B@/03_site/$1 || exit 1
[ -s run.out ] && { echo "### run.out zaten var -> DUR"; exit 1; }
mpirun pw.x -npool 2 -in run.in > run.out
"""

def uret():
    T = open(g.TPL).read()
    print('%-8s %-5s %4s %4s %5s' % ('ad', 'calc', 'ecut', 'nat', 'nelec'))
    for s in SITE:
        n = 'ga%s_rel' % s
        t, ne, kay, sp = g.uret1(T, n, s, 'relax', 80, 2)
        g.yaz('%s/03_site/%s' % (B, n), t)
        print('%-8s %-5s %4d %4d %5d' % (n, 'relax', 80, 65, ne))
    open(B + '/job_site.sh', 'w').write(JOB.replace('@B@', B))
    print('yazildi: 03_site/ga{A,B,C}_rel/run.in, job_site.sh (her site AYRI is: sbatch -J gaB_X job_site.sh gaX_rel)')

def konum(t):
    k = t.rfind('ATOMIC_POSITIONS'); X = []
    if k < 0: return X
    for l in t[k:].splitlines()[1:]:
        w = l.split()
        if len(w) < 4 or w[0] not in ('Bi', 'Ga', 'O'): break
        X.append((w[0], [float(q) for q in w[1:4]]))
    return X

def fark(x, y): return [((p - q) + .5) % 1 - .5 for p, q in zip(x, y)]
def boy(d, a): return a * math.sqrt(sum(c * c for c in d))

def komsu(X, i, a, n=5):
    return sorted(boy(fark(x, X[i][1]), a) for s, x in X if s == 'O')[:n]

def geo(X, a):
    gi = [i for i, (s, x) in enumerate(X) if s == 'Ga']
    if len(gi) != 1: g.dur('Ga sayisi %d' % len(gi))
    bi = [i for i, (s, x) in enumerate(X) if s == 'Bi']
    bt = min(bi, key=lambda i: sum(c * c for c in fark(X[i][1], [.5] * 3)))
    u = fark(X[bt][1], [.5] * 3); v = fark(X[gi[0]][1], [0.] * 3)
    par = lambda d: a * sum(d) / math.sqrt(3)
    return dict(bt=bt, kay=par(u), kaymod=boy(u, a), gak=par(v), gamod=boy(v, a),
                btO=komsu(X, bt, a), gaO=komsu(X, gi[0], a))

def oku():
    s_ = 10.183 / 10.3625
    R = {}
    print('%-8s %2s %2s %4s %17s %7s %9s %10s' % ('ad', 'JD', 'BF', 'adim', 'E(Ry)', 'P', 'F(Ry/b)', 'wall'))
    for s in SITE:
        r = g.oku1('%s/03_site/ga%s_rel/run.out' % (B, s))
        if r is None: print('ga%s_rel  YOK' % s); continue
        R[s] = r
        print('%-8s %2d %2d %4d %17.8f %7.2f %9.5f %10s' % ('ga%s_rel' % s, r['jd'], r['bf'], r['ns'], r['E'], r['P'], r['F'], r['w']))
    ok = {s for s in R if R[s]['jd'] and R[s]['bf']}
    if 'B' in ok:
        for s in ('A', 'C'):
            if s in ok: print('dE(%s-B) = %+.4f eV   (Fe, Makale A: A +0.636, C +0.822)' % (s, (R[s]['E'] - R['B']['E']) * RY))
            else: print('dE(%s-B) : %s JD/bfgs yok -> HESAPLANMADI (K1)' % (s, s))
    else: print('### B JD/bfgs yok -> dE hesaplanmadi (K1)')
    a = 19.58225862 * 0.529177
    print('geometri (a = %.4f A; olcekli = x%.4f -> a_exp %.3f):' % (a, s_, A_EXP))
    for s in SITE:
        if s not in R: continue
        X = konum(R[s]['t'])
        if len(X) != 65: print('ga%s_rel: %d atom okundu -> geometri atlandi' % (s, len(X))); continue
        q = geo(X, a); f = '' if s in ok else '  [YAKINSAMAMIS - kullanma]'
        print('%s: Bi_tet(#%d) kayma||[111] %+.3f (|u| %.3f; olcekli %+.3f) | Bi_tet-O %s%s' % (
            s, q['bt'] + 1, q['kay'], q['kaymod'], q['kay'] * s_, ' '.join('%.3f' % d for d in q['btO']), f))
        bg = [d for d in q['gaO'] if d < 2.3]
        print('   Ga kayma||[111] %+.3f (|u| %.3f) | Ga-O %s | <2.3 A: %d bag, olcekli ort %.3f' % (
            q['gak'], q['gamod'], ' '.join('%.3f' % d for d in q['gaO']), len(bg), s_ * sum(bg) / max(len(bg), 1)))
    print('Radaev 1989 (Ga, notron): Bi_T kayma %.4f A, Ga-O %.3f A; Fe/PBE B: kayma +0.321 (olcekli +0.315)' % (EXP['kayma'], EXP['GaO']))

if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] not in ('uret', 'oku'): g.dur('kullanim: ga_f2.py uret | oku')
    uret() if sys.argv[1] == 'uret' else oku()
