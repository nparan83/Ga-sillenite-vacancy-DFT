#!/usr/bin/env python3
# ga_f2a.py : SilenitB F2a - A bosluk, TERS DONMEMIS baslangic (gaS_rel.pos - Fe posA boslugundaki O; v2 13.11)
# Kullanim (dizin /arf/scratch/nparan/SilenitB):  python3 99_analiz/ga_f2a.py uret | oku
# Ayarlar ga_f2 ile ayni (safrel.in, 80/320, k2, nspin=1). Is: sbatch -J gaB_A0 job_site.sh gaA0_rel
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
import ga_f2 as g2
B, RY, N = g.B, 13.605693, 'gaA0_rel'
A_PBE = 19.58225862 * 0.529177

def oku_pos(f):
    X = []
    for l in open(f):
        w = l.split()
        if len(w) >= 4 and w[0] in ('Bi', 'Ga', 'Fe', 'O'): X.append((w[0], [float(q) for q in w[1:4]]))
    return X

def uret():
    S = oku_pos(B + '/02_O40/gaS_rel.pos')
    if len(S) != 66: g.dur('gaS_rel.pos %d atom' % len(S))
    gi = [i for i, (s, x) in enumerate(S) if s == 'Ga']
    if len(gi) != 1: g.dur('Ga sayisi %d' % len(gi))
    xg = S[gi[0]][1]
    if g2.boy(g2.fark(xg, [0.] * 3), A_PBE) > 0.05: g.dur('Ga orijinde degil')
    c = sorted((g2.boy(g2.fark(x, xg), A_PBE), i) for i, (s, x) in enumerate(S) if s == 'O')[:4]
    Fa = oku_pos(g.F + '/08_paw_U/posA_relax.pos')
    FO = [x for s, x in Fa if s == 'O']
    if len(FO) != 39: g.dur('Fe posA O sayisi %d' % len(FO))
    print('Ga tetrahedronu O adaylari (gaS_rel):')
    R = []
    for d, i in c:
        u = g2.fark(S[i][1], xg); par = A_PBE * sum(u) / math.sqrt(3)
        mf = min(g2.boy(g2.fark(S[i][1], y), A_PBE) for y in FO)
        R.append((mf, par, d, i)); print('  O#%-3d Ga-O %.3f  ||[111] %+.3f  Fe-posA en yakin O %.3f A' % (i + 1, d, par, mf))
    R.sort(reverse=True); k1 = R[0][3]
    if R[0][0] < 0.5 or R[1][0] > 0.5: g.dur('Fe-A boslugu net degil (tek aday yok)')
    if abs(abs(R[0][1]) - R[0][2]) > 0.05: g.dur('cikarilacak O Ga uzerindeki [111] ekseninde degil')
    if R[0][1] < 0: g.dur('bosluk -[111] yonunde: oku() icindeki -0.753 isaret varsayimi gecersiz')
    P = [(s, x) for i, (s, x) in enumerate(S) if i != k1]
    d = B + '/99_analiz/_f2a/08_paw_U'; os.makedirs(d, exist_ok=True)
    # uret1 yalniz Bi/Fe/O okur ve Fe->Ga cevirir: gecici dosyada Ga 'Fe' etiketiyle yazilir
    open(d + '/posA0_relax.pos', 'w').write(''.join('%-3s %14.10f %14.10f %14.10f\n' % ('Fe' if s == 'Ga' else s, *x) for s, x in P))
    g.POS['A0'] = 'A0'; F0 = g.F; g.F = B + '/99_analiz/_f2a'
    t, ne, kay, sp = g.uret1(open(g.TPL).read(), N, 'A0', 'relax', 80, 2)
    g.F = F0
    g.yaz(B + '/03_site/' + N, t)
    print('cikarilan O#%d | %s relax 80/320 k2 nat 65 nelec %d | yazildi: 03_site/%s/run.in' % (k1 + 1, N, ne, N))

def oku():
    r = g.oku1(B + '/03_site/%s/run.out' % N)
    if r is None: g.dur(N + ' run.out yok')
    print('%s: JD %d bfgs %d adim %d E %.8f Ry P %.2f F %.5f wall %s' % (N, r['jd'], r['bf'], r['ns'], r['E'], r['P'], r['F'], r['w']))
    ra, rb = g.oku1(B + '/03_site/gaA_rel/run.out'), g.oku1(B + '/03_site/gaB_rel/run.out')
    X = g2.konum(r['t'])
    if len(X) != 65: g.dur('%d atom okundu' % len(X))
    q = g2.geo(X, A_PBE)
    print('Ga ||[111] %+.3f (|u| %.3f) | Ga-O %s' % (q['gak'], q['gamod'], ' '.join('%.3f' % x for x in q['gaO'])))
    print('Bi_tet ||[111] %+.3f | Bi_tet-O %s' % (q['kay'], ' '.join('%.3f' % x for x in q['btO'])))
    if not (r['jd'] and r['bf']): print('### JD/bfgs yok -> KARAR VERILMEZ (K1)'); return
    dA = (r['E'] - ra['E']) * RY * 1000; dB = (r['E'] - rb['E']) * RY
    print('E(A0) - E(gaA_rel) = %+.1f meV | E(A0) - E(gaB_rel) = %+.4f eV  (gaA_rel: Ga -0.753, A-B +0.7436)' % (dA, dB))
    if abs(dA) <= 5 and abs(q['gak'] + 0.753) <= 0.05: print('KARAR R1: ters donme KENDILIGINDEN (ayni minimum)')
    elif dA < -5: print('KARAR R3: ALARM - gaA_rel global minimum DEGIL; A-B degisir, Fe-A da sinanmali')
    elif dA > 5 and abs(q['gak']) < 0.3: print('KARAR R2: iki minimum; ters donmus olan %.0f meV daha dusuk' % dA)
    else: print('KARAR: R1/R2/R3 disi -> geometri elle incelenecek')

if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] not in ('uret', 'oku'): g.dur('kullanim: ga_f2a.py uret | oku')
    uret() if sys.argv[1] == 'uret' else oku()
