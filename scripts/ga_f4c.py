#!/usr/bin/env python3
# ga_f4c.py : SilenitB F4c - B yapisi s=1.005, iki baslangic: c = merkezli (gaB_s1.000_f1.010 son konumlari),
# o = merkez-disi (gaB_s1.010 son konumlari, Ga -0.104). + merkezli kolun BM3 fiti (0.94-0.99, 1.00 en dusuk, 1.005 merkezliyse).
# Kullanim (dizin /arf/scratch/nparan/SilenitB): python3 99_analiz/ga_f4c.py uret | oku
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
import ga_f2 as g2
import ga_f4 as g4
B, RY, A0, B2A = g.B, 13.605693, 19.58225862, 0.529177
S5 = 1.005
BAS = {'c': 'gaB_s1.000_f1.010', 'o': 'gaB_s1.010'}
def ad(k): return 'gaB_s%.3f_%s' % (S5, k)
def V(s): return (A0 * s * B2A) ** 3

def son(n):
    r = g.oku1('%s/04_eos/%s/run.out' % (B, n))
    if r is None or not (r['jd'] and r['bf']): g.dur(n + ': run.out yok ya da JD/bfgs yok (K1)')
    return r

def uret():
    L = open('%s/04_eos/gaB_s1.000/run.in' % B).read().splitlines()
    i = [q for q, l in enumerate(L) if l.strip().upper().startswith('ATOMIC_POSITIONS')]
    if len(i) != 1 or 'crystal' not in L[i[0]].lower(): g.dur('ATOMIC_POSITIONS (crystal) tek degil')
    i = i[0]; j = i + 1
    while j < len(L) and L[j].split() and L[j].split()[0] in ('Bi', 'Ga', 'O'): j += 1
    if j - i - 1 != 65: g.dur('sablonda %d konum' % (j - i - 1))
    for k, n in BAS.items():
        X = g2.konum(son(n)['t'])
        if len(X) != 65: g.dur('%s: %d atom' % (n, len(X)))
        t = '\n'.join(L[:i + 1] + ['%-3s %14.10f %14.10f %14.10f' % (e, *x) for e, x in X] + L[j:]) + '\n'
        t = g.setk(g.setk(t, 'celldm(1)', '%.8f' % (A0 * S5)), 'prefix', "'%s'" % ad(k))
        q = g2.geo(X, A0 * S5 * B2A)
        g.yaz('%s/04_eos/%s' % (B, ad(k)), t)
        print('%s: a %.5f A, baslangic Ga ||[111] %+.3f, Bi_tet %+.3f (kaynak %s) | yazildi' % (ad(k), A0 * S5 * B2A, q['gak'], q['kay'], n))
    if 'temizlik' not in open(B + '/job_eos.sh').read(): print('### UYARI: job_eos.sh temizliksiz -> SilenitB_kodlari.md 9. bolumdeki surumu yaz')

def oku():
    R = {}
    for k in BAS:
        r = g.oku1('%s/04_eos/%s/run.out' % (B, ad(k)))
        if r is None: print(ad(k), 'YOK'); continue
        q = g2.geo(g2.konum(r['t']), A0 * S5 * B2A); R[k] = (r, q)
        print('%s: JD %d bfgs %d adim %d E %.8f Ry P %.2f GPa | Ga %+.3f Bi_tet %+.3f | Ga-O %s | %s' % (ad(k), r['jd'], r['bf'], r['ns'], r['E'], r['P'] / 10,
              q['gak'], q['kay'], ' '.join('%.3f' % d for d in q['gaO'][:4]), r['w']))
    G = {k: v for k, v in R.items() if v[0]['jd'] and v[0]['bf']}
    if len(G) == 2:
        dE = (G['o'][0]['E'] - G['c'][0]['E']) * RY * 1000; dg = G['o'][1]['gak'] - G['c'][1]['gak']
        print('E(o) - E(c) = %+.2f meV | Ga(o) - Ga(c) = %+.3f A' % (dE, dg))
        if abs(dE) <= 1 and abs(dg) < 0.02: print('KARAR C1: tek minimum (iki baslangic ayni yere)')
        elif dE > 1: print('KARAR C2: merkez-disi YUKSEK (yarikararli) -> merkezli kol gecerli')
        elif dE < -1: print('KARAR C3: merkez-disi DUSUK -> gecis s = 1.000-1.005 arasinda')
        else: print('KARAR C?: elle incelenecek')
    P = [(s, V(s), son('gaB_s%.3f' % s)['E'] * RY, son('gaB_s%.3f' % s)['P'] / 10) for s in (0.94, 0.96, 0.98, 0.99)]
    r1 = min((son('gaB_s1.000'), son('gaB_s1.000_f1.010')), key=lambda r: r['E'])
    P.append((1.0, V(1.0), r1['E'] * RY, r1['P'] / 10))
    if G:
        k = min(G, key=lambda k: G[k][0]['E'])
        if abs(G[k][1]['gak']) < 0.05: P.append((S5, V(S5), G[k][0]['E'] * RY, G[k][0]['P'] / 10)); print('kola alindi: s=1.005 (%s)' % k)
        else: print('s=1.005 en dusuk durum merkez-disi -> merkezli kola ALINMADI')
    f = g4.bm3([x[1] for x in P], [x[2] for x in P])
    if f is None: print('### BM3 minimum yok'); return
    dP = max(abs(f['P'](x[1]) - x[3]) for x in P); ic = min(x[1] for x in P) < f['V0'] < max(x[1] for x in P); neg = min(x[3] for x in P)
    print('MERKEZLI KOL BM3 (%d nokta): V0 %.2f A3 (a0 %.4f) | B0 %.1f GPa | B0p %.2f | artik %.2f meV | maks |dP| %.2f GPa | V0 icte %s | min P %.2f GPa' % (
          len(P), f['V0'], f['a0'], f['B0'], f['Bp'], f['res'], dP, ic, neg))
    ok = [f['res'] <= 2, dP <= 0.3, ic and neg <= -0.3, len(P) >= 6]
    print("E2' %s  E3' %s  E4' %s  n>=6 %s  -> %s" % (*('OK' if x else 'X' for x in ok), 'B0 ve B0p RAPORLANIR' if all(ok) else 'yalniz B0 (B0p yok)'))

if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] not in ('uret', 'oku'): g.dur('kullanim: ga_f4c.py uret | oku')
    uret() if sys.argv[1] == 'uret' else oku()
