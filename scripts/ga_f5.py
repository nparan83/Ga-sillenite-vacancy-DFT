#!/usr/bin/env python3
# ga_f5.py : SilenitB F5 - eksik O siteleri (P23'te 6 esdegersiz O; A, B, C test edildi): D1 = O(2) 4e Ga tarafi (x~0.806),
# D2 = O(2) 4e Bi_tet tarafi (x~0.3145), Cp = ikinci 12j yorungesi; Br = B kontrolu (ayni baslangic ve sarsma).
# Baslangic: 02_O40/gaS_rel son konumlari - 1 O, tum atomlar +-0.03 A rastgele sarsma (tohum 20261124; simetri tuzagini kirmak icin).
# Sablon 03_site/gaB_rel/run.in (relax, 80/320, k2, nspin=1). Kullanim: python3 99_analiz/ga_f5.py uret | oku
import sys, os, re, math, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
import ga_f2 as g2
import ga_f2c as gc
B, RY = g.B, 13.605693
A = 19.58225862 * 0.529177
KAY = ('D1', 'D2', 'Cp', 'Br')
OPS = [lambda x, y, z: (x, y, z), lambda x, y, z: (-x, -y, z), lambda x, y, z: (-x, y, -z), lambda x, y, z: (x, -y, -z),
       lambda x, y, z: (z, x, y), lambda x, y, z: (z, -x, -y), lambda x, y, z: (-z, -x, y), lambda x, y, z: (-z, x, -y),
       lambda x, y, z: (y, z, x), lambda x, y, z: (-y, z, -x), lambda x, y, z: (y, -z, -x), lambda x, y, z: (-y, -z, x)]   # P23 nokta islemleri
C_SITE = (0.0117, 0.2484, 0.3657)                         # C (Makale A / F2), S cercevesinde (Ş1 analizi, 23.11)
def ad(k): return 'f5_%s' % k
def d(x, y): return g2.boy(g2.fark(x, y), A)

def siniflar(X):
    iO = [i for i, (e, _) in enumerate(X) if e == 'O']; ga = [x for e, x in X if e == 'Ga'][0]
    bt = min((d(x, [.5] * 3), i) for i, (e, x) in enumerate(X) if e == 'Bi')[1]
    S = {'A': [i for i in iO if d(X[i][1], ga) < 2.3], 'B': [i for i in iO if d(X[i][1], X[bt][1]) < 2.3]}
    e4 = [i for i in iO if i not in S['A'] + S['B'] and max(abs(c) for c in g2.fark(X[i][1], [0.] * 3)) - min(abs(c) for c in g2.fark(X[i][1], [0.] * 3)) < 0.01]
    S['D1'] = [i for i in e4 if d(X[i][1], ga) < d(X[i][1], X[bt][1])]; S['D2'] = [i for i in e4 if i not in S['D1']]
    c0 = min(iO, key=lambda i: d(X[i][1], C_SITE))
    yor = lambda j: {min(iO, key=lambda i: d(X[i][1], [q % 1 for q in op(*X[j][1])])) for op in OPS}
    S['C'] = sorted(yor(c0)); kalan = [i for i in iO if not any(i in v for v in S.values())]
    S['Cp'] = sorted(yor(kalan[0])) if kalan else []
    n = {k: len(v) for k, v in S.items()}
    if n != dict(A=4, B=4, D1=4, D2=4, C=12, Cp=12) or sorted(sum(S.values(), [])) != iO: g.dur('O yorunge sayimi %s' % n)
    sec = {'D1': max(S['D1'], key=lambda i: sum(X[i][1])), 'D2': min(S['D2'], key=lambda i: sum(X[i][1])), 'Cp': S['Cp'][0],
           'Br': max(S['B'], key=lambda i: sum(X[i][1]))}           # D1 (x,x,x) x~0.806; D2 (x,x,x) x~0.3145; Br = B'nin +[111] O'su
    return S, sec, bt

def uret():
    J = open(B + '/job_site.sh').read()
    if 'temizlik' not in J or 'wfc*.dat' in J: g.dur('job_site.sh temizliksiz ya da eski desen (wfc*.dat; HDF5 dalga fonksiyonlari kalir)')
    r = g.oku1(B + '/02_O40/gaS_rel/run.out')
    if r is None or not (r['jd'] and r['bf']): g.dur('gaS_rel JD/bfgs yok')
    X = g2.konum(r['t'])
    if len(X) != 66: g.dur('gaS_rel %d atom' % len(X))
    S, sec, bt = siniflar(X)
    T = open(B + '/03_site/gaB_rel/run.in').read()
    if not re.search(r"calculation\s*=\s*'relax'", T) or not re.search(r'nat\s*=\s*65\b', T): g.dur('sablon relax/nat=65 degil')
    L = T.splitlines(); i0 = [k for k, l in enumerate(L) if l.strip().upper().startswith('ATOMIC_POSITIONS')]
    if len(i0) != 1: g.dur('ATOMIC_POSITIONS tek degil')
    i0 = i0[0]; j0 = i0 + 1
    while j0 < len(L) and L[j0].split() and L[j0].split()[0] in ('Bi', 'Ga', 'O'): j0 += 1
    if j0 - i0 - 1 != 65: g.dur('sablonda %d konum' % (j0 - i0 - 1))
    for k in KAY:
        R = random.Random(20261124); v = sec[k]
        Y = [(e, [(q + R.uniform(-0.03, 0.03) / A) % 1 for q in x]) for m, (e, x) in enumerate(X) if m != v]
        t = '\n'.join(L[:i0 + 1] + ['%-3s %14.10f %14.10f %14.10f' % (e, *x) for e, x in Y] + L[j0:]) + '\n'
        t = g.setk(t, 'prefix', "'%s'" % ad(k))
        g.yaz('%s/03_site/%s' % (B, ad(k)), t)
        print('%s: bosaltilan O #%d %s | Ga %.2f A, Bi_tet %.2f A | en yakin Bi %s A | yazildi' % (ad(k), v, ' '.join('%.4f' % q for q in X[v][1]),
              d(X[v][1], [x for e, x in X if e == 'Ga'][0]), d(X[v][1], X[bt][1]), ' '.join('%.2f' % q for q in sorted(d(X[v][1], x) for e, x in X if e == 'Bi')[:3])))
    print('yorungeler (O indeksleri):', {k: v for k, v in S.items()})

def oku():
    rS = g.oku1(B + '/02_O40/gaS_rel/run.out'); X0 = g2.konum(rS['t']); S, sec, bt = siniflar(X0)
    ref = {k: g.oku1('%s/03_site/ga%s_rel/run.out' % (B, k)) for k in ('A', 'B', 'C')}
    EB = ref['B']['E']; print('baglam (F2): A %+.4f  C %+.4f eV (gaB_rel = 0)' % ((ref['A']['E'] - EB) * RY, (ref['C']['E'] - EB) * RY))
    E = {}
    for k in KAY:
        r = g.oku1('%s/03_site/%s/run.out' % (B, ad(k)))
        if r is None: print(ad(k), 'YOK'); continue
        X = g2.konum(r['t']); v = sec[k]; X0v = [x for m, x in enumerate(X0) if m != v]
        dmax = max((d(x, y[1]), e) for (e, x), y in zip(X, X0v)) if len(X) == 65 else (0, '-')
        dol = min(d(x, X0[v][1]) for e, x in X if e == 'O') if X else 0
        cev = gc.gacevre(X) if len(X) == 65 else []
        bts = d(X[bt - (1 if bt > v else 0)][1], [.5] * 3) if len(X) == 65 else 0
        print('%s: JD %d bfgs %d adim %d F %.5f P %+.2f kbar | dE(B) %+.4f eV | V_O yerine en yakin O %.2f A | maks kayma %.2f A (%s) | Ga %s | Bi_tet merkezden %.2f A | %s' % (
            ad(k), r['jd'], r['bf'], r['ns'], r['F'] or 0, r['P'] or 0, (r['E'] - EB) * RY, dol, dmax[0], dmax[1],
            ' '.join('%+.2f/%d' % (c[1], len(c[2])) for c in cev), bts, r['w']))
        if r['jd'] and r['bf']: E[k] = r['E']
    if len(E) < 4: print('### K1: 4/4 JD + bfgs gerekir'); return
    dBr = (E['Br'] - EB) * RY * 1000; Eref = max(EB, E['Br'])            # muhafazakar: marji kucuk yapan B referansi
    print('K2 kontrol: E(Br) - E(gaB_rel) = %+.1f meV -> %s' % (dBr, 'baslangic/sarsma bagimsiz' if abs(dBr) <= 10 else 'BAGIMLI: karar icin yuksek B referansi (marji kucuk olan) kullanilir'))
    for k in ('D1', 'D2', 'Cp'): print('   %s: E - E(gaB_rel) %+.3f eV | E - E(Br, ayni protokol) %+.3f eV' % (k, (E[k] - EB) * RY, (E[k] - E['Br']) * RY))
    m = min((E[k] - Eref) * RY for k in ('D1', 'D2', 'Cp'))
    print('K3: min(D1, D2, Cp) - B = %+.3f eV -> %s' % (m, 'B tam site setinde taban durum' if m >= 0.3 else ('B en dusuk, marj < 0.3 eV (metin zayiflar)' if m >= 0 else 'DUR: B taban durum DEGIL, ana iddia yeniden degerlendirilir')))

if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] not in ('uret', 'oku'): g.dur('kullanim: ga_f5.py uret | oku')
    uret() if sys.argv[1] == 'uret' else oku()
