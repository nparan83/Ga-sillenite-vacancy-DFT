#!/usr/bin/env python3
# ga_f5b.py : SilenitB F5-b - f5_Cp (C' bosluk) 16 saatte yakinsamadi (144 adim, TIMEOUT) ve oku ciktisi
# V_O yerine bir O'nun girdigini (0.06 A) ve dE(B) = +1.3 meV oldugunu gosterdi -> bosluk baska siteye goc etmis olabilir.
#   tani : f5_Cp ve Br'nin son geometrisinde bos siteyi bulur, yorungesini (A/B/C/Cp/D1/D2) ve gaB_rel'e simetri
#          eslesmeli RMS farkini verir; adim adim bos site + enerji (goc hangi adimda, enerji artisi var mi). 0 ç-saat.
#   uret : tani + 03_site/f5_Cp2 (f5_Cp son konumlarindan devam, ayni ayarlar). oku : K1' ve K4 karari.
# Onceden yazilan kriterler (25.11, sonuctan once): K1' JD + bfgs converged. K4: bos site B yorungesinde VE
# |E - E(gaB_rel)| <= 5 meV -> "C' yerel minimum degil, bosluk B'ye gecer"; B yorungesi ama E < E_B - 5 meV veya baska
# yorunge ve E < E_B -> DUR (B referansi degisir); diger durumlar ayri konfigurasyon olarak raporlanir.
# RMS(gaB_rel) karar olcutu DEGIL; Br'nin RMS'i kalibrasyon olarak yaninda basilir. Kullanim: python3 99_analiz/ga_f5b.py tani|uret|oku
import sys, os, re, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
import ga_f2 as g2
import ga_f5 as f5
B, RY, A, d = g.B, 13.605693, f5.A, f5.d
S3 = B + '/03_site'

def bloklar(t):
    out = []
    for m in re.finditer(r'ATOMIC_POSITIONS', t):
        X = []
        for l in t[m.start():].splitlines()[1:]:
            w = l.split()
            if len(w) < 4 or w[0] not in ('Bi', 'Ga', 'O'): break
            X.append((w[0], [float(q) for q in w[1:4]]))
        out.append(X)
    return out

def bos(X, X0, S):
    iO = [i for i, (e, _) in enumerate(X0) if e == 'O']
    b = [i for i in iO if min(d(x, X0[i][1]) for e, x in X if e == 'O') > 1.0]
    return [(i, [k for k, v in S.items() if i in v][0]) for i in b]

def eslestir(Xa, Xb, X0, S):
    ba, bb = bos(Xa, X0, S), bos(Xb, X0, S)
    if len(ba) != 1 or len(bb) != 1: return None
    pa, pb, en = X0[ba[0][0]][1], X0[bb[0][0]][1], None
    for op in f5.OPS:
        if d([q % 1 for q in op(*pa)], pb) > 0.5: continue
        Y = [(e, [q % 1 for q in op(*x)]) for e, x in Xa]
        r = [min(d(y, z) for f, z in Xb if f == e) for e, y in Y]
        s = (math.sqrt(sum(q * q for q in r) / len(r)), max(r))
        en = s if en is None or s[0] < en[0] else en
    return en

def baglam():
    X0 = g2.konum(g.oku1(B + '/02_O40/gaS_rel/run.out')['t']); S, sec, bt = f5.siniflar(X0)
    rB = g.oku1(S3 + '/gaB_rel/run.out'); XB = g2.konum(rB['t'])
    if bos(XB, X0, S) == [] or bos(XB, X0, S)[0][1] != 'B': g.dur('gaB_rel bos sitesi B degil: %s' % bos(XB, X0, S))
    return X0, S, sec, rB['E'], XB

def tani(k, X0, S, sec, EB, XB):
    r = g.oku1('%s/%s/run.out' % (S3, k))
    if r is None: print(k, 'YOK'); return None
    Bl = bloklar(r['t']); X = Bl[-1]; b = bos(X, X0, S); es = eslestir(X, XB, X0, S)
    print('%s: JD %d bfgs %d adim %d F %.5f | son ! enerji dE(B) %+.1f meV | bos site %s | RMS/maks fark (gaB_rel, simetri) %s' % (
        k, r['jd'], r['bf'], r['ns'], r['F'] or 0, (r['E'] - EB) * RY * 1000,
        ' '.join('#%d %s (%s)' % (i, c, ' '.join('%.3f' % q for q in X0[i][1])) for i, c in b) or 'YOK',
        '%.3f/%.3f A' % es if es else '-'))
    return r, Bl, b

def uret(yaz=True):
    X0, S, sec, EB, XB = baglam(); v = sec['Cp']
    print("C' bosaltilan: #%d %s | B yorungesi: %s" % (v, ' '.join('%.3f' % q for q in X0[v][1]), S['B']))
    tani('f5_Br', X0, S, sec, EB, XB); r, Bl, b = tani('f5_Cp', X0, S, sec, EB, XB)
    En = [float(q) for q in re.findall(r'^!\s+total energy\s+=\s+(-?[\d.]+)', r['t'], re.M)]
    gec = [(n, bos(X, X0, S)) for n, X in enumerate(Bl)]; ilk = next((n for n, q in gec if [i for i, c in q] != [v]), None)
    yuk = [(n + 1, (En[n + 1] - En[n]) * RY * 1000) for n in range(len(En) - 1) if En[n + 1] - En[n] > 1e-3 / RY]
    print('enerji (meV, gaB_rel = 0) adim 1/10/25/50/100/son: %s' % ' '.join('%+.0f' % ((En[n] - EB) * RY * 1000) for n in (0, 9, 24, 49, 99, -1) if n < len(En)))
    print('bos site ilk degistigi blok: %s -> %s | >1 meV enerji artisi olan adim: %d %s' % (ilk, gec[ilk][1] if ilk is not None else '-', len(yuk), yuk[:6]))
    if not yaz: return
    T = open(S3 + '/f5_Cp/run.in').read(); L = T.splitlines()
    i0 = [n for n, l in enumerate(L) if l.strip().upper().startswith('ATOMIC_POSITIONS')]
    if len(i0) != 1 or len(Bl[-1]) != 65: g.dur('f5_Cp girdi/son blok uygun degil (%d atom)' % len(Bl[-1]))
    i0 = i0[0]; j0 = i0 + 1
    while j0 < len(L) and L[j0].split() and L[j0].split()[0] in ('Bi', 'Ga', 'O'): j0 += 1
    t = '\n'.join(L[:i0 + 1] + ['%-3s %14.10f %14.10f %14.10f' % (e, *x) for e, x in Bl[-1]] + L[j0:]) + '\n'
    g.yaz(S3 + '/f5_Cp2', g.setk(t, 'prefix', "'f5_Cp2'")); print('03_site/f5_Cp2/run.in yazildi (f5_Cp son blogu, %d. adim sonrasi)' % len(Bl))

def oku():
    X0, S, sec, EB, XB = baglam(); tani('f5_Br', X0, S, sec, EB, XB); q = tani('f5_Cp2', X0, S, sec, EB, XB)
    if q is None: return
    r, Bl, b = q; dE = (r['E'] - EB) * RY * 1000
    if not (r['jd'] and r['bf']): print("### K1': f5_Cp2 JD + bfgs yok -> karar verilmez"); return
    if len(b) == 1 and b[0][1] == 'B' and abs(dE) <= 5: print("K4: C' yerel minimum DEGIL; bosluk B'ye gecer (dE %+.1f meV). K3 tamamlanan siteler: min(D1, D2) - B ile karar." % dE)
    elif (len(b) == 1 and b[0][1] == 'B' and dE < -5) or dE < 0: print('### K4: B referansindan dusuk konfigurasyon (%+.1f meV, bos %s) -> DUR' % (dE, b))
    else: print('K4: ayri konfigurasyon, bos %s, dE %+.1f meV -> raporlanir' % (b, dE))

if __name__ == '__main__':
    m = sys.argv[1] if len(sys.argv) == 2 else ''
    if m not in ('tani', 'uret', 'oku'): g.dur('kullanim: ga_f5b.py tani | uret | oku')
    uret(False) if m == 'tani' else uret() if m == 'uret' else oku()
