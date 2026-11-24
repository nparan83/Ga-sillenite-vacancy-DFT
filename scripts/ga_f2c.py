#!/usr/bin/env python3
# ga_f2c.py : SilenitB F2c - Bi24Ga2O39 (ICDD 00-042-0180 bileşimi; iki M sitesi de Ga, Bi_tet yok, 1 O boslugu/hucre).
# Baslangic: Ga silenit A/B/C relaks son konumlari, Bi_tet -> Ga.  2A: bosluk Ga1 tetrahedronunda (gaA_rel'den),
# 2B: bosluk Ga2 tetrahedronunda (gaB_rel'den; 2A ile ayni site turu -> yol testi), 2C: kafes O boslugu (gaC_rel'den).
# Kullanim (dizin /arf/scratch/nparan/SilenitB): python3 99_analiz/ga_f2c.py uret | oku     Is: sbatch -J g2A job_site.sh ga2A_rel (B, C ayni)
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
import ga_f2 as g2
import ga_elek as ge
B, RY, B2A = g.B, 13.605693, 0.529177
A = 19.58225862 * B2A; SC = 10.183 / 10.3625
KAY = {'A': 'gaA_rel', 'B': 'gaB_rel', 'C': 'gaC_rel'}
E_B39, E_GA2O3 = None, -358.67210951           # Ry; Ga2O3: 06_ref/bGa2O3_scf (80/320) / 4 f.u.
B0 = 45.6                                       # GPa, Bi25GaO39 sikistirma kolu (durum 3f-ek3)
KART = {'Bi24Ga2O39 (00-042-0180)': 10.1760, 'Bi25GaO40 (00-046-0414)': 10.1689, 'Bi25GaO39 Radaev 1989': 10.183}

def uret():
    for k, src in KAY.items():
        r = g.oku1('%s/03_site/%s/run.out' % (B, src))
        if r is None or not (r['jd'] and r['bf']): g.dur(src + ': JD/bfgs yok')
        X = g2.konum(r['t'])
        if len(X) != 65: g.dur('%s: %d atom' % (src, len(X)))
        bt = g2.geo(X, A)['bt']; X = [('Ga' if i == bt else e, x) for i, (e, x) in enumerate(X)]
        n = {e: sum(1 for x, _ in X if x == e) for e in ('Bi', 'Ga', 'O')}
        if n != {'Bi': 24, 'Ga': 2, 'O': 39}: g.dur('%s: sayim %s' % (k, n))
        T = open('%s/03_site/%s/run.in' % (B, src)).read()
        t = g.setk(ge.konum_koy(T, X), 'prefix', "'ga2%s_rel'" % k)
        ne = int(round(sum(g.zval(__import__('re').search(r'(?im)^[ \t]*%s[ \t]+[0-9.]+[ \t]+(\S+\.upf)[ \t]*$' % e, t).group(1)) * n[e] for e in n)))
        if ne % 2: g.dur('tek elektron')
        g.yaz('%s/03_site/ga2%s_rel' % (B, k), t)
        print('ga2%s_rel: kaynak %s, Bi_tet #%d -> Ga | Bi24 Ga2 O39, N = %d | yazildi' % (k, src, bt + 1, ne))

def gacevre(X):
    out = []
    for i, (e, x) in enumerate(X):
        if e != 'Ga': continue
        c = [0.] * 3 if g2.boy(g2.fark(x, [0.] * 3), A) < g2.boy(g2.fark(x, [.5] * 3), A) else [.5] * 3
        d = g2.fark(x, c); par = A * sum(d) / math.sqrt(3)
        O = sorted(g2.boy(g2.fark(y, x), A) for s, y in X if s == 'O')
        out.append(('orijin' if c[0] == 0 else 'merkez', par, [o for o in O if o < 2.3], O[:5]))
    return out

def oku():
    rb = g.oku1('%s/03_site/gaB_rel/run.out' % B); R = {}
    for k in KAY:
        r = g.oku1('%s/03_site/ga2%s_rel/run.out' % (B, k))
        if r is None: print('ga2%s_rel YOK' % k); continue
        print('ga2%s_rel: JD %d bfgs %d adim %d E %.8f Ry P %.2f kbar F %.5f | %s' % (k, r['jd'], r['bf'], r['ns'], r['E'], r['P'], r['F'], r['w']))
        X = g2.konum(r['t'])
        if len(X) == 65:
            for s, par, b, O in gacevre(X): print('   Ga (%s): ||[111] %+.3f | Ga-O<2.3: %d bag, ort %.3f | en yakin O %s' % (s, par, len(b), sum(b) / max(len(b), 1), ' '.join('%.3f' % o for o in O)))
        if r['jd'] and r['bf']: R[k] = r
    if not R: return
    k0 = min(R, key=lambda k: R[k]['E'])
    print('en dusuk: ga2%s | ' % k0 + ' '.join('dE(%s) %+.4f eV' % (k, (R[k]['E'] - R[k0]['E']) * RY) for k in sorted(R)))
    if 'A' in R and 'B' in R: print('yol testi (2A ve 2B ayni site turu): |dE| = %.1f meV (> 5 meV ise yol bagimliligi)' % (abs(R['A']['E'] - R['B']['E']) * RY * 1000))
    dEd = (R[k0]['E'] - 24 / 25 * rb['E'] - 13 / 25 * E_GA2O3) * RY
    print('ayrisma: Bi24Ga2O39 -> 24/25 Bi25GaO39 + 13/25 Ga2O3 : dE = %+.3f eV/f.u. (%+.1f meV/atom; > 0 ise Bi24Ga2O39 kararsiz, PBE oksit->oksit, duzeltmesiz)' % (dEd, dEd / 65 * 1000))
    dP = (R[k0]['P'] - rb['P']) / 10; da = A * dP / (3 * B0)
    print('kafes: P(2) - P(Bi25GaO39) = %+.2f GPa (ayni hucre, 80 Ry) -> da ~ %+.4f A (PBE), olcekli %+.4f A (B0 %.1f GPa, birinci mertebe)' % (dP, da, da * SC, B0))
    print('deney a: ' + ' | '.join('%s %.4f' % kv for kv in KART.items()) + ' -> kart farki 42-0180 - 46-0414 = %+.4f A' % (KART['Bi24Ga2O39 (00-042-0180)'] - KART['Bi25GaO40 (00-046-0414)']))

if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] not in ('uret', 'oku'): g.dur('kullanim: ga_f2c.py uret | oku')
    uret() if sys.argv[1] == 'uret' else oku()
