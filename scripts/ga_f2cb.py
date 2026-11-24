#!/usr/bin/env python3
# ga_f2cb.py : SilenitB F2c-b - Bi24Ga2O39 (ga2A) kafes parametresi sağlaması: sabit hacim relaks, 90/360, s = 0.97-1.00.
# a0 = P(V)=0 noktasi (4 noktaya ikinci derece P(V) fiti); ayni yontem Bi25GaO39 (gaB) EOS noktalarina da uygulanir.
# Kullanim (dizin /arf/scratch/nparan/SilenitB): python3 99_analiz/ga_f2cb.py uret | oku     Is: sbatch -J e2A0.97 job_eos.sh ga2A_s0.970 (4 nokta)
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
import ga_f2 as g2
import ga_f4 as g4
import ga_elek as ge
import ga_f2c as gc
B, RY, B2A, A0 = g.B, 13.605693, 0.529177, 19.58225862
SC = 10.183 / 10.3625
OL = (0.97, 0.98, 0.99, 1.00)
def ad(s): return 'ga2A_s%.3f' % s
def V(s): return (A0 * s * B2A) ** 3

def uret():
    if 'temizlik' not in open(B + '/job_eos.sh').read(): g.dur('job_eos.sh temizliksiz (SilenitB_kodlari 9)')
    r = g.oku1(B + '/03_site/ga2A_rel/run.out')
    if r is None or not (r['jd'] and r['bf']): g.dur('ga2A_rel JD/bfgs yok')
    X = g2.konum(r['t'])
    if len(X) != 65 or sum(1 for e, _ in X if e == 'Ga') != 2: g.dur('ga2A_rel atom/Ga sayimi')
    T = ge.konum_koy(open(B + '/03_site/ga2A_rel/run.in').read(), X)
    for s in OL:
        t = T
        for k, v in (('celldm(1)', '%.8f' % (A0 * s)), ('ecutwfc', '90.0'), ('ecutrho', '360.0'), ('prefix', "'%s'" % ad(s)), ('outdir', "'./out'"), ('tstress', '.true.')):
            t = g.setk(t, k, v)
        g.yaz('%s/04_eos/%s' % (B, ad(s)), t); print('%s: a %.5f A, 90/360 relax | yazildi' % (ad(s), A0 * s * B2A))

def p0(pts):                                      # P(V) = c0 + c1 v + c2 v^2, v = V - Vm ; kok
    Vm = sum(v for v, _ in pts) / len(pts); x = [v - Vm for v, _ in pts]
    c = g4.coz([[sum(q ** (a + b) for q in x) for b in range(3)] for a in range(3)], [sum(p * q ** a for (v, p), q in zip(pts, x)) for a in range(3)])
    D = c[1] ** 2 - 4 * c[2] * c[0]; r = [(-c[1] + z * D ** .5) / (2 * c[2]) for z in (1, -1)] if abs(c[2]) > 1e-12 else [-c[0] / c[1]]
    v0 = min(r, key=abs) + Vm; res = max(abs(c[0] + c[1] * q + c[2] * q * q - p) for (v, p), q in zip(pts, x))
    return v0, res

def oku():
    G = []
    for s in OL:
        r = g.oku1('%s/04_eos/%s/run.out' % (B, ad(s)))
        if r is None: print(ad(s), 'YOK'); continue
        X = g2.konum(r['t']); cev = gc.gacevre(X) if len(X) == 65 else []
        print('%s: JD %d bfgs %d adim %d P %.2f GPa | %s | %s' % (ad(s), r['jd'], r['bf'], r['ns'], r['P'] / 10,
              ' '.join('Ga(%s) %+.3f/%d bag' % (c[0], c[1], len(c[2])) for c in cev), r['w']))
        if r['jd'] and r['bf']: G.append((V(s), r['P'] / 10))
    Bp = []
    for s, ns in ((0.96, ['gaB_s0.960']), (0.98, ['gaB_s0.980']), (0.99, ['gaB_s0.990', 'gaB_s0.990_f1.010']), (1.00, ['gaB_s1.000', 'gaB_s1.000_f1.010'])):
        R = [x for x in (g.oku1('%s/04_eos/%s/run.out' % (B, n)) for n in ns) if x and x['jd'] and x['bf']]
        if R: Bp.append((V(s), min(R, key=lambda x: x['E'])['P'] / 10))
    if len(G) < 4 or len(Bp) < 4: print('### en az 4 + 4 gecerli nokta gerekir (K1)'); return
    vg, rg = p0(G); vb, rb = p0(Bp); ag, ab = vg ** (1 / 3), vb ** (1 / 3)
    print('Bi24Ga2O39 (2A): V0 %.2f A3, a0 %.4f A (P fit maks artik %.3f GPa)' % (vg, ag, rg))
    print('Bi25GaO39  (B) : V0 %.2f A3, a0 %.4f A (P fit maks artik %.3f GPa; BM3 sikistirma kolu V0 1111.40)' % (vb, ab, rb))
    d = (ag - ab) * SC
    print('Da = a0(Bi24Ga2O39) - a0(Bi25GaO39) = %+.4f A (PBE), %+.4f A olcekli | birinci mertebe tahmin -0.0307 A' % (ag - ab, d))
    print('KARAR: %s' % ('|Da| > 0.02 A -> kart 42-0180 a degeri tam Ga2 bilesimiyle UYUMSUZ' if abs(d) > 0.02 else '|Da| <= 0.02 A -> kafes parametresi bilesimi AYIRT ETMEZ'))
    print('beklenen a(Bi24Ga2O39): Radaev Bi25GaO39 10.183 -> %.3f | 46-0414 Bi25Ga 10.1689 -> %.3f | kart 42-0180: 10.1760' % (10.183 + d, 10.1689 + d))

if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] not in ('uret', 'oku'): g.dur('kullanim: ga_f2cb.py uret | oku')
    uret() if sys.argv[1] == 'uret' else oku()
