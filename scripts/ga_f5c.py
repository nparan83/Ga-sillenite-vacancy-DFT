#!/usr/bin/env python3
# ga_f5c.py : SilenitB F5-c - karsi-olgusal Bi24Ga2O39'da eksik O(2) sitesi. Iki M de Ga iken P23 yorungeleri I23'te
# birlesir: (A,B) tet. 8c, (C,C') 24f, (D1,D2) O(2) 8c. 2A/2B/2C test edildi; O(2) edilmedi. 2D1 (f5_D1'den) ve 2D2
# (f5_D2'den, 2A/2B gibi yol testi), Bi_tet -> Ga; ayarlar kaynak run.in (relax 80/320 k2). Kullanim: python3 99_analiz/ga_f5c.py uret | oku
# Onceden yazilan (26.11): K1 JD+bfgs; K5a |E(2D1)-E(2D2)| <= 5 meV; K5b D = min(2D) - E(2A):
#   D >= -10 meV -> karsi-olgusalda A/B/C'nin altinda site yok, iddia yazildigi gibi; -100 <= D < -10 meV -> zayif secicilik,
#   ifade revize; D < -100 meV -> DUR (karsi-olgusal O(2)'yi seciyor; baslik/iddia yeniden kurulur).
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
import ga_f2c as gc
B, RY = g.B, 13.605693
gc.KAY = {'D1': 'f5_D1', 'D2': 'f5_D2'}

def oku():
    R = {}
    for k in ('A', 'B', 'C', 'D1', 'D2'):
        r = g.oku1('%s/03_site/ga2%s_rel/run.out' % (B, k))
        if r is None: print('ga2%s_rel YOK' % k); continue
        print('ga2%s_rel: JD %d bfgs %d adim %d F %.5f P %+.2f kbar | %s' % (k, r['jd'], r['bf'], r['ns'], r['F'] or 0, r['P'] or 0, r['w']))
        if r['jd'] and r['bf']: R[k] = r['E']
    if not {'A', 'D1', 'D2'} <= set(R): print('### K1: ga2A, ga2D1, ga2D2 JD + bfgs gerekir'); return
    print('dE (ga2A = 0, meV): ' + ' '.join('%s %+.1f' % (k, (R[k] - R['A']) * RY * 1000) for k in sorted(R)))
    p = abs(R['D1'] - R['D2']) * RY * 1000; D = (min(R['D1'], R['D2']) - R['A']) * RY * 1000
    print('K5a yol testi |2D1 - 2D2| = %.1f meV -> %s' % (p, 'gecti' if p <= 5 else 'YOL BAGIMLI (dusuk olan kullanilir, metinde belirtilir)'))
    print('K5b D = min(2D) - 2A = %+.1f meV -> %s' % (D, 'A/B/C altinda site yok; iddia yazildigi gibi' if D >= -10 else
          ('zayif secicilik; ifade revize' if D >= -100 else 'DUR: karsi-olgusal O(2) seciyor; iddia yeniden kurulur')))

if __name__ == '__main__':
    m = sys.argv[1] if len(sys.argv) == 2 else ''
    if m not in ('uret', 'oku'): g.dur('kullanim: ga_f5c.py uret | oku')
    gc.uret() if m == 'uret' else oku()
