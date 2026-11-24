#!/usr/bin/env python3
# ga_f4_tani.py : F4 EOS tanisi (0 c-saat, login dugumu). Her EOS kosumu icin relaks yolu: adim sayisi, relaksasyon enerjisi,
# Bi_tet kaymasi ve Ga kaymasi (ilk -> son, her ~10 adimda), son Bi_tet-O ve Ga-O mesafeleri.
# Kullanim (dizin /arf/scratch/nparan/SilenitB): python3 99_analiz/ga_f4_tani.py
import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
import ga_f2 as g2
B, RY, A0, B2A = g.B, 13.605693, 19.58225862, 0.529177
OL = (0.94, 0.96, 0.98, 0.99, 1.00, 1.01, 1.02)

def bloklar(t):
    X = []
    for m in re.finditer(r'ATOMIC_POSITIONS[^\n]*\n', t):
        b = []
        for l in t[m.end():].splitlines():
            w = l.split()
            if len(w) < 4 or w[0] not in ('Bi', 'Ga', 'O'): break
            b.append((w[0], [float(q) for q in w[1:4]]))
        X.append(b)
    return X

for y in ('B', 'S'):
    for s in OL:
        n = 'ga%s_s%.3f' % (y, s); f = '%s/04_eos/%s/run.out' % (B, n)
        if not os.path.exists(f): print(n, 'YOK'); continue
        t = open(f, errors='ignore').read(); a = A0 * s * B2A
        E = [float(x) * RY for x in re.findall(r'^!\s+total energy\s+=\s+(-?[\d.]+)', t, re.M)]
        X = [b for b in bloklar(t) if len(b) in (65, 66)]
        T = [open('%s/04_eos/%s/run.in' % (B, n)).read()]
        X0 = bloklar(T[0])[0]; Q = [g2.geo(b, a) for b in [X0] + X]
        iz = ' '.join('%d:%+.3f/%+.3f' % (k, Q[k]['kay'], Q[k]['gak']) for k in sorted(set(list(range(0, len(Q), 10)) + [len(Q) - 1])))
        print('%s adim %d  E_ilk-E_son %.1f meV | Bi_tet/Ga kaymasi (adim:Bt/Ga) %s' % (n, len(E), (E[0] - E[-1]) * 1000 if E else float('nan'), iz))
        print('   son Bi_tet-O %s | Ga-O %s' % (' '.join('%.3f' % d for d in Q[-1]['btO']), ' '.join('%.3f' % d for d in Q[-1]['gaO'])))
