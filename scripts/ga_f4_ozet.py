#!/usr/bin/env python3
# ga_f4_ozet.py : F4 nihai EOS ozeti (0 c-saat, login dugumu). B ve S icin: noktalar (her hacimde en dusuk enerjili kosum), yerel B,
# -dE/dV <-> P tutarliligi, BM3 fitleri (7 nokta / sikistirma kolu 0.94-1.00 / B icin +1.005 duyarlilik). Veri: 99_analiz/eos_ozet.dat
# Kullanim (dizin /arf/scratch/nparan/SilenitB): python3 99_analiz/ga_f4_ozet.py | tee 99_analiz/ga_f4_ozet.log
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
import ga_f2 as g2
import ga_f4 as g4
B, RY, A0, B2A, GPA = g.B, 13.605693, 19.58225862, 0.529177, 160.21766
OL7 = (0.94, 0.96, 0.98, 0.99, 1.00, 1.01, 1.02)

def nokta(n, s):
    r = g.oku1('%s/04_eos/%s/run.out' % (B, n))
    if r is None or not (r['jd'] and r['bf']): return None
    q = g2.geo(g2.konum(r['t']), A0 * s * B2A)
    return dict(s=s, V=(A0 * s * B2A) ** 3, E=r['E'] * RY, P=r['P'] / 10, ga=q['gak'], bt=q['kay'], ad=n)

def kume(y):
    K = {}
    for s in OL7:
        ad = ['ga%s_s%.3f' % (y, s)] + (['gaB_s%.3f_f1.010' % s] if y == 'B' and s in (0.99, 1.0) else [])
        X = [p for p in (nokta(n, s) for n in ad) if p]
        if X: K[s] = min(X, key=lambda p: p['E'])
    if y == 'B':
        X = [p for p in (nokta('gaB_s1.005_%s' % k, 1.005) for k in 'co') if p]
        if X: K[1.005] = min(X, key=lambda p: p['E'])
    return [K[s] for s in sorted(K)]

def fit(L, lab):
    if len(L) < 5: print('  %-30s n=%d < 5 -> fit yok' % (lab, len(L))); return
    f = g4.bm3([p['V'] for p in L], [p['E'] for p in L])
    if f is None: print('  %-30s BM3 minimum yok' % lab); return
    dP = max(abs(f['P'](p['V']) - p['P']) for p in L)
    print('  %-30s n=%d V0 %.2f a0 %.4f | B0 %.1f GPa B0p %.2f | artik %.2f meV | maks|dP| %.2f GPa | P %.2f..%.2f' % (
          lab, len(L), f['V0'], f['a0'], f['B0'], f['Bp'], f['res'], dP, min(p['P'] for p in L), max(p['P'] for p in L)))

o = open(B + '/99_analiz/eos_ozet.dat', 'w'); o.write('# yapi s V_A3 E_eV P_GPa Ga111_A Bitet111_A kosum\n')
for y in 'BS':
    L = kume(y); E1 = [p['E'] for p in L if p['s'] == 1.0]; E1 = E1[0] if E1 else L[0]['E']
    print('== %s (%s) ==' % (y, 'Bi25GaO39' if y == 'B' else 'Bi25GaO40'))
    for p in L:
        print('  s %.3f V %8.2f E-E1 %8.4f eV P %6.2f GPa Ga %+.3f Bt %+.3f  %s' % (p['s'], p['V'], p['E'] - E1, p['P'], p['ga'], p['bt'], p['ad']))
        o.write('%s %.4f %.4f %.6f %.4f %.4f %.4f %s\n' % (y, p['s'], p['V'], p['E'], p['P'], p['ga'], p['bt'], p['ad']))
    for a, b in zip(L, L[1:]):
        Vm = (a['V'] + b['V']) / 2
        print('  %.3f-%.3f yerel B %5.1f GPa | P_ort %6.2f | -dE/dV %6.2f GPa' % (a['s'], b['s'], -Vm * (b['P'] - a['P']) / (b['V'] - a['V']),
              (a['P'] + b['P']) / 2, -(b['E'] - a['E']) / (b['V'] - a['V']) * GPA))
    fit([p for p in L if p['s'] in OL7], '7 nokta (s=1.00 en dusuk)')
    fit([p for p in L if p['s'] <= 1.0], 'sikistirma kolu 0.94-1.00')
    if y == 'B': fit([p for p in L if p['s'] <= 1.005], 'duyarlilik: 0.94-1.005')
print('yazildi: 99_analiz/eos_ozet.dat')
