#!/usr/bin/env python3
# ga_f4b.py : SilenitB F4b - Ga merkez-disi kaymasi (off-centering) testi. B yapisinda hedef hacmi, KAYNAK hacmin son
# konumlarindan (Ga -0.104 A kaymis, s=1.010) yeniden relaks eder; orijinal EOS noktasiyla enerji/geometri kiyasi.
# Kullanim (dizin /arf/scratch/nparan/SilenitB): python3 99_analiz/ga_f4b.py uret | oku
# Is: sbatch -J h1.000 job_eos.sh gaB_s1.000_f1.010 ; sbatch -J h0.990 job_eos.sh gaB_s0.990_f1.010
import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
import ga_f2 as g2
B, RY, A0, B2A = g.B, 13.605693, 19.58225862, 0.529177
CIFT = ((1.000, 1.010), (0.990, 1.010))      # (hedef s, kaynak s)
def ad(s, k=None): return 'gaB_s%.3f' % s + ('' if k is None else '_f%.3f' % k)

def son(n):
    r = g.oku1('%s/04_eos/%s/run.out' % (B, n))
    if r is None or not (r['jd'] and r['bf']): g.dur(n + ': run.out yok ya da JD/bfgs yok (K1)')
    return r

def uret():
    for s, k in CIFT:
        X = g2.konum(son(ad(k))['t'])
        if len(X) != 65: g.dur('%s: %d atom' % (ad(k), len(X)))
        T = open('%s/04_eos/%s/run.in' % (B, ad(s))).read(); L = T.splitlines()
        i = [q for q, l in enumerate(L) if l.strip().upper().startswith('ATOMIC_POSITIONS')]
        if len(i) != 1 or 'crystal' not in L[i[0]].lower(): g.dur('ATOMIC_POSITIONS (crystal) tek degil')
        i = i[0]; j = i + 1
        while j < len(L) and L[j].split() and L[j].split()[0] in ('Bi', 'Ga', 'O'): j += 1
        if j - i - 1 != 65: g.dur('sablonda %d konum' % (j - i - 1))
        t = '\n'.join(L[:i + 1] + ['%-3s %14.10f %14.10f %14.10f' % (e, *x) for e, x in X] + L[j:]) + '\n'
        t = g.setk(t, 'prefix', "'%s'" % ad(s, k))
        q = g2.geo(X, A0 * s * B2A)
        g.yaz('%s/04_eos/%s' % (B, ad(s, k)), t)
        print('%s: hedef a %.5f A, baslangic Ga ||[111] %+.3f, Bi_tet %+.3f | yazildi' % (ad(s, k), A0 * s * B2A, q['gak'], q['kay']))

def oku():
    for s, k in CIFT:
        r = g.oku1('%s/04_eos/%s/run.out' % (B, ad(s, k))); r0 = son(ad(s)); a = A0 * s * B2A
        if r is None: print(ad(s, k), 'YOK'); continue
        q, q0 = g2.geo(g2.konum(r['t']), a), g2.geo(g2.konum(r0['t']), a)
        dE = (r['E'] - r0['E']) * RY * 1000
        print('%s: JD %d bfgs %d adim %d P %.2f GPa wall %s' % (ad(s, k), r['jd'], r['bf'], r['ns'], r['P'] / 10, r['w']))
        print('   Ga ||[111] %+.3f (orijinal %+.3f) | Bi_tet %+.3f (%+.3f) | Ga-O %s' % (q['gak'], q0['gak'], q['kay'], q0['kay'], ' '.join('%.3f' % d for d in q['gaO'][:4])))
        print('   E(yeni) - E(orijinal EOS noktasi) = %+.2f meV | P(orijinal) %.2f GPa' % (dE, r0['P'] / 10))
        if not (r['jd'] and r['bf']): print('   ### JD/bfgs yok -> KARAR VERILMEZ (K1)'); continue
        if dE < -1 and abs(q['gak']) >= 0.05: print('   KARAR H-OFF: merkez-disi Ga bu hacimde DAHA DUSUK -> EOS kolu ve F2 gaB_rel yeniden')
        elif abs(dE) <= 1 and abs(q['gak']) < 0.05: print('   KARAR H-GERI: merkeze dondu -> kayma bu hacmin USTUNDE basliyor')
        else: print('   KARAR H-ARA: elle incelenecek')

if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] not in ('uret', 'oku'): g.dur('kullanim: ga_f4b.py uret | oku')
    uret() if sys.argv[1] == 'uret' else oku()
