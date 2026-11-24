#!/usr/bin/env python3
# ga_hse.py : SilenitB F2b-4 - HSE06 saglamasi (S = Bi25GaO40, B = Bi25GaO39; SG15 relakse geometriler, 80/320, ecutfock 160, nq 1x1x1, k 2x2x2).
# Ayni kosumun ilk (PBE) ve son (HSE) "highest occupied, lowest unoccupied level" satirlarindan gap; k4 PBE gap'ine kayma (makas) olarak uygulanir.
# Kullanim (dizin /arf/scratch/nparan/SilenitB): python3 99_analiz/ga_hse.py uret S [B] | oku      Is: sbatch -J hse_S job_hse.sh S
import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
import ga_elek as ge
B = g.B
H7 = B + '/07_hse'
DENEY = (2.68, 2.99)          # eV: Wang 2012 (2.68), Kargin derlemesi (2.99); Lin 2008 2.95 (durum 2b)

def uret(sec):
    for y in sec:
        T, X, sv, p = ge.kaynak(y)
        ne = int(round(sum(g.zval(re.search(r'(?im)^[ \t]*%s[ \t]+[0-9.]+[ \t]+(\S+\.upf)[ \t]*$' % e, T).group(1)) * sum(1 for x, _ in X if x == e) for e in ('Bi', 'Ga', 'O'))))
        t = re.sub(r'(?is)&(ions|cell)\b[^/]*/[ \t]*\n?', '', ge.konum_koy(T, X))
        for k, v in (('calculation', "'scf'"), ('prefix', "'ga%s_hse'" % y), ('outdir', "'./out'"), ('tprnfor', '.false.'), ('tstress', '.false.'), ('occupations', "'fixed'")): t = g.setk(t, k, v)
        t = re.sub(r'(?i)[ \t]*\b(smearing|degauss)\s*=\s*[^,\n/]+,?', '', t)
        t = re.sub(r'(?i)(&system)', r"\1\n  nbnd = %d, input_dft = 'hse', ecutfock = 160.0, nqx1 = 1, nqx2 = 1, nqx3 = 1," % (ne // 2 + 8), t, count=1)
        if re.search(r'(?i)smearing|degauss', t): g.dur('smearing kaldirilamadi')
        g.yaz('%s/ga%s_hse' % (H7, y), ge.kprep(t, 2))
        print('%s: N %d, nbnd %d, HSE06 80/320 ecutfock 160, nq 1, k 2x2x2, occupations fixed | yazildi: 07_hse/ga%s_hse/run.in' % (y, ne, ne // 2 + 8, y))
    open(B + '/job_hse.sh', 'w').write(g.HEAD.replace('@J@', 'gaB_hse').replace('@T@', '16:00:00') + JOB.replace('@B@', B))
    print('yazildi: job_hse.sh (16 sa, temizlikli)')

JOB = """cd @B@/07_hse/ga${1}_hse || exit 1
[ -s run.out ] && { echo "### run.out zaten var -> DUR"; exit 1; }
mpirun pw.x -npool 4 -in run.in > run.out
grep -q 'JOB DONE' run.out && rm -f out/*.wfc[0-9]* out/*.save/wfc* && echo "temizlik hse $1"
"""

def oku():
    R = {}
    for y in ('S', 'B'):
        r = g.oku1('%s/ga%s_hse/run.out' % (H7, y))
        if r is None: print('ga%s_hse YOK' % y); continue
        L = [tuple(map(float, m)) for m in re.findall(r'highest occupied, lowest unoccupied level \(ev\):\s*(-?[\d.]+)\s+(-?[\d.]+)', r['t'])]
        dx = re.findall(r'est\. exchange err \(dexx\)\s*=\s*([\d.eE+-]+)', r['t'])
        nc = 'convergence NOT achieved' in r['t']
        print('ga%s_hse: JD %d | %d scf blogu | son dexx %s | NOT achieved %s | E %.8f Ry | %s' % (y, r['jd'], len(L), dx[-1] if dx else '-', nc, r['E'], r['w']))
        if not (r['jd'] and len(L) >= 2 and not nc): print('   ### tamamlanmamis -> KARAR VERILMEZ (K1)'); continue
        gp, gh = L[0][1] - L[0][0], L[-1][1] - L[-1][0]
        b = ge.xml_bant('%s/ga%s_nscf/out/%s.save/data-file-schema.xml' % (ge.E5, y, ge.YAP[y][1]))
        R[y] = (gp, gh, b['gap'] + gh - gp)
        print('   k2: PBE gap %.3f eV (PAW k2: S 0.930 / B 2.239) | HSE gap %.3f eV | kayma %+.3f | PBE k4 %.3f -> HSE (makas) %.3f eV' % (gp, gh, gh - gp, b['gap'], R[y][2]))
    if len(R) == 2:
        d = R['B'][2] - R['S'][2]; uz = lambda x: 0. if DENEY[0] <= x <= DENEY[1] else min(abs(x - DENEY[0]), abs(x - DENEY[1]))
        print('HSE: gap(B) - gap(S) = %.3f eV | deney araligina uzaklik B %.2f, S %.2f eV' % (d, uz(R['B'][2]), uz(R['S'][2])))
        print('KARAR: %s' % ('H-KORUNUR (fark >= 1.0 eV ve B deneye daha yakin)' if d >= 1.0 and uz(R['B'][2]) < uz(R['S'][2]) else 'H-ZAYIF (ayrim PBE e bagli olabilir)'))

if __name__ == '__main__':
    a = sys.argv[1:]
    if len(a) >= 2 and a[0] == 'uret' and set(a[1:]) <= {'S', 'B'}: uret(a[1:])
    elif a == ['oku']: oku()
    else: g.dur('kullanim: ga_hse.py uret S [B] | oku')
