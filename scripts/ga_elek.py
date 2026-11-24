#!/usr/bin/env python3
# ga_elek.py : SilenitB F2b-2 - Ga silenit elektronik yapi (S = Bi25GaO40, B = Bi25GaO39; SG15 relakse geometriler).
#  (1) SG15 nscf k 4x4x4 (relaks yuk yogunlugundan) -> bant sayimiyla gap (VBM = bant N/2, CBM = bant N/2+1).
#  (2) PAW scf 60/480 k 2x2x2 (Makale A 2.5: pslibrary 1.0.0, tek nokta) + projwfc (degauss 0.01 Ry, DeltaE 0.02 eV) -> PDOS.
# Kullanim (dizin /arf/scratch/nparan/SilenitB): python3 99_analiz/ga_elek.py uret | oku     Is: sbatch -J el_S job_elek.sh S ; sbatch -J el_B job_elek.sh B
import sys, os, re, glob, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
import ga_f2 as g2
B, F, HA, B2A = g.B, g.F, 27.211386, 0.529177
E5 = B + '/05_elek'
YAP = {'S': ('02_O40/gaS_rel', 'gaS_rel', 66), 'B': ('03_site/gaB_rel', 'gaB_rel', 65)}
PAW = {'Bi': 'Bi.pbe-dn-kjpaw_psl.1.0.0.UPF', 'Ga': 'Ga.pbe-dnl-kjpaw_psl.1.0.0.UPF', 'O': 'O.pbe-n-kjpaw_psl.1.0.0.UPF'}
A0 = 19.58225862 * B2A

def kaynak(y):
    d, p, n = YAP[y]; r = g.oku1('%s/%s/run.out' % (B, d))
    if r is None or not (r['jd'] and r['bf']): g.dur(d + ': JD/bfgs yok (K1)')
    X = g2.konum(r['t'])
    if len(X) != n: g.dur('%s: %d atom' % (d, len(X)))
    return open('%s/%s/run.in' % (B, d)).read(), X, '%s/%s/out/%s.save' % (B, d, p), p

def konum_koy(T, X):
    L = T.splitlines(); i = [q for q, l in enumerate(L) if l.strip().upper().startswith('ATOMIC_POSITIONS')][0]; j = i + 1
    while j < len(L) and L[j].split() and L[j].split()[0] in ('Bi', 'Ga', 'O'): j += 1
    return '\n'.join(L[:i + 1] + ['%-3s %14.10f %14.10f %14.10f' % (e, *x) for e, x in X] + L[j:]) + '\n'

def kesme(f):
    t = open(f, errors='ignore').read(); m = re.findall(r'Suggested minimum cutoff for (wavefunctions|charge density):\s*([\d.]+)', t)
    return {k[:4]: float(v) for k, v in m}

def kprep(t, k):
    m = re.search(r'(?im)^[ \t]*K_POINTS[^\n]*\n[ \t]*\d+[ \t]+\d+[ \t]+\d+[ \t]+(\d+[ \t]+\d+[ \t]+\d+)', t)
    return t[:m.start()] + 'K_POINTS automatic\n  %d %d %d %s' % (k, k, k, m.group(1)) + t[m.end():]

def uret():
    for e, u in PAW.items():
        h = '%s/00_pseudo/%s' % (B, u)
        if not os.path.exists(h):
            s = '%s/00_pseudo/%s' % (F, u)
            if not os.path.exists(s): g.dur('PAW pseudo yok: ' + u)
            shutil.copy(s, h); print('kopyalandi:', u)
        c = kesme(h); print('%-32s onerilen kesme: %s' % (u, c))
        if c.get('wave', 0) > 60 or c.get('char', 0) > 480: g.dur(u + ': onerilen kesme 60/480 ustunde')
    for y in YAP:
        T, X, sv, p = kaynak(y); ne = int(round(sum(g.zval(re.search(r'(?im)^[ \t]*%s[ \t]+[0-9.]+[ \t]+(\S+\.upf)[ \t]*$' % e, T).group(1)) * sum(1 for x, _ in X if x == e) for e in ('Bi', 'Ga', 'O'))))
        if not any(os.path.exists(sv + '/charge-density.' + x) for x in ('dat', 'hdf5')): g.dur(sv + ': charge-density.{dat,hdf5} yok')
        d = '%s/ga%s_nscf' % (E5, y); dst = '%s/out/%s.save' % (d, p)
        if os.path.exists(d + '/run.out'): g.dur(d + ': run.out var')
        if not os.path.exists(dst): shutil.copytree(sv, dst, ignore=shutil.ignore_patterns('wfc*'))
        t = re.sub(r'(?is)&(ions|cell)\b[^/]*/[ \t]*\n?', '', konum_koy(T, X))
        t = g.setk(g.setk(g.setk(g.setk(g.setk(t, 'calculation', "'nscf'"), 'outdir', "'./out'"), 'tprnfor', '.false.'), 'tstress', '.false.'), 'prefix', "'%s'" % p)
        t = re.sub(r'(?i)(&system)', r'\1\n  nbnd = %d,' % (ne // 2 + 24), t, count=1)
        g.yaz(d, kprep(t, 4))
        t = konum_koy(T, X); t = re.sub(r'(?is)&(ions|cell)\b[^/]*/[ \t]*\n?', '', t)
        for k, v in (('calculation', "'scf'"), ('prefix', "'ga%s_paw'" % y), ('outdir', "'./out'"), ('ecutwfc', '60.0'), ('ecutrho', '480.0'), ('tprnfor', '.false.'), ('tstress', '.false.')): t = g.setk(t, k, v)
        t = re.sub(r'(?i)(&system)', r'\1\n  nbnd = %d,' % (ne // 2 + 40), t, count=1)
        for e, u in PAW.items():
            t, c = re.subn(r'(?im)^([ \t]*%s[ \t]+[0-9.]+[ \t]+)\S+\.upf[ \t]*$' % e, r'\g<1>' + u, t)
            if c != 1: g.dur('ATOMIC_SPECIES %s %d kez' % (e, c))
        g.yaz('%s/ga%s_paw' % (E5, y), kprep(t, 2))
        open('%s/ga%s_paw/proj.in' % (E5, y), 'w').write("&projwfc\n  prefix='ga%s_paw', outdir='./out', degauss=0.01, DeltaE=0.02, filpdos='pdos/ga%s'\n/\n" % (y, y))
        q = g2.geo(X, A0); print('%s: N = %d, nscf nbnd %d (k4), PAW nbnd %d (k2), Bi_tet atom #%d | yazildi: 05_elek/ga%s_{nscf,paw}' % (y, ne, ne // 2 + 24, ne // 2 + 40, q['bt'] + 1, y))
    open(B + '/job_elek.sh', 'w').write(g.HEAD.replace('@J@', 'gaB_elek').replace('@T@', '06:00:00') + JOB.replace('@B@', B))
    print('yazildi: job_elek.sh (temizlikli)')

JOB = """y=$1; E=@B@/05_elek
cd $E/ga${y}_nscf || exit 1
[ -s run.out ] || mpirun pw.x -npool 4 -in run.in > run.out
grep -q 'JOB DONE' run.out && rm -f out/*.wfc[0-9]* out/*.save/wfc* && echo "temizlik nscf $y"
cd $E/ga${y}_paw || exit 1
[ -s run.out ] || mpirun pw.x -npool 2 -in run.in > run.out
grep -q 'JOB DONE' run.out || { echo "### PAW scf $y bitmedi -> projwfc yok"; exit 1; }
mkdir -p pdos; mpirun projwfc.x -npool 2 -in proj.in > proj.out
grep -q 'JOB DONE' proj.out && rm -f out/*.wfc[0-9]* out/*.save/wfc* && echo "temizlik paw $y"
"""

def xml_bant(f):
    t = open(f, errors='ignore').read()
    ne = float(re.search(r'<nelec>\s*([\d.eE+-]+)', t).group(1))
    K = [[float(v) * HA for v in m.split()] for m in re.findall(r'<eigenvalues[^>]*>(.*?)</eigenvalues>', t, re.S)]
    n = int(round(ne)) // 2
    vb = max(range(len(K)), key=lambda i: K[i][n - 1]); cb = min(range(len(K)), key=lambda i: K[i][n])
    return dict(ne=ne, nk=len(K), vbm=K[vb][n - 1], cbm=K[cb][n], gap=K[cb][n] - K[vb][n - 1], kv=vb, kc=cb,
                dir=min(k[n] - k[n - 1] for k in K))

def pdos(d, y, vbm, cbm, bt):
    T = {}
    for f in glob.glob('%s/pdos/ga%s.pdos_atm#*' % (d, y)):
        m = re.search(r'atm#(\d+)\((\w+)\)_wfc#\d+\((\w)', f); a, e, l = int(m.group(1)), m.group(2), m.group(3)
        key = ('Bi_tet' if a == bt else e) + ' ' + l
        for row in open(f):
            if row.startswith('#'): continue
            w = row.split(); E, v = round(float(w[0]), 3), float(w[1]); T.setdefault(E, {}); T[E][key] = T[E].get(key, 0.) + v
    if not T: return None
    ks = sorted({k for v in T.values() for k in v}); Es = sorted(T)
    def pay(a, b):
        s = {k: sum(T[E].get(k, 0.) for E in Es if a <= E <= b) for k in ks}; z = sum(s.values()) or 1.
        return ' '.join('%s %.0f%%' % (k, 100 * s[k] / z) for k in sorted(s, key=lambda k: -s[k])[:6])
    with open('%s/pdos_ga%s.dat' % (E5, y), 'w') as o:
        o.write('# E-VBM(eV) ' + ' '.join(k.replace(' ', '_') for k in ks) + '\n')
        for E in Es: o.write('%.3f ' % (E - vbm) + ' '.join('%.5f' % T[E].get(k, 0.) for k in ks) + '\n')
    return pay(vbm - 1., vbm), pay(cbm, cbm + 1.)

def oku():
    for y in YAP:
        T, X, sv, p = kaynak(y); q = g2.geo(X, A0)
        rn, rp = g.oku1('%s/ga%s_nscf/run.out' % (E5, y)), g.oku1('%s/ga%s_paw/run.out' % (E5, y))
        print('== %s (%s) ==' % (y, 'Bi25GaO40' if y == 'S' else 'Bi25GaO39'))
        if rn and rn['jd']:
            b = xml_bant('%s/ga%s_nscf/out/%s.save/data-file-schema.xml' % (E5, y, p))
            print('  SG15 nscf k4: N %.0f, %d k | VBM %.3f CBM %.3f eV | gap %.3f eV (k#%d -> k#%d) | en kucuk dogrudan %.3f eV | %s' % (b['ne'], b['nk'], b['vbm'], b['cbm'], b['gap'], b['kv'] + 1, b['kc'] + 1, b['dir'], rn['w']))
        else: print('  nscf YOK/JD yok')
        if not (rp and rp['jd']): print('  PAW scf YOK/JD yok'); continue
        b = xml_bant('%s/ga%s_paw/out/ga%s_paw.save/data-file-schema.xml' % (E5, y, y))
        po = open('%s/ga%s_paw/proj.out' % (E5, y), errors='ignore').read() if os.path.exists('%s/ga%s_paw/proj.out' % (E5, y)) else ''
        sp = re.findall(r'Spilling Parameter:\s*([\d.]+)', po)
        print('  PAW scf k2: gap %.3f eV | E %.6f Ry | spilling %s | %s' % (b['gap'], rp['E'], sp[-1] if sp else '-', rp['w']))
        r = pdos('%s/ga%s_paw' % (E5, y), y, b['vbm'], b['cbm'], q['bt'] + 1)
        if r: print('  VB tepesi [VBM-1, VBM]: %s\n  CB tabani [CBM, CBM+1]: %s\n  yazildi: 05_elek/pdos_ga%s.dat' % (r[0], r[1], y))
        else: print('  pdos dosyasi yok')

if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] not in ('uret', 'oku'): g.dur('kullanim: ga_elek.py uret | oku')
    uret() if sys.argv[1] == 'uret' else oku()
