#!/usr/bin/env python3
# ga_ref.py : SilenitB F2b-1 - Ga referans fazlari (alfa-Ga, beta-Ga2O3): COD CIF -> QE vc-relax (90/360) -> scf (80/320) + okuma.
# CIF: 06_ref/cif/9008085.cif (alfa-Ga, Sharma & Donohue, Z. Kristallogr. 117 (1962) 293)
#      06_ref/cif/2004987.cif (beta-Ga2O3, Ahman, Svensson, Albertsson, Acta Cryst. C52 (1996) 1336)
# Kullanim (dizin /arf/scratch/nparan/SilenitB): python3 99_analiz/ga_ref.py uret | scf <aGa|bGa2O3> | oku
import sys, os, re, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g
B, RY, B2A = g.B, 13.605693, 0.529177
R = B + '/06_ref'
FAZ = {'aGa': dict(cif='9008085', nat=8, k='12 7 12', occ="occupations='smearing', smearing='mv', degauss=0.02", np=8),
       'bGa2O3': dict(cif='2004987', nat=20, k='3 10 5', occ="occupations='smearing', smearing='gaussian', degauss=0.005", np=4)}
SP = {'Ga': ('69.723', 'Ga_ONCV_PBE-1.0.upf'), 'O': ('15.999', 'O_ONCV_PBE-1.0.upf')}
E_O2 = -63.61024670   # Ry, Bi25FeO39/06_referans/O2 (SG15, 80/320, nspin=2; 12.11 teyit)

def sayi(s): return float(re.sub(r'\(.*\)', '', s))
def islem(op):                                   # '1/2+x' -> (katsayilar, sabit) ; eval kullanmaz
    out = []
    for e in op.replace(' ', '').lower().split(','):
        e = e if e[0] in '+-' else '+' + e; c = [0., 0., 0.]; t = 0.
        for tok in re.findall(r'[+-][^+-]+', e):
            sg = -1. if tok[0] == '-' else 1.; b = tok[1:]
            if b[-1] in 'xyz':
                k = b[:-1].rstrip('*'); c['xyz'.index(b[-1])] += sg * (eval_k(k) if k else 1.)
            else: t += sg * eval_k(b)
        out.append((c, t))
    return out
def eval_k(s):
    if not re.fullmatch(r'\d+(\.\d*)?(/\d+)?|\.\d+', s): g.dur('simetri ifadesi okunamadi: ' + s)
    p = s.split('/'); return float(p[0]) / (float(p[1]) if len(p) > 1 else 1.)

def cif(f):
    T = open(f, errors='ignore').read()
    cel = [sayi(re.search(r'_cell_%s\s+(\S+)' % k, T).group(1)) for k in ('length_a', 'length_b', 'length_c', 'angle_alpha', 'angle_beta', 'angle_gamma')]
    L = T.splitlines(); ops, ats = [], []
    i = 0
    while i < len(L):
        if L[i].strip() == 'loop_':
            h = []; i += 1
            while i < len(L) and L[i].strip().startswith('_'): h.append(L[i].strip().split()[0]); i += 1
            rows = []
            while i < len(L) and L[i].strip() and not L[i].strip().startswith(('_', 'loop_', '#')): rows.append(L[i].strip()); i += 1
            if any('symop_operation_xyz' in x or 'equiv_pos_as_xyz' in x for x in h):
                for r in rows:
                    m = re.findall(r"'([^']*)'", r); ops.append(m[0] if m else r.split()[-1])
            if '_atom_site_fract_x' in h:
                lab = h.index('_atom_site_type_symbol') if '_atom_site_type_symbol' in h else h.index('_atom_site_label')
                ix = [h.index('_atom_site_fract_%s' % c) for c in 'xyz']
                for r in rows:
                    w = r.split(); el = re.match(r'[A-Z][a-z]?', w[lab]).group(0)
                    ats.append((el, [sayi(w[j]) for j in ix]))
        else: i += 1
    if not ops or not ats: g.dur(f + ': simetri ya da atom dongusu yok')
    X = []
    for el, p in ats:
        for op in ops:
            q = [(sum(c[j] * p[j] for j in range(3)) + t) % 1. for c, t in islem(op)]
            if not any(e2 == el and all(min(abs(a - b), 1 - abs(a - b)) < 1e-3 for a, b in zip(q, p2)) for e2, p2 in X): X.append((el, q))
    return cel, X

def vek(a, b, c, al, be, ga):
    al, be, ga = (math.radians(x) for x in (al, be, ga))
    cx = c * math.cos(be); cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    return [[a, 0, 0], [b * math.cos(ga), b * math.sin(ga), 0], [cx, cy, math.sqrt(c * c - cx * cx - cy * cy)]]

def girdi(ad, calc, ec, V, X, F):
    els = sorted(set(e for e, _ in X), key=lambda e: ('Ga', 'O').index(e))
    s = "&control\n  calculation='%s', prefix='%s', outdir='./out', pseudo_dir='%s/00_pseudo'\n  tprnfor=.true., tstress=.true., forc_conv_thr=1.0d-4\n/\n" % (calc, ad, B)
    s += "&system\n  ibrav=0, nat=%d, ntyp=%d, ecutwfc=%.1f, ecutrho=%.1f, %s\n/\n&electrons\n  conv_thr=1.0d-10\n/\n" % (len(X), len(els), ec, 4 * ec, F['occ'])
    if calc == 'vc-relax': s += "&ions\n/\n&cell\n  cell_dofree='all', press_conv_thr=0.2\n/\n"
    s += 'ATOMIC_SPECIES\n' + ''.join('  %s %s %s\n' % (e, *SP[e]) for e in els)
    s += 'CELL_PARAMETERS angstrom\n' + ''.join('  %14.9f %14.9f %14.9f\n' % tuple(v) for v in V)
    s += 'ATOMIC_POSITIONS crystal\n' + ''.join('  %-3s %13.9f %13.9f %13.9f\n' % (e, *p) for e, p in X)
    return s + 'K_POINTS automatic\n  %s 0 0 0\n' % F['k']

def komsu(V, X, n=4):                            # her farkli ortam icin en yakin n komsu (bilesikte yalniz farkli element)
    tek = len(set(e for e, _ in X)) == 1; out = []
    for e, p in X:
        d = sorted(math.sqrt(sum(sum(((q[j] - p[j] + .5) % 1 - .5 + sh[j]) * V[j][k] for j in range(3)) ** 2 for k in range(3)))
                   for e2, q in X if tek or e2 != e for sh in [(i, jj, kk) for i in (-1, 0, 1) for jj in (-1, 0, 1) for kk in (-1, 0, 1)])
        d = tuple(round(x, 3) for x in d if x > 1e-3)[:n]
        if (e, d) not in out: out.append((e, d))
    return out

def uret():
    for e in SP: g.zval(SP[e][1])                 # pseudo 00_pseudo'da yoksa DUR
    if not os.path.exists(B + '/job_eos.sh') or 'temizlik' not in open(B + '/job_eos.sh').read(): print('(not: job_eos.sh temizliksiz; job_ref.sh temizlikli yazilir)')
    for f, F in FAZ.items():
        cel, X = cif('%s/cif/%s.cif' % (R, F['cif']))
        if len(X) != F['nat']: g.dur('%s: %d atom (beklenen %d)' % (f, len(X), F['nat']))
        V = vek(*cel)
        g.yaz('%s/%s_vc' % (R, f), girdi(f + '_vc', 'vc-relax', 90, V, X, F))
        print('%s: CIF %s a,b,c %.4f %.4f %.4f beta %.2f | %d atom | vc-relax 90/360 k %s' % (f, F['cif'], cel[0], cel[1], cel[2], cel[4], len(X), F['k']))
        for e, d in komsu(V, X): print('   %s komsu: %s' % (e, ' '.join('%.3f' % x for x in d)))
    open(B + '/job_ref.sh', 'w').write(g.HEAD.replace('@J@', 'gaB_ref').replace('@T@', '04:00:00') + JREF.replace('@B@', B))
    print('yazildi: 06_ref/{aGa,bGa2O3}_vc/run.in, job_ref.sh')

JREF = """cd @B@ || exit 1
for f in aGa bGa2O3; do
  np=$([ $f = aGa ] && echo 8 || echo 4)
  ( cd 06_ref/${f}_vc && [ ! -s run.out ] && mpirun pw.x -npool $np -in run.in > run.out )
  grep -q 'JOB DONE' 06_ref/${f}_vc/run.out && grep -q 'bfgs converged' 06_ref/${f}_vc/run.out || { echo "### $f vc-relax yakinsamadi -> scf yok"; continue; }
  python3 99_analiz/ga_ref.py scf $f || continue
  ( cd 06_ref/${f}_scf && mpirun pw.x -npool $np -in run.in > run.out )
  if grep -q 'JOB DONE' 06_ref/${f}_scf/run.out; then rm -f 06_ref/${f}_vc/out/*.wfc[0-9]* 06_ref/${f}_vc/out/*.save/wfc*.dat 06_ref/${f}_scf/out/*.wfc[0-9]* 06_ref/${f}_scf/out/*.save/wfc*.dat; echo "temizlik $f"; fi
done
"""

def son_yapi(t):
    k = t.rfind('Begin final coordinates')
    if k < 0: g.dur('Begin final coordinates yok')
    u = t[k:]; m = re.search(r'CELL_PARAMETERS\s*\(([^)]*)\)\s*\n(.*)\n(.*)\n(.*)\n', u)
    h = m.group(1).strip(); V = [[float(x) for x in m.group(i).split()] for i in (2, 3, 4)]
    f = float(h.split('=')[1]) * B2A if h.startswith('alat') else (B2A if h == 'bohr' else 1.)
    V = [[x * f for x in v] for v in V]; X = []
    for l in u[u.index('ATOMIC_POSITIONS'):].splitlines()[1:]:
        w = l.split()
        if len(w) < 4 or w[0] not in SP: break
        X.append((w[0], [float(q) for q in w[1:4]]))
    return V, X

def scf(f):
    F = FAZ[f]; V, X = son_yapi(open('%s/%s_vc/run.out' % (R, f)).read())
    if len(X) != F['nat']: g.dur('%s: son yapida %d atom' % (f, len(X)))
    g.yaz('%s/%s_scf' % (R, f), girdi(f + '_scf', 'scf', 80, V, X, F)); print('yazildi: 06_ref/%s_scf/run.in' % f)

def oku():
    E = {}
    for f, F in FAZ.items():
        rv, rs = g.oku1('%s/%s_vc/run.out' % (R, f)), g.oku1('%s/%s_scf/run.out' % (R, f))
        if rv is None: print(f, 'vc YOK'); continue
        cel, X0 = cif('%s/cif/%s.cif' % (R, F['cif']))
        V, X = son_yapi(rv['t']); ln = [math.sqrt(sum(x * x for x in v)) for v in V]
        be = math.degrees(math.acos(sum(a * b for a, b in zip(V[0], V[2])) / ln[0] / ln[2]))
        print('%s vc: JD %d bfgs %d adim %d P %.2f kbar | a,b,c %.4f %.4f %.4f (CIF %.4f %.4f %.4f; %+.2f/%+.2f/%+.2f %%) beta %.2f (CIF %.2f) | %s' % (
              f, rv['jd'], rv['bf'], rv['ns'], rv['P'], *ln, cel[0], cel[1], cel[2], *[100 * (x / y - 1) for x, y in zip(ln, cel[:3])], be, cel[4], rv['w']))
        for e, d in komsu(V, X): print('   %s komsu: %s' % (e, ' '.join('%.3f' % x for x in d)))
        if rs is None or not rs['jd']: print('   %s scf YOK/JD yok' % f); continue
        n = {e: sum(1 for x, _ in X if x == e) for e in SP}; E[f] = (rs['E'], n)
        print('   scf 80/320: E %.8f Ry  P %.2f kbar  | E/birim: %s' % (rs['E'], rs['P'], ('%.8f Ry/Ga' % (rs['E'] / n['Ga'])) if f == 'aGa' else ('%.8f Ry/Ga2O3' % (rs['E'] / (n['Ga'] / 2)))))
    if len(E) == 2:
        mGa = E['aGa'][0] / E['aGa'][1]['Ga']; eG = E['bGa2O3'][0] / (E['bGa2O3'][1]['Ga'] / 2)
        dH = (eG - 2 * mGa - 1.5 * E_O2) * RY
        print('DeltaH_f(beta-Ga2O3) = E(Ga2O3) - 2 mu(Ga) - 3/2 E(O2) = %.3f eV/f.u. (PBE, O2 duzeltmesiz; deney kiyasi kaynakla yapilacak)' % dH)

if __name__ == '__main__':
    a = sys.argv[1:]
    if a == ['uret']: uret()
    elif len(a) == 2 and a[0] == 'scf' and a[1] in FAZ: scf(a[1])
    elif a == ['oku']: oku()
    else: g.dur('kullanim: ga_ref.py uret | scf aGa|bGa2O3 | oku')
