#!/usr/bin/env python3
# ga_f1.py : SilenitB F1b/F1 - Ga silenit girdi uretici + okuyucu (SG15 ONCV, nspin=1)
# Kullanim (dizin /arf/scratch/nparan/SilenitB):  python3 99_analiz/ga_f1.py uret | oku
# Sablon: Bi25FeO39/03_saf_relax/safrel.in ; konumlar: 08_paw_U/pos{Saf,A,B}_relax.pos (Fe -> Ga)
import sys, re, os
B = os.environ.get('SILB', '/arf/scratch/nparan/SilenitB')
F = os.environ.get('SILF', '/arf/scratch/nparan/Bi25FeO39')
TPL, PSD = F + '/03_saf_relax/safrel.in', B + '/00_pseudo'
POS = {'S': 'Saf', 'A': 'A', 'B': 'B'}
GA = 'Ga_ONCV_PBE-1.0.upf'
# (yapi, ecutwfc Ry, k) ; ecutrho = 4 x ecutwfc ; sira = is icindeki kosum sirasi
YAK = [('S', 80, 2), ('S', 70, 2), ('S', 90, 2), ('S', 100, 2), ('S', 80, 1), ('S', 80, 3),
       ('B', 80, 2), ('A', 80, 2), ('B', 100, 2), ('A', 100, 2)]
RY = 13.605693

def dur(m): sys.exit('### ' + m + ' -> DUR')
def ad(y, e, k): return '%s_e%03dk%d' % (y, e, k)

def setk(t, key, val):
    pat = r"(\b" + re.escape(key) + r"\s*=\s*)('[^']*'|[^\s,/!]+)"
    t2, n = re.subn(pat, lambda m: m.group(1) + val, t, count=1, flags=re.I)
    if n != 1: dur('sablonda anahtar yok: ' + key)
    return t2

def zval(f):
    p = os.path.join(PSD, f)
    if not os.path.exists(p): dur('pseudo yok: ' + p)
    m = re.search(r'z_valence\s*=\s*"\s*([0-9.eE+-]+)', open(p, errors='ignore').read())
    if not m: dur('z_valence okunamadi: ' + f)
    return float(m.group(1))

def uret1(T, name, yap, calc, ec, k):
    P = []
    for l in open(F + '/08_paw_U/pos%s_relax.pos' % POS[yap]):
        w = l.split()
        if len(w) >= 4 and w[0] in ('Bi', 'Fe', 'O'):
            P.append('%-3s %14.10f %14.10f %14.10f' % ('Ga' if w[0] == 'Fe' else w[0], *map(float, w[1:4])))
    n = {e: sum(1 for l in P if l.split()[0] == e) for e in ('Bi', 'Ga', 'O')}
    if n['Bi'] != 25 or n['Ga'] != 1 or n['O'] != (40 if yap == 'S' else 39): dur('%s atom sayimi %s' % (yap, n))
    t = T
    for p in (r'starting_magnetization\s*\(\s*\d+\s*\)\s*=\s*[-+0-9.dDeE]+', r'tot_magnetization\s*=\s*[-+0-9.dDeE]+',
              r'nspin\s*=\s*\d+'):
        t = re.sub(p + r'[ \t]*,?', '', t, flags=re.I)
    t, c = re.subn(r'(?im)^[ \t]*Fe[ \t]+[0-9.]+[ \t]+\S+[ \t]*$', '  Ga   69.723  ' + GA, t)
    if c != 1: dur('ATOMIC_SPECIES Fe satiri %d kez bulundu' % c)
    L = t.splitlines()
    i = [q for q, l in enumerate(L) if l.strip().upper().startswith('ATOMIC_POSITIONS')]
    if len(i) != 1 or 'crystal' not in L[i[0]].lower(): dur('ATOMIC_POSITIONS (crystal) tek degil/yok')
    i = i[0]; j = i + 1
    while j < len(L) and L[j].split() and L[j].split()[0] in ('Bi', 'Fe', 'O'): j += 1
    if j - i - 1 < 60: dur('sablon konum blogu okunamadi (%d satir)' % (j - i - 1))
    t = '\n'.join(L[:i + 1] + P + L[j:]) + '\n'
    m = re.search(r'(?im)^[ \t]*K_POINTS[^\n]*automatic[^\n]*\n[ \t]*(\d+)[ \t]+(\d+)[ \t]+(\d+)[ \t]+(\d+)[ \t]+(\d+)[ \t]+(\d+)', t)
    if not m: dur('K_POINTS automatic yok')
    kay = ' '.join(m.group(4, 5, 6))
    t = t[:m.start(1)] + '%d %d %d %s' % (k, k, k, kay) + t[m.end(6):]
    for key, val in (('calculation', "'%s'" % calc), ('prefix', "'%s'" % name), ('outdir', "'./out'"),
                     ('pseudo_dir', "'%s'" % PSD), ('ecutwfc', '%d.0' % ec), ('ecutrho', '%d.0' % (4 * ec)),
                     ('nat', str(len(P)))):
        t = setk(t, key, val)
    for key in ('tprnfor', 'tstress'):
        if re.search(r'\b' + key + r'\s*=', t, re.I): t = setk(t, key, '.true.')
        else: t = re.sub(r'(?i)(&control)', r'\1\n  %s = .true.\n' % key, t, count=1)
    if calc == 'scf': t = re.sub(r'(?is)&(ions|cell)\b[^/]*/[ \t]*\n?', '', t)
    elif not re.search(r'(?i)&ions', t): dur('relax icin &IONS yok')
    kal = [l.strip() for l in t.splitlines() if re.search(r'\bFe\b', l) or re.search(r'magnetiz|nspin', l, re.I)]
    if kal: dur('Fe/spin kalintisi: ' + ' | '.join(kal[:3]))
    sp = dict(re.findall(r'(?im)^[ \t]*(Bi|Ga|O)[ \t]+[0-9.]+[ \t]+(\S+)[ \t]*$', t))
    if set(sp) != {'Bi', 'Ga', 'O'}: dur('ATOMIC_SPECIES eksik: %s' % sp)
    ne = sum(zval(sp[e]) * n[e] for e in n)
    if abs(ne - round(ne)) > 1e-6 or int(round(ne)) % 2: dur('%s elektron sayisi %.2f cift degil' % (name, ne))
    return t, int(round(ne)), kay, sp

HEAD = """#!/bin/bash
#SBATCH -A nparan
#SBATCH -p hamsi
#SBATCH -J @J@
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=@T@
#SBATCH -o slurm-%j.out
export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0
"""
JYAK = HEAD.replace('@J@', 'gaB_yak').replace('@T@', '06:00:00') + """R=@B@/01_yakinsama
for n in @N@; do
  [ -s $R/$n/run.in ] || continue
  np=2; case $n in *k1) np=1;; esac
  ( cd $R/$n && { mpirun pw.x -npool $np -in run.in > run.out || mpirun pw.x -npool 1 -in run.in > run.out; } \\
    && grep -q 'JOB DONE' run.out && rm -f out/*.save/wfc* )
done
"""
JREL = HEAD.replace('@J@', 'gaB_S40').replace('@T@', '08:00:00') + """cd @B@/02_O40/gaS_rel || exit 1
mpirun pw.x -npool 2 -in run.in > run.out || { grep -q '^!' run.out || mpirun pw.x -npool 1 -in run.in > run.out; }
"""

def yaz(d, t):
    os.makedirs(d, exist_ok=True)
    if os.path.exists(d + '/run.out'): dur('run.out zaten var, uzerine yazilmaz: ' + d)
    open(d + '/run.in', 'w').write(t)

def uret():
    T = open(TPL).read()
    if not re.search(r"calculation\s*=\s*'relax'", T, re.I): dur("sablon 'relax' degil")
    print('%-10s %-5s %4s %5s %3s %4s %5s' % ('ad', 'calc', 'ecut', 'rho', 'k', 'nat', 'nelec'))
    for y, e, k in YAK:
        t, ne, kay, sp = uret1(T, ad(y, e, k), y, 'scf', e, k)
        yaz('%s/01_yakinsama/%s' % (B, ad(y, e, k)), t)
        print('%-10s %-5s %4d %5d %3d %4d %5d' % (ad(y, e, k), 'scf', e, 4 * e, k, 66 if y == 'S' else 65, ne))
    t, ne, kay, sp = uret1(T, 'gaS_rel', 'S', 'relax', 80, 2)
    yaz(B + '/02_O40/gaS_rel', t)
    print('%-10s %-5s %4d %5d %3d %4d %5d' % ('gaS_rel', 'relax', 80, 320, 2, 66, ne))
    open(B + '/job_yak.sh', 'w').write(JYAK.replace('@B@', B).replace('@N@', ' '.join(ad(*x) for x in YAK)))
    open(B + '/job_gaS.sh', 'w').write(JREL.replace('@B@', B))
    print('k kaymasi (sablondan): %s | turler: %s' % (kay, ' '.join(sp[e] for e in ('Bi', 'Ga', 'O'))))
    print('yazildi: 01_yakinsama/ (%d scf), 02_O40/gaS_rel, job_yak.sh, job_gaS.sh' % len(YAK))

def son(p, t):
    m = re.findall(p, t, re.M)
    return float(m[-1]) if m else float('nan')

def oku1(f):
    if not os.path.exists(f): return None
    t = open(f, errors='ignore').read()
    return dict(t=t, jd=int('JOB DONE' in t), E=son(r'^!\s+total energy\s+=\s+(-?[\d.]+)', t),
                P=son(r'P=\s*(-?[\d.]+)', t), F=son(r'Total force =\s+([\d.]+)', t),
                EF=son(r'the Fermi energy is\s+(-?[\d.]+)', t), nit=son(r'convergence has been achieved in\s+(\d+)', t),
                ns=len(re.findall(r'^!', t, re.M)), bf=int('bfgs converged' in t),
                w=(re.findall(r'PWSCF\s*:.*CPU\s+(.*?)\s*WALL', t) or ['?'])[-1].replace(' ', ''))

def dmin(X, i, a):
    r = []
    for s, x in X:
        if s == 'O':
            d = [x[c] - X[i][1][c] for c in range(3)]
            r.append(a * sum((q - round(q)) ** 2 for q in d) ** .5)
    return ' '.join('%.3f' % q for q in sorted(r)[:4])

def oku():
    R = {ad(*x): oku1('%s/01_yakinsama/%s/run.out' % (B, ad(*x))) for x in YAK}
    print('%-10s %2s %16s %9s %8s %9s %8s %4s %9s' % ('ad', 'JD', 'E(Ry)', 'dE/at_meV', 'P(kbar)', 'F(Ry/b)', 'EF(eV)', 'nit', 'wall'))
    for x in YAK:
        n, r = ad(*x), R[ad(*x)]
        if r is None: print('%-10s YOK' % n); continue
        ref = R.get('S_e080k3' if n in ('S_e080k1', 'S_e080k3') else 'S_e100k2') if x[0] == 'S' else None
        de = (r['E'] - ref['E']) * RY * 1000 / 66 if ref else float('nan')
        print('%-10s %2d %16.8f %9.3f %8.2f %9.5f %8.4f %4.0f %9s' % (n, r['jd'], r['E'], de, r['P'], r['F'], r['EF'], r['nit'], r['w']))
    if R.get('S_e080k2') and R.get('S_e080k3'):
        print('k: dE/at(k2-k3) = %.3f meV (S_e080k1/k3 satirlari k3 referansli; ecut satirlari e100 referansli)'
              % ((R['S_e080k2']['E'] - R['S_e080k3']['E']) * RY * 1000 / 66))
    for e in (80, 100):
        a, b = R.get(ad('A', e, 2)), R.get(ad('B', e, 2))
        if a and b and a['jd'] and b['jd']: print('dE(A-B) e%d = %+.4f eV  (Fe Makale A, relakse: +0.636; burada Fe geometrisi, RELAKSE DEGIL)' % (e, (a['E'] - b['E']) * RY))
        else: print('dE(A-B) e%d : eksik/JD yok' % e)
    r = oku1(B + '/02_O40/gaS_rel/run.out')
    if r is None: print('gaS_rel YOK'); return
    print('gaS_rel: JD %d bfgs %d adim %d E %.8f Ry P %.2f kbar F %.5f wall %s' % (r['jd'], r['bf'], r['ns'], r['E'], r['P'], r['F'], r['w']))
    t = r['t']; a = son(r'celldm\(1\)=\s*([\d.]+)', t) * 0.529177; k = t.rfind('ATOMIC_POSITIONS'); X = []
    for l in t[k:].splitlines()[1:]:
        w = l.split()
        if len(w) < 4 or w[0] not in ('Bi', 'Ga', 'O'): break
        X.append((w[0], [float(q) for q in w[1:4]]))
    if len(X) != 66 or not r['bf']: print('### son konumlar okunamadi ya da bfgs yakinsamadi (%d atom) -> .pos YAZILMADI' % len(X)); return
    open(B + '/02_O40/gaS_rel.pos', 'w').write(''.join('%-3s %14.10f %14.10f %14.10f\n' % (s, *x) for s, x in X))
    ga = [i for i, (s, x) in enumerate(X) if s == 'Ga'][0]
    bt = min((i for i, (s, x) in enumerate(X) if s == 'Bi'), key=lambda i: sum(((X[i][1][c] - .5) - round(X[i][1][c] - .5)) ** 2 for c in range(3)))
    print('a = %.4f A (kubik varsayim) | Ga-O: %s (Fe-O 4x1.917) | Bi_tet-O: %s (Fe hucresinde 4x2.031)' % (a, dmin(X, ga, a), dmin(X, bt, a)))
    print('yazildi: 02_O40/gaS_rel.pos')

if __name__ == '__main__':
    if len(sys.argv) != 2 or sys.argv[1] not in ('uret', 'oku'): dur('kullanim: ga_f1.py uret | oku')
    uret() if sys.argv[1] == 'uret' else oku()
