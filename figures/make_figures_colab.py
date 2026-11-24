# ===== SilenitB (Ga sillenite) – Figures 1–6 | data: sekil_veri_SilenitB.tar.gz (ga_sekil_veri.py, md5 9dee52ca) =====
# All plotted numbers are read from files in the package; literature constants are marked [lit] with their source.
import os, re, glob, tarfile, numpy as np, matplotlib, matplotlib.pyplot as plt
LOCAL = os.environ.get('SILB_LOCAL')                      # local test only; in Colab leave unset
KLASOR = '/content/drive/MyDrive/Research/DFT/Bi25FeO39/SillenitB'   # Drive: Research/DFT/Bi25FeO39/SillenitB
PAKET = 'sekil_veri_SilenitB.tar.gz'
if LOCAL:
    D, OUT = LOCAL, os.environ.get('SILB_OUT', '.')
else:
    from google.colab import drive
    if not os.path.isdir('/content/drive/MyDrive'): drive.mount('/content/drive')
    tgz = os.path.join(KLASOR, PAKET)
    if not os.path.exists(tgz):
        print('BULUNAMADI:', tgz)
        for p in ('/content/drive/MyDrive', '/content/drive/MyDrive/Research', '/content/drive/MyDrive/Research/DFT/Bi25FeO39'):
            print('  ', p, '->', sorted(os.listdir(p))[:15] if os.path.isdir(p) else 'YOK (farkli Google hesabi ile mi baglandi?)')
        alt = glob.glob('/content/drive/MyDrive/Research/DFT/Bi25FeO39/Sil*/' + PAKET) + glob.glob('/content/drive/*/*/Research/DFT/Bi25FeO39/Sil*/' + PAKET)
        if alt: tgz = alt[0]; print('bulundu:', tgz)
        else:
            from google.colab import files
            print('Paketi elle yukleyin:'); up = files.upload(); assert PAKET in up, 'yuklenen dosya adi ' + PAKET + ' olmali'; tgz = '/content/' + PAKET
    with tarfile.open(tgz) as tf:
        try: tf.extractall('/content', filter='data')
        except TypeError: tf.extractall('/content')
    D = '/content/sekil_veri'; OUT = os.path.join(os.path.dirname(tgz) if tgz.startswith('/content/drive') else '/content', 'sekiller'); os.makedirs(OUT, exist_ok=True)
print('data:', D, '| output:', OUT)

# ---------- style: thin marks, recessive axes; categorical slots 1–3 (validated all-pairs) + neutral gray ----------
C1, C2, C3, GR, INK, INK2 = '#2a78d6', '#eb6834', '#1baf7a', '#8a8984', '#0b0b0b', '#52514e'
matplotlib.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.linewidth': 0.6, 'axes.edgecolor': INK2,
    'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'lines.linewidth': 1.4, 'legend.frameon': False, 'legend.fontsize': 7,
    'savefig.dpi': 600, 'savefig.bbox': 'tight', 'mathtext.default': 'regular'})
W1, W2 = 3.35, 7.0                                          # single / double column (in)
def save(fig, name):
    for ext in ('pdf', 'png'): fig.savefig(f'{OUT}/{name}.{ext}')
    print('saved', name)
    if not LOCAL:                                           # Colab: show each figure under the cell
        from IPython.display import display; display(fig)
    plt.close(fig)

# ---------- read numbers ----------
T = open(f'{D}/ga_tablolar.txt').read()
mu = re.search(r'aralik\s+([\d.]+)\s*eV', T); DMU = float(mu.group(1))
T1 = {m[0]: dict(dE=float(m[1]), rich=float(m[2]), poor=float(m[3])) for m in re.findall(r'^(B|A|C)\s+([+-][\d.]+)\s+([+-][\d.]+)\s+([+-][\d.]+)', T, re.M)}
T3 = {m[0]: float(m[1]) for m in re.findall(r'^ga2(A|B|C)\s+dE\s+([+-][\d.]+)\s*eV', T, re.M)}
gS = re.search(r'S \(Bi25GaO40\): PBE SG15 k4 ([\d.]+).*?HSE06 k2: gaS_hse ([\d.]+)', T); gB = re.search(r'B \(Bi25GaO39\): PBE SG15 k4 ([\d.]+).*?HSE06 k2: gaB_hse ([\d.]+)', T)
GAP = {'S': (float(gS[1]), float(gS[2])), 'B': (float(gB[1]), float(gB[2]))}
L2 = open(f'{D}/ga_f2cb_oku.log').read(); DA = float(re.search(r'\(PBE\),\s*([+-][\d.]+)\s*A olcekli', L2).group(1))
assert set(T1) == {'A', 'B', 'C'} and set(T3) == {'A', 'B', 'C'}, 'ga_tablolar.txt parse'
print('T1', T1, '\nT3', T3, '\nmu range', DMU, '| gaps', GAP, '| Da', DA)

def dat(f):
    R = [l.split() for l in open(f'{D}/{f}') if l.strip() and not l.startswith('#')]
    return np.array([[float(x) for x in r[:5]] for r in R])
def bm3(V, E):                                              # E = sum c_k t^k, t = (V/Vm)^(-2/3) - 1 (same as ga_f4.bm3)
    Vm = V.mean(); t = (V / Vm) ** (-2 / 3) - 1; c = np.polyfit(t, E, 3)
    r = [q.real for q in np.roots(np.polyder(c)) if abs(q.imag) < 1e-12 and np.polyval(np.polyder(c, 2), q.real) > 0]
    q0 = min(r, key=abs); V0 = Vm * (1 + q0) ** -1.5
    def P(v):                                               # GPa
        tv = (v / Vm) ** (-2 / 3) - 1; return -np.polyval(np.polyder(c), tv) * (-2 / 3 / Vm) * (v / Vm) ** (-5 / 3) * 160.21766
    h = 1e-3 * V0; B0 = -V0 * (P(V0 + h) - P(V0 - h)) / (2 * h)
    return dict(V0=V0, B0=B0, E0=np.polyval(c, q0), E=lambda v: np.polyval(c, (v / Vm) ** (-2 / 3) - 1), P=P)

# ---------- Fig. 2: vacancy-site energies, counterfactual ----------
try:
    fig, ax = plt.subplots(figsize=(W1, 2.3))
    cats = ['body-centre\ntetrahedron', 'origin\ntetrahedron', 'lattice O(1)']
    v25 = [T1['B']['dE'], T1['A']['dE'], T1['C']['dE']]       # Bi25GaO39: B = Bi_tet (centre), A = Ga (origin), C = O(1)
    v24 = [T3['B'], T3['A'], T3['C']]                          # Bi24Ga2O39: 2B = Ga (centre), 2A = Ga (origin), 2C = O(1)
    x = np.arange(3); w = 0.36
    for xs, v, col, lab in ((x - w / 2 - 0.01, v25, C1, r'Bi$_{25}$GaO$_{39}$ (Bi$_{tet}$ at centre)'), (x + w / 2 + 0.01, v24, C2, r'Bi$_{24}$Ga$_2$O$_{39}$ (Ga at centre)')):
        ax.bar(xs, v, w, color=col, label=lab, zorder=2)
        for xi, vi in zip(xs, v): ax.text(xi, vi + 0.02, '0' if vi == 0 else (f'{1000 * vi:.1f} meV' if vi < 0.1 else f'{vi:.2f} eV'), ha='center', va='bottom', fontsize=6.5, color=INK)
    ax.set_xticks(x, cats); ax.tick_params(axis='x', length=0); ax.set_ylabel('E – E$_{min}$ (eV)'); ax.set_ylim(0, 1.0); ax.legend(loc='upper left')
    ax.yaxis.grid(True, lw=0.4, color='#e4e3df', zorder=0); ax.spines[['top', 'right']].set_visible(False)
    save(fig, 'Fig2_site_energies')
except Exception:
    import traceback; print("### HATA: Fig. 2"); traceback.print_exc()

# ---------- Fig. 3: formation energy vs Delta mu_O ----------
try:
    fig, ax = plt.subplots(figsize=(W1, 2.4)); dm = np.array([-DMU, 0.])
    for k, col, nm, dy, va in (('B', C1, r'V$_O$ at Bi$_{tet}$ (B)', 0.08, 'bottom'), ('A', C2, r'V$_O$ at Ga (A)', -0.10, 'top'), ('C', C3, 'V$_O$ at O(1) (C)', 0.08, 'bottom')):
        ef = T1[k]['rich'] + dm; ax.plot(dm, ef, color=col)
        ax.text(-DMU + 0.05, ef[0] + dy, nm, color=INK, fontsize=7, ha='left', va=va, rotation=0)
    ax.axhline(0, color=GR, lw=0.6, ls='--')
    ax.set_xlim(-DMU, 0); ax.set_xlabel(r'$\Delta\mu_O$ (eV)'); ax.set_ylabel(r'E$_{form}$ (eV)')
    ax.text(-DMU + 0.03, -2.5, 'Bi + Ga metal limit', fontsize=6.5, color=INK2, ha='left', va='bottom'); ax.text(-0.03, -2.5, r'O-rich ($\frac{1}{2}$O$_2$)', fontsize=6.5, color=INK2, ha='right', va='bottom')
    ax.set_ylim(-2.6, 1.1); ax.spines[['top', 'right']].set_visible(False)
    save(fig, 'Fig3_formation_energy')
except Exception:
    import traceback; print("### HATA: Fig. 3"); traceback.print_exc()

# ---------- Fig. 4: P(V) for four sillenites + Ga/Bi_tet [111] displacement ----------
try:
    oz = [l.split() for l in open(f'{D}/eos_ozet.dat') if l[0] == 'B']
    Bv, BE, BP, Bga, Bbt = (np.array([float(r[i]) for r in oz]) for i in (2, 3, 4, 5, 6)); Bs = np.array([float(r[1]) for r in oz])
    cen = Bs <= 1.0                                             # centred branch (F4c: s >= 1.005 off-centre, excluded)
    fB = bm3(Bv[cen], BE[cen]); S = dat('eos_gaS.dat'); fS = bm3(S[:, 1], S[:, 2]); Si = dat('eos_si.dat'); fSi = bm3(Si[:, 1], Si[:, 2])
    A2 = dat('eos_ga2A.dat'); cq = np.polyfit(A2[:, 1] - A2[:, 1].mean(), A2[:, 3], 2)
    print('check B0: Bi25GaO39 %.1f  Bi25GaO40 %.1f  Bi12SiO20 %.1f GPa' % (fB['B0'], fS['B0'], fSi['B0']))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(W2, 2.6), gridspec_kw=dict(width_ratios=[1.35, 1], wspace=0.3))
    for V, P, f, col, mk, nm in ((Si[:, 1], Si[:, 3], fSi, GR, 'D', r'Bi$_{12}$SiO$_{20}$'), (S[:, 1], S[:, 3], fS, C3, '^', r'Bi$_{25}$GaO$_{40}$'), (Bv[cen], BP[cen], fB, C1, 'o', r'Bi$_{25}$GaO$_{39}$')):
        vv = np.linspace(V.min() - 5, V.max() + 5, 200); a1.plot(vv, f['P'](vv), color=col, lw=1.2)
        a1.plot(V, P, mk, color=col, ms=4.5, mec='white', mew=0.6, label=f'{nm}  (B$_0$ = {f["B0"]:.1f} GPa)')
    a1.plot(Bv[~cen], BP[~cen], 'o', mfc='white', mec=C1, ms=4.5, mew=1.0, label=r'Bi$_{25}$GaO$_{39}$, Ga off-centre (not fitted)')
    vv = np.linspace(A2[:, 1].min() - 5, A2[:, 1].max() + 5, 100); a1.plot(vv, np.polyval(cq, vv - A2[:, 1].mean()), color=C2, lw=1.2, ls='--')
    a1.plot(A2[:, 1], A2[:, 3], 's', color=C2, ms=4.5, mec='white', mew=0.6, label=r'Bi$_{24}$Ga$_2$O$_{39}$ (2A)')
    a1.axhline(0, color=GR, lw=0.6); a1.set_xlabel('V (Å$^3$)'); a1.set_ylabel('P (GPa)'); a1.legend(loc='upper right'); a1.text(-0.13, 1.03, '(a)', transform=a1.transAxes, va='bottom', fontweight='bold')
    o = np.argsort(Bv); a2.plot(Bv[o], Bga[o], '-o', color=C1, ms=4, mec='white', mew=0.6, label='Ga'); a2.plot(Bv[o], Bbt[o], '-s', color=C2, ms=4, mec='white', mew=0.6, label=r'Bi$_{tet}$')
    a2.axvline(fB['V0'], color=GR, lw=0.6, ls='--'); a2.text(fB['V0'] + 2, 0.42, 'V$_0$', color=INK2, fontsize=7)
    a2.axhline(0, color=GR, lw=0.6); a2.set_xlabel('V (Å$^3$)'); a2.set_ylabel('displacement along [111] (Å)'); a2.legend(loc='center left'); a2.text(-0.17, 1.03, '(b)', transform=a2.transAxes, va='bottom', fontweight='bold')
    for a in (a1, a2): a.spines[['top', 'right']].set_visible(False)
    save(fig, 'Fig4_EOS')
except Exception:
    import traceback; print("### HATA: Fig. 4"); traceback.print_exc()

# ---------- Fig. 5: PDOS, S (Bi5+ ordered, O40) vs B (neutron model, O39) ----------
try:
    fig, axs = plt.subplots(2, 1, figsize=(W1, 3.6), sharex=True, gridspec_kw=dict(hspace=0.08))
    for ax, k, ttl in ((axs[0], 'S', r'S: Bi$_{25}$GaO$_{40}$'), (axs[1], 'B', r'B: Bi$_{25}$GaO$_{39}$')):
        h = open(f'{D}/pdos_ga{k}.dat').readline().split()[1:]; d = np.loadtxt(f'{D}/pdos_ga{k}.dat'); col = {n: d[:, i] for i, n in enumerate(h)}
        E = col['E-VBM(eV)']; tot = sum(v for n, v in col.items() if n != 'E-VBM(eV)')
        ax.fill_between(E, tot, color='#e4e3df', lw=0, label='total')
        for y, c, nm in ((col['O_p'], C1, 'O 2p'), (col['Bi_s'] + col['Bi_p'], C2, 'Bi 6s+6p'), (col['Bi_tet_s'] + col['Bi_tet_p'], C3, r'Bi$_{tet}$ 6s+6p')):
            ax.plot(E, y, color=c, lw=1.0, label=nm)
        pbe, hse = GAP[k]; ax.axvspan(0, pbe, color='#f3f2ee', zorder=0); ax.axvline(0, color=GR, lw=0.6, ls='--')
        ax.text(0.98, 0.93, f'{ttl}\nPBE {pbe:.2f} eV | HSE06 {hse:.2f} eV', transform=ax.transAxes, ha='right', va='top', fontsize=7, color=INK)
        sel = (E > -6) & (E < 5); ax.set_ylim(0, tot[sel].max() * 1.15); ax.set_ylabel('PDOS (states/eV)'); ax.spines[['top', 'right']].set_visible(False)
    axs[0].legend(loc='upper left', fontsize=6.5); axs[1].set_xlim(-6, 5); axs[1].set_xlabel('E – E$_{VBM}$ (eV)')
    save(fig, 'Fig5_PDOS')
except Exception:
    import traceback; print("### HATA: Fig. 5"); traceback.print_exc()

# ---------- Fig. 6: experimental lattice parameters vs DFT expectation ----------
try:
    A_R, A_46, A_42 = 10.183, 10.1689, 10.1760    # [lit] Radaev 1989 (neutron, Bi25GaO39); PDF 00-046-0414 (EDS Bi25Ga); PDF 00-042-0180 (nominal Bi24Ga2O39; cell from Levin & Roth 1964)
    fig, ax = plt.subplots(figsize=(W1, 2.0))
    for y, a0, nm in ((1, A_R, 'Radaev 1989\n(single crystal)'), (0, A_46, 'PDF 00-046-0414\n(hydrothermal)')):
        ax.annotate('', xy=(a0 + DA, y), xytext=(a0, y), arrowprops=dict(arrowstyle='->', color=INK2, lw=0.8))
        ax.plot(a0, y, 'o', color=C1, ms=6, mec='white', mew=0.6); ax.plot(a0 + DA, y, 'o', mfc='white', mec=C2, ms=6, mew=1.2)
    ax.axvline(A_42, color=GR, lw=0.8, ls='--'); ax.text(A_42 + 0.001, 1.75, 'PDF 00-042-0180\n"Bi$_{24}$Ga$_2$O$_{39}$"', fontsize=6.5, color=INK2, va='top')
    ax.plot([], [], 'o', color=C1, label=r'measured, Bi$_{25}$Ga sillenite'); ax.plot([], [], 'o', mfc='white', mec=C2, label=fr'+ $\Delta a$ = {DA:.3f} Å (DFT) → Bi$_{{24}}$Ga$_2$O$_{{39}}$')
    ax.set_ylim(-1.5, 1.8); ax.set_yticks([0, 1], ['PDF 00-046-0414\n(hydrothermal)', 'Radaev 1989\n(single crystal)'], fontsize=6.5); ax.tick_params(axis='y', length=0)
    ax.set_xlim(10.130, 10.192); ax.set_xlabel('a (Å)'); ax.legend(loc='lower center', fontsize=6.3, ncol=1)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    save(fig, 'Fig6_lattice_parameter')
except Exception:
    import traceback; print("### HATA: Fig. 6"); traceback.print_exc()

# ---------- Fig. 1: structure (XSF files of the package; orthographic ball-and-stick, no external program) ----------
try:
    from matplotlib.patches import Circle, Polygon
    from matplotlib.collections import LineCollection
    def rxsf(f):
        L = open(f'{D}/{f}').read().split('\n'); i = L.index('PRIMVEC'); a = float(L[i + 1].split()[0]); j = L.index('PRIMCOORD'); n = int(L[j + 1].split()[0])
        w = [l.split() for l in L[j + 2:j + 2 + n]]; return a, np.array([x[0] for x in w]), np.array([[float(v) for v in x[1:4]] for x in w])
    mic = lambda d, a: d - a * np.round(d / a)
    COL = {'Bi': '#8e6bb8', 'Ga': '#3a9e5c', 'Si': '#c8a052', 'O': '#d9453b'}; RAD = {'Bi': 0.42, 'Ga': 0.36, 'O': 0.26}
    def rot(az, el):
        a, e = np.radians(az), np.radians(el)
        Rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]]); Rx = np.array([[1, 0, 0], [0, np.cos(e), -np.sin(e)], [0, np.sin(e), np.cos(e)]])
        return Rx @ Rz                                          # screen x = row0, screen y = row2, depth = row1 (larger = farther)
    R = rot(-28, 18)
    def P(x): q = R @ x; return q[0], q[2], q[1]
    def ball(ax, x, el, z0, r=None, alpha=1.0, ls='-', fc=None):
        u, v, d = P(x); r = r or RAD[el]; c = fc or COL[el]
        ax.add_patch(Circle((u, v), r, fc=c, ec='white' if ls == '-' else c, lw=0.5, ls=ls, alpha=alpha, zorder=z0 - d))
        if ls == '-' and fc != 'none': ax.add_patch(Circle((u - 0.3 * r, v + 0.3 * r), 0.35 * r, fc='white', ec='none', alpha=0.35 * alpha, zorder=z0 - d + 0.001))
    def stick(ax, x1, x2, z0, c='#b9b8b3', lw=0.6):
        (u1, v1, d1), (u2, v2, d2) = P(x1), P(x2); ax.plot([u1, u2], [v1, v2], color=c, lw=lw, zorder=z0 - (d1 + d2) / 2 - 0.5, solid_capstyle='round')
    def tet(ax, Xs, fc, z0, alpha=0.35):
        import itertools
        for f in itertools.combinations(range(len(Xs)), 3):
            Q = [P(Xs[k]) for k in f]; ax.add_patch(Polygon([(q[0], q[1]) for q in Q], closed=True, fc=fc, ec=fc, lw=0.6, alpha=alpha, zorder=z0 - np.mean([q[2] for q in Q]) - 0.3))

    aS, eS, XS = rxsf('Bi25GaO40_S.xsf'); sh = np.full(3, aS / 4)          # shift: Ga -> a/4, Bi_tet -> 3a/4 (both tetrahedra inside the drawn cell)
    W = lambda X: (X + sh) % aS
    iGa = int(np.where(eS == 'Ga')[0][0]); iBt = int(np.argmin([np.linalg.norm(mic(x - aS / 2, aS)) if e == 'Bi' else 99 for e, x in zip(eS, XS)]))
    iO = np.where(eS == 'O')[0]
    def nb(X, el, i, n, rmax=2.6):                              # nearest O around atom i, unwrapped next to it
        d = np.array([mic(X[j] - X[i], aS) for j in np.where(el == 'O')[0]]); k = np.argsort(np.linalg.norm(d, axis=1))[:n]
        return [X[i] + d[m] for m in k if np.linalg.norm(d[m]) < rmax]
    # vacancy sites (S frame): A = O of the Ga tetrahedron on +[111]; B = O of the Bi_tet tetrahedron on +[111]; C = O missing in Bi25GaO39_C (12j)
    def missing(f):
        _, e2, X2 = rxsf(f); O2 = X2[e2 == 'O']; dd = [np.linalg.norm(mic(O2 - XS[j], aS), axis=1).min() for j in iO]
        return iO[int(np.argmax(dd))], dd
    site = {'A': missing('Bi25GaO39_A.xsf')[0], 'B': missing('Bi25GaO39_B.xsf')[0]}
    _, ddC = missing('Bi25GaO39_C.xsf'); fr = (XS[iO] / aS) % 1
    site['C'] = iO[int(np.argmin([np.linalg.norm(mic(f - np.array([0.01055, 0.24769, 0.36709]), 1)) for f in fr]))]   # Makale A / F2 C definition
    print('vacancy sites (frac):', {k: (XS[v] / aS % 1).round(4) for k, v in site.items()})

    fig = plt.figure(figsize=(W2, 3.1)); ax = fig.add_axes([0, 0, 0.56, 1]); ax.set_aspect('equal'); ax.axis('off')
    X = W(XS); Z = 100
    for c in [np.array(p) for p in np.ndindex(2, 2, 2)]:          # cell edges
        for k in range(3):
            if c[k] == 0: e2 = c.copy(); e2[k] = 1; stick(ax, c * aS, e2 * aS, Z, c='#cfcec9', lw=0.6)
    for i, e in enumerate(eS):                                     # framework Bi-O sticks (< 2.35 A)
        if e != 'Bi' or i == iBt: continue
        for j in iO:
            d = mic(X[j] - X[i], aS)
            if np.linalg.norm(d) < 2.35 and np.all((X[i] + d >= -0.01) & (X[i] + d <= aS + 0.01)): stick(ax, X[i], X[i] + d, Z)
    tet(ax, nb(X, eS, iGa, 4), COL['Ga'], Z); tet(ax, nb(X, eS, iBt, 4), COL['Bi'], Z)
    for i, e in enumerate(eS):
        if i in site.values(): continue
        ball(ax, X[i], e, Z, r=RAD[e] * (1.15 if i in (iGa, iBt) else 1.0))
    for k, c in (('A', C2), ('B', C1), ('C', C3)):
        x = X[site[k]]; ball(ax, x, 'O', Z, r=0.30, fc=c); u, v, d = P(x)
        ax.add_patch(Circle((u, v), 0.55, fc='none', ec=c, lw=1.0, ls=(0, (2, 1.5)), zorder=Z - d + 0.01))
        ax.text(u + 0.55, v + 0.45, k, color=INK, fontsize=8, fontweight='bold', zorder=Z + 50)
    for i, lab, du, dv in ((iGa, 'Ga', -0.3, -1.0), (iBt, r'Bi$_{tet}$', 0.9, -0.2)):
        u, v, d = P(X[i]); ax.text(u + du, v + dv, lab, fontsize=7, color=INK, ha='center', zorder=Z + 50, bbox=dict(fc='white', ec='none', alpha=0.7, pad=0.5))
    for e, nm in (('Bi', 'Bi'), ('Ga', 'Ga'), ('O', 'O')): ax.plot([], [], 'o', color=COL[e], mec='white', ms=6 if e != 'O' else 4.5, label=nm)
    ax.plot([], [], 'o', mfc='none', mec=INK2, ms=7, ls='none', label='vacancy site (A, B, C)')
    ax.legend(loc='lower left', bbox_to_anchor=(0.0, -0.02), fontsize=6.5, ncol=4, handletextpad=0.3, columnspacing=0.8)
    ax.autoscale_view(); ax.text(0.02, 0.97, '(a)', transform=ax.transAxes, fontweight='bold', va='top')

    def local(axl, f, idx, n, vac, lab, col, shift_lbl):
        a2, e2, X2 = rxsf(f); c0 = (XS[idx] + sh); x0 = c0 + mic(X2[idx] - XS[idx], aS)  # relaxed cation, in shifted frame
        Os = [x0 + mic(o - X2[idx], aS) for o in [X2[j] for j in np.where(e2 == 'O')[0]]]; Os = sorted(Os, key=lambda o: np.linalg.norm(o - x0))[:n]
        vx = c0 + mic(XS[vac] - XS[idx], aS)
        tet(axl, Os + [x0] if n == 3 else Os, col, Z, alpha=0.25)
        for o in Os: stick(axl, x0, o, Z, c='#8a8984', lw=0.8); ball(axl, o, 'O', Z)
        ball(axl, x0, e2[idx], Z, r=RAD[e2[idx]] * 1.15)
        u, v, d = P(vx); axl.add_patch(Circle((u, v), 0.30, fc='none', ec=INK2, lw=0.9, ls=(0, (2, 1.5)), zorder=Z + 40)); axl.text(u + 0.35, v + 0.2, r'V$_O$', fontsize=7, color=INK)
        u0, v0, _ = P(c0); axl.plot(u0, v0, '+', color=INK, ms=5, mew=0.8, zorder=Z + 41)
        dist = sorted(np.linalg.norm(o - x0) for o in Os)
        axl.set_title(lab, fontsize=7.5, color=INK, pad=2); axl.set_aspect('equal'); axl.axis('off'); axl.autoscale_view()
        axl.text(0.5, -0.04, shift_lbl + '\n' + ' + '.join(f'{d:.2f}' for d in dist) + ' Å', transform=axl.transAxes, ha='center', va='top', fontsize=6.3, color=INK2)
        return np.linalg.norm(mic(X2[idx] - XS[idx], aS))
    ab = fig.add_axes([0.60, 0.20, 0.18, 0.65]); aa = fig.add_axes([0.80, 0.20, 0.18, 0.65])
    sB = local(ab, 'Bi25GaO39_B.xsf', iBt, 3, site['B'], r'B: Bi$_{tet}$O$_3$ + V$_O$', COL['Bi'], 'Bi$_{tet}$ shift')
    sA = local(aa, 'Bi25GaO39_A.xsf', iGa, 4, site['A'], r'A: inverted GaO$_4$', COL['Ga'], 'Ga shift')
    ab.texts[-1].set_text(ab.texts[-1].get_text().replace('shift', f'shift {sB:.2f} Å')); aa.texts[-1].set_text(aa.texts[-1].get_text().replace('shift', f'shift {sA:.2f} Å'))
    fig.text(0.60, 0.95, '(b)', fontweight='bold', va='top'); fig.text(0.80, 0.95, '(c)', fontweight='bold', va='top')
    save(fig, 'Fig1_structure')
except Exception:
    import traceback; print("### HATA: Fig. 1"); traceback.print_exc()

print('done ->', OUT)
