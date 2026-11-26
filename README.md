# Ga-sillenite-vacancy-DFT

Data and scripts for the article

> N. Paran, *Oxygen-vacancy site selectivity in gallium sillenite requires tetrahedral bismuth: a first-principles counterfactual study of Bi₂₅GaO₃₉ and Bi₂₄Ga₂O₃₉* (in preparation).

**Status:** work in progress. The complete oxygen-site set of Bi₂₅GaO₃₉ (stages F5, F5-b) is included. The O(2)-site relaxations of the counterfactual Bi₂₄Ga₂O₃₉ (stage F5-c; inputs `ga2D1_rel`, `ga2D2_rel`) are running and will be added to `calc/03_site.tar.gz`. A versioned release with a Zenodo DOI will be made at submission; until then file contents may change.

Author: Nejdet Paran, Department of Physics, Faculty of Science, Erciyes University, 38039 Kayseri, Türkiye · ORCID [0000-0003-0317-0438](https://orcid.org/0000-0003-0317-0438)

## Contents

| path | content |
|---|---|
| `calc/<stage>.tar.gz` | Quantum ESPRESSO inputs, outputs and `data-file-schema.xml` of every calculation, one archive per stage (see below). Archives larger than 9 MB are split, e.g. `cat 04_eos.tar.gz.part_* > 04_eos.tar.gz` (same for `03_site_f5`). Wavefunctions, charge densities and pseudopotential files are not included. |
| `data/` | Derived data files read by the figure script, the summary logs written by the analysis scripts (`*_oku.log`), the table file `ga_tablolar.txt`, and the relaxed structures (`*.xsf`, Cartesian Å). |
| `scripts/` | Analysis scripts that generate the inputs (`uret`) and read the outputs (`oku`) of each stage on the cluster; `scripts_md5.txt` lists their MD5 sums. `ga_paket.sh` builds the `calc/` archives. `scripts/jobs/` holds the SLURM job scripts (TRUBA, 56 cores per job). |
| `figures/make_figures_colab.py` | Produces Figs. 1–6 from `data/` (Google Colab or local Python). |
| `MANIFEST.md5` | MD5 sums of all files; check with `md5sum -c MANIFEST.md5`. |

### Calculation stages (`calc/`)

| archive | stage | content |
|---|---|---|
| `01_yakinsama` | F1b | cutoff and k-point convergence (stoichiometric and two vacancy cells) |
| `02_O40` | F1 | stoichiometric Bi₂₅GaO₄₀ (Bi⁵⁺-ordered model), relaxed |
| `03_site` | F2, F2a, F2c | vacancy configurations A, B, C of Bi₂₅GaO₃₉; A from an un-inverted start (`gaA0_rel`); counterfactual Bi₂₄Ga₂O₃₉ (2A, 2B, 2C) |
| `03_site_f5` | F5, F5-b | remaining oxygen sites D1, D2, C′ and control B (`f5_Br`), all started from the stoichiometric cell minus one O with ±0.03 Å random displacements; `f5_Cp` stopped at the 16 h limit and was continued as `f5_Cp2` |
| `04_eos` | F4, F4b, F4c, F2c-b, F4-Si | fixed-volume relaxations for the equations of state of Bi₂₅GaO₃₉ (B), Bi₂₅GaO₄₀, Bi₂₄Ga₂O₃₉ (2A) and Bi₁₂SiO₂₀ |
| `05_elek` | F2b-2 | 4 × 4 × 4 non-self-consistent gaps; PAW projected densities of states (summary `pdos_ga{S,B}.dat`) |
| `06_ref` | F2b-1 | α-Ga and β-Ga₂O₃ reference phases (vc-relax + scf); source CIFs |
| `07_hse` | F2b-4 | HSE06 calculations (exact-exchange cutoff 160 Ry; 320 Ry test for Bi₂₅GaO₄₀) |

## Tables and figures → files

| item | data | produced by |
|---|---|---|
| Table 1 (oxygen sites) | `data/Bi25GaO40_S.xsf` | `scripts/ga_f5.py` (`siniflar`) |
| Sec. 3.2 (D1, D2, C′) | `calc/03_site_f5`, `data/ga_f5_oku.log`, `data/ga_f5b_tani.log`, `data/ga_f5b_oku.log` | `scripts/ga_f5.py oku`, `scripts/ga_f5b.py tani`/`oku` |
| Table 2 (neutron model) | `calc/03_site` (`gaB_rel`), `data/ga_f2_oku.log` | `scripts/ga_f2.py oku` |
| Table 3 (site and formation energies) | `data/ga_tablolar.txt` (T1), `data/ga_ref_oku.log`, `data/ga_f5_oku.log` (D1, D2: E_form = E_form(B) + ΔE) | `scripts/ga_tablo.py`, `scripts/ga_ref.py oku`, `scripts/ga_f5.py oku` |
| Table 4 (Bi₂₄Ga₂O₃₉) | `data/ga_f2c_oku.log`, `data/ga_tablolar.txt` (T3); 2D1/2D2 pending | `scripts/ga_f2c.py oku`, `scripts/ga_f5c.py oku` |
| Table 5 (gaps) | `data/ga_elek_oku.log`, `data/ga_hse_oku.log`, `data/ga_tablolar.txt` (T4) | `scripts/ga_elek.py oku`, `scripts/ga_hse.py oku` |
| Table 6 (equations of state) | `data/eos_*.dat`, `data/ga_f4c_oku.log`, `data/ga_f2cb_oku.log`, `data/ga_f4si_oku.log` | `scripts/ga_f4*.py`, `scripts/ga_f2cb.py`, `scripts/ga_f4si.py` |
| Δa (Sec. 3.3, Fig. 6) | `data/ga_f2cb_oku.log`, `data/eos_ga2A.dat` | `scripts/ga_f2cb.py oku` |
| Fig. 1 | `data/*.xsf` | `figures/make_figures_colab.py` |
| Figs. 2, 3 | `data/ga_tablolar.txt` | idem |
| Fig. 4 | `data/eos_ozet.dat`, `eos_gaS.dat`, `eos_si.dat`, `eos_ga2A.dat` | idem |
| Fig. 5 | `data/pdos_gaS.dat`, `pdos_gaB.dat`, `data/ga_tablolar.txt` | idem |
| Fig. 6 | `data/ga_f2cb_oku.log` (Δa); measured lattice parameters from the literature (in the script, with sources) | idem |

Reproduce the figures locally:

```bash
SILB_LOCAL=data SILB_OUT=figures/out python3 figures/make_figures_colab.py
```

## Computational details (summary)

Quantum ESPRESSO 7.2; PBE; SG15 ONCV pseudopotentials v1.0 (`Bi_ONCV_PBE-1.0.upf`, `Ga_ONCV_PBE-1.0.upf` [MD5 72d1ce27…], `O_ONCV_PBE-1.0.upf`, `Si_ONCV_PBE-1.0.upf` [MD5 2de5a83d…]); 80/320 Ry for relaxations and site energies, 90/360 Ry for equations of state; 2 × 2 × 2 Monkhorst–Pack mesh; Gaussian smearing 0.005 Ry; spin-unpolarized; `etot_conv_thr = 1e-5` Ry, `forc_conv_thr = 1e-4` Ry/bohr; ordered P23 cell, a = 10.3625 Å. PAW projections: pslibrary 1.0.0 (`Bi.pbe-dn-kjpaw_psl.1.0.0.UPF`, `Ga.pbe-dnl-kjpaw_psl.1.0.0.UPF` [MD5 b6829151…], `O.pbe-n-kjpaw_psl.1.0.0.UPF`), 60/480 Ry. HSE06: `ecutfock` 160 Ry, q-mesh 1 × 1 × 1, fixed occupations. Full details in the article.

Starting geometries of the vacancy configurations A, B and C are the relaxed iron-sillenite structures of a companion study (N. Paran, submitted to J. Solid State Chem., JSSC-26-2217) with Fe replaced by Ga; those iron structures belong to that study and are not redistributed here. External structures: α-Ga (COD 9008085), β-Ga₂O₃ (COD 2004987), Bi₁₂SiO₂₀ (COD 1533225).

The scripts use Turkish identifiers and messages (`uret` = generate inputs, `oku` = read results, `DUR` = stop). They are written for the cluster directory layout `/arf/scratch/nparan/SilenitB/` (path set near the top of each script; `ga_f1.py` also reads the environment variables `SILB` and `SILF`); edit that line to run them elsewhere.

## License

Scripts: MIT. Data: CC BY 4.0. See `LICENSE`.
