#!/usr/bin/env python3
# ga_f5c_ga.py : Ga environment of the counterfactual O(2) configurations ga2D1_rel, ga2D2_rel (Table 4, 2D rows).
# Same analysis as ga_f2c.py oku (gacevre: Ga displacement along [111] and Ga-O distances) applied to the F5-c outputs.
# Usage (cluster layout or SILB=<root containing 03_site/>): python3 scripts/ga_f5c_ga.py > data/ga_f5c_ga.log
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ga_f1 as g, ga_f2 as g2, ga_f2c as gc
for k in ('D1', 'D2'):
    r = g.oku1('%s/03_site/ga2%s_rel/run.out' % (g.B, k)); X = g2.konum(r['t'])
    print('ga2%s_rel: JD %d bfgs %d adim %d E %.8f Ry P %.2f kbar F %.5f | %s' % (k, r['jd'], r['bf'], r['ns'], r['E'], r['P'], r['F'], r['w']))
    for s, par, b, O in gc.gacevre(X):
        print('   Ga (%s): ||[111] %+.3f | Ga-O<2.3: %d bag, ort %.3f | en yakin O %s' % (s, par, len(b), sum(b) / max(len(b), 1), ' '.join('%.3f' % o for o in O)))
