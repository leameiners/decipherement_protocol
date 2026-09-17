"""Regression check against the Proto-Elamite Dossier's published numbers,
covering P1, P2, P3 (should report not-applicable), P5, P6, P7, P9 and the
P10 totaling-tablet helper -- the fullest exercise of the toolkit so far.
"""
import sys, pickle
sys.path.insert(0, '/tmp/claude-0/-home-user-claude-tests/dd3b3458-a783-5e4a-9d0b-609f2f7ec756/scratchpad/protocol_toolkit')
from decipherment_protocol import Document, Corpus, tests

SCRATCH = '/tmp/claude-0/-home-user-claude-tests/dd3b3458-a783-5e4a-9d0b-609f2f7ec756/scratchpad/'
with open(SCRATCH + 'pe_parsed.pkl', 'rb') as f:
    D = pickle.load(f)
cat, parsed = D['cat'], D['parsed']

docs = []
for pnum, lines in parsed.items():
    segments = [l['signs'] for l in lines]  # keep empty-sign lines too, for segment_numerals alignment
    seg_nums = [l['nums'] for l in lines]
    meta = cat[pnum]
    docs.append(Document(doc_id=pnum, segments=[s for s in segments if s],
                          site=meta['provenience'], artifact_type=meta['object_type'],
                          segment_numerals=seg_nums))

corpus = Corpus(name='Proto-Elamite', documents=docs)
print(f"documents: {len(docs)} (dossier: 1,558)")

p1 = tests.p1_census(corpus)
print(f"\nP1: occurrences={p1['total_occurrences']} (dossier: 14,711)  "
      f"inventory={p1['inventory_size']} (dossier: 608)  hapax={p1['hapax_legomena']} (dossier: 205)")

p2 = tests.p2_rank_frequency(corpus, cutoffs=(30, 60, 100, None))
print("P2 zipf (dossier: s=0.72/0.80/0.81/1.60):")
for cutoff, r in p2['zipf_by_cutoff'].items():
    print(f"  cutoff={cutoff}: s={r['exponent']:.3f} r2={r['r2']:.3f}")
print(f"MPL exponent = -{p2['mpl_exponent']:.3f} (dossier: -1.099)  r2={p2['mpl_r2']:.3f}")

p3 = tests.p3_site_specialization(corpus)
print(f"\nP3 applicable: {p3['applicable']} (dossier: False, 94.7% one site) -- reason: {p3.get('reason')}")

p5 = tests.p5_positional_classes(corpus)
print(f"\nP5 class counts (dossier: initial 59 / final 39 / free 123): {p5['class_counts']}")
m288 = next((r for r in p5['rows'] if r['sign'] == 'M288'), None)
print(f"M288 check: {m288['class'] if m288 else 'not found'} "
      f"frac_final={m288['frac_final']:.2f} n={m288['n']} (dossier: final-preferring, 94%, n=462)")

p6 = tests.p6_cooccurrence_network(corpus)
print(f"\nP6: {len(p6['communities'])} communities, Q={p6['modularity_q']:.4f} (dossier: 5 communities, Q=0.1433)")

p7 = tests.p7_segment_length(corpus)
print(f"\nP7: mean={p7['mean_length']:.3f} median={p7['median_length']} n={p7['n_segments']} "
      f"(dossier: mean=1.784, median=1, n=8,247)")

p9 = tests.p9_periodicity(corpus, p5['classes'], lags=(1, 2, 3), n_perm=500)
print(f"\nP9 (n_perm=500, dossier used 2000 -- expect close but not identical z-scores):")
for lag, r in p9['lag_results'].items():
    print(f"  lag={lag}: z={r['z']:.2f} (dossier: lag1=+18.94, lag2=-0.04, lag3=-9.55)")

p10 = tests.p10_totaling_tablet_test(corpus)
print(f"\nP10 totaling-tablet: tested={p10['tested']} (dossier: 926) hit_rate={p10['hit_rate']*100:.1f}% "
      f"(dossier: 17.7%) null_rate={p10['null_rate']*100:.1f}% (dossier: 13.1%)")
