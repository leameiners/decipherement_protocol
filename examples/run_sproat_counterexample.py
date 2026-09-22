"""An explicit implementation of Sproat's (2010) synthetic counterexample to
Rao et al. (2009), run through the P1/P2/P7/P11 protocol instead of only
being cited.

Sproat (2010) argued that Rao et al.'s (2009) entropy result did not, on its
own, distinguish Indus from a non-linguistic system, and backed the claim
with a constructed counterexample: an artificial sign system with roughly
Indus's inventory size (~400 signs), a Zipf-Mandelbrot frequency
distribution with an exponent in Indus's own range (~1.5), and sequences
built with no dependency beyond what that shared marginal frequency
distribution already implies -- i.e. conditionally independent draws, not a
Markov chain with any real memory. Rao et al.'s (2010) reply argued that
Sproat's construction only matches Indus's low-order (0/1) entropy, and
that extending the same method to block entropies at orders up to 6 shows
real scripts' entropy keeps SCALING DOWN with more context in a way a
conditionally-independent construction cannot.

This module builds that counterexample directly rather than taking either
side's characterization of it on faith, so this project's own P11 test
(decipherment_protocol.tests.p11_conditional_entropy, now run to order 6 by
default -- see its docstring) can be pointed at it exactly as it is pointed
at TALLYGRAM and the seven real corpora. TALLYGRAM (run_synthetic_control.py)
is a positive non-linguistic control: real combinatorial/bookkeeping
structure, no language. This is the complementary NEGATIVE control: no
structure of any kind beyond a shared marginal frequency table, order by
order. If P11's real-vs-null gap (gap_by_order) collapses to ~0 at every
order here while it does not for the real corpora, that is evidence the
extended test has the specificity Rao et al.'s reply claims; if this
construction also shows a nontrivial gap at some order, that is evidence
Sproat's original challenge survives the order-6 extension too, and is
reported as found either way.

Construction, made explicit (Sproat's own paper does not publish exact
generation code, so this is this project's own reasonable reading of the
description above, not a verbatim reproduction):
  - 400 signs (Sproat's own stated inventory size, close to Indus's ~400-700
    depending on allograph-consolidation convention -- see this project's P4
    discussion).
  - Zipf-Mandelbrot marginal with exponent 1.5 (the exponent Sproat reports
    as matching Indus).
  - Each sign in each segment drawn i.i.d. from that fixed marginal --
    NO conditioning on position, on the preceding sign, or on anything else.
    This is deliberately the simplest construction consistent with Sproat's
    description ("conditional independence for bigrams" implies no bigram
    dependency at all, which i.i.d. draws trivially satisfy) and is exactly
    what stats.iid_resample_entropy already implements as this package's own
    null -- so this corpus is, in effect, that null made into a first-class
    citable object in its own right, at the same scale as this project's
    real corpora, rather than only an internal control inside another test.
  - Segment lengths drawn from a geometric distribution with mean 5,
    matching the commonly reported ~5-sign average length of an Indus
    inscription (Wells 2015; see this project's own P7 results on Indus,
    which land in the same range) -- a documented approximation, not a
    resampling of this project's actual Indus segment-length histogram.
  - 4,500 segments, matching the order of magnitude of this project's own
    Indus corpus (see run_indus.py).
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from decipherment_protocol import Document, Corpus, tests

random.seed(2010)

N_SIGNS = 400
ZIPF_EXPONENT = 1.5
N_SEGMENTS = 4500
MEAN_SEGMENT_LENGTH = 5.0

SIGN_VOCAB = [f"SP{i:03d}" for i in range(1, N_SIGNS + 1)]
WEIGHTS = [1.0 / (r ** ZIPF_EXPONENT) for r in range(1, N_SIGNS + 1)]


def geometric_length(mean):
    p = 1.0 / mean
    length = 1
    while random.random() > p:
        length += 1
    return length


docs = []
for doc_i in range(N_SEGMENTS):
    length = geometric_length(MEAN_SEGMENT_LENGTH)
    seg = random.choices(SIGN_VOCAB, weights=WEIGHTS, k=length)
    docs.append(Document(doc_id=f"SP{doc_i:05d}", segments=[seg]))

corpus = Corpus(name="Sproat (2010) counterexample (synthetic, i.i.d. Zipf(1.5), 400 signs)", documents=docs)
total_occ = sum(len(d.flat_signs) for d in docs)
print(f"synthetic documents: {len(docs)}  total occurrences: {total_occ}")

p1 = tests.p1_census(corpus)
print(f"\nP1: occurrences={p1['total_occurrences']} inventory={p1['inventory_size']} "
      f"hapax={p1['hapax_legomena']} ({p1['hapax_share']*100:.1f}%)")

p2 = tests.p2_rank_frequency(corpus, cutoffs=(30, 60, 100, None))
print("\nP2 zipf (sanity check -- should recover ~1.5, the exponent this corpus was built with):")
for cutoff, r in p2['zipf_by_cutoff'].items():
    print(f"  cutoff={cutoff}: s={r['exponent']:.3f} r2={r['r2']:.3f}")
print(f"MPL exponent = {p2['mpl_exponent']:.3f}  r2={p2['mpl_r2']:.3f}")

p7 = tests.p7_segment_length(corpus)
print(f"\nP7: mean={p7['mean_length']:.3f} median={p7['median_length']} n={p7['n_segments']} "
      f"(built from a geometric distribution with target mean {MEAN_SEGMENT_LENGTH})")

# ---------- P11: conditional entropy, orders 0-6 (see module docstring and
# decipherment_protocol.tests.p11_conditional_entropy's docstring) ----------
p11 = tests.p11_conditional_entropy(corpus, max_order=6)
print(f"\nP11 entropy (bits), orders 0-6 (Rao et al. 2010's own block-entropy extension):")
for _k in range(7):
    _row = p11['by_order'][_k]
    print(f"  order{_k}: entropy={_row['entropy']:.3f}  n_contexts={_row['n_contexts']}  n_obs={_row['n_observations']}")
print(f"  within-segment shuffle null (order1): {p11['within_segment_shuffle_order1']['entropy']:.3f} "
      f"(weak null -- keeps each segment's own sign multiset, only randomizes order)")
print("  i.i.d. resample from marginal null, orders 0-6 (correct null -- shares this corpus's own sparsity bias;")
print("  since this corpus IS an i.i.d. construction, this null should track the real values closely at every order):")
for _k in range(7):
    _row = p11['iid_resample_by_order'][_k]
    print(f"    order{_k}: entropy={_row['entropy']:.3f}  n_contexts={_row['n_contexts']}  n_obs={_row['n_observations']}")
print("  real vs i.i.d.-null gap by order (this corpus's own negative-control result -- should sit near 0 bits")
print("  at every order, since it was built with no dependency beyond the shared marginal at all):")
for _k in sorted(p11['gap_by_order']):
    print(f"    order{_k}: {p11['gap_by_order'][_k]:+.3f} bits")
