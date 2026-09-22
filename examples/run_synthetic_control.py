"""A synthetic, compositionally-structured but non-linguistic control corpus
for the P1-P11 protocol.

Raghavendra (2026) builds SIGIL, a purpose-built generative emblem system
with explicit compositional meaning and no phonological value, and shows it
reproduces the same repetition, directional-asymmetry, and lexical-
distribution regularities the undeciphered-script literature treats as
evidence of encoded language -- meaning those regularities are not, on their
own, specific to language. This module is this project's own answer to the
same question, independently designed rather than a reproduction of SIGIL:
does the full P1-P11 battery this paper runs on seven real scripts also find
"structure" on a system built to have structure but no language in it at
all?

The generator, TALLYGRAM: a small synthetic administrative-tally system with
no linguistic content whatsoever. Each "record" (this corpus's segment,
directly comparable to one line-item on a real Bronze Age tally) is three
fixed-position slots -- AGENT, COMMODITY, QUANTITY -- filled from three
disjoint closed vocabularies (40, 25, and 12 marks respectively, chosen to
sit in the same order of magnitude as this paper's real sign inventories).
None of these marks has a sound value, a meaning beyond its own slot
identity, or any grammar; the only "rules" governing the system are:
  - Slot order is fixed (AGENT always first, COMMODITY always second,
    QUANTITY always third) -- exactly the kind of positional regularity P5
    is built to detect, present here by pure combinatorial construction.
  - AGENT and COMMODITY popularity are Zipf-weighted (rank ~ 1/rank), not
    uniform -- realistic organizational skew (a few common trading partners
    and commodities, many rare ones), not a linguistic property.
  - Each AGENT has a preferred subset of 3-5 commodities drawn 70% of the
    time (30% uniform over the full commodity list) -- real co-occurrence
    structure for P6 to find, again with no linguistic content.
  - A "document" (tablet analog) is 1-10 records; 30% of documents with 3+
    records close with a TOTAL record (TOTAL-mark + the true arithmetic sum
    of that document's quantities) -- exactly the ledger structure P10 is
    built to detect, present here as pure bookkeeping, not language.
  - Each document carries a synthetic DEPARTMENT (site analog, 4 values)
    correlated with which agents appear in it, and an artifact-type analog
    (single-line vs. multi-line record) -- giving P3 a real site x type
    table to test, again by construction, not organic administrative
    practice.

If this system's P1/P2/P5/P6/P7/P9/P11 profile is difficult to distinguish
from this paper's seven real corpora, that is itself the honest result: it
would mean this battery, like the ones Raghavendra surveys, detects
organization generally rather than language specifically. If it is clearly
distinguishable on some tests, that is evidence those specific tests do
carry some real specificity -- also reported as found, not adjusted to fit
either expectation in advance.
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from decipherment_protocol import Document, Corpus, tests

random.seed(2026)

N_AGENTS, N_COMMODITIES, N_QUANTITIES = 40, 25, 12
N_DOCUMENTS = 4000
DEPARTMENTS = ['dept-A', 'dept-B', 'dept-C', 'dept-D']

AGENT_VOCAB = [f"AG{i:02d}" for i in range(1, N_AGENTS + 1)]
COMMODITY_VOCAB = [f"CM{i:02d}" for i in range(1, N_COMMODITIES + 1)]
QUANTITY_VOCAB = [f"Q{i:02d}" for i in range(1, N_QUANTITIES + 1)]  # values 1..12
TOTAL_MARK = 'TOTAL'

def zipf_weights(n):
    return [1.0 / r for r in range(1, n + 1)]

agent_weights = zipf_weights(N_AGENTS)
commodity_weights = zipf_weights(N_COMMODITIES)
quantity_weights = zipf_weights(N_QUANTITIES)  # small quantities more common, like real tallies

# each agent has a preferred subset of 3-5 commodities (drawn 70% of the time)
agent_preferred_commodities = {
    a: random.sample(COMMODITY_VOCAB, random.randint(3, 5)) for a in AGENT_VOCAB
}
# departments correlate with a subset of agents (10 agents "belong" to each department)
shuffled_agents = AGENT_VOCAB[:]
random.shuffle(shuffled_agents)
department_agents = {d: shuffled_agents[i * 10:(i + 1) * 10] for i, d in enumerate(DEPARTMENTS)}
agent_department = {a: d for d, agents in department_agents.items() for a in agents}

docs = []
for doc_i in range(N_DOCUMENTS):
    n_records = random.randint(1, 10)
    dept = random.choice(DEPARTMENTS)
    # 80% of a document's agents come from its own department, 20% from anywhere
    dept_agents = department_agents[dept]
    segments, seg_nums = [], []
    quantities_this_doc = []
    for _ in range(n_records):
        agent = random.choice(dept_agents) if random.random() < 0.8 else random.choices(AGENT_VOCAB, weights=agent_weights, k=1)[0]
        if random.random() < 0.7:
            commodity = random.choice(agent_preferred_commodities[agent])
        else:
            commodity = random.choices(COMMODITY_VOCAB, weights=commodity_weights, k=1)[0]
        q_idx = random.choices(range(1, N_QUANTITIES + 1), weights=quantity_weights, k=1)[0]
        quantity_sign = f"Q{q_idx:02d}"
        segments.append([agent, commodity, quantity_sign])
        seg_nums.append([(float(q_idx), 'qty')])
        quantities_this_doc.append(q_idx)
    closes_with_total = n_records >= 3 and random.random() < 0.3
    if closes_with_total:
        segments.append([TOTAL_MARK])
        seg_nums.append([(float(sum(quantities_this_doc)), 'qty')])
    artifact_type = 'multi-line-ledger' if n_records >= 4 else 'single-line'
    docs.append(Document(doc_id=f"TG{doc_i:05d}", segments=segments, segment_numerals=seg_nums,
                          site=dept, artifact_type=artifact_type))

corpus = Corpus(name='TALLYGRAM (synthetic non-linguistic control)', documents=docs)
total_records = sum(len(d.segments) for d in docs)
print(f"synthetic documents: {len(docs)}  total records (segments): {total_records}")

p1 = tests.p1_census(corpus)
print(f"\nP1: occurrences={p1['total_occurrences']} inventory={p1['inventory_size']} "
      f"hapax={p1['hapax_legomena']} ({p1['hapax_share']*100:.1f}%)")

p2 = tests.p2_rank_frequency(corpus, cutoffs=(30, 60, 100, None))
print("\nP2 zipf:")
for cutoff, r in p2['zipf_by_cutoff'].items():
    print(f"  cutoff={cutoff}: s={r['exponent']:.3f} r2={r['r2']:.3f}")
print(f"MPL exponent = {p2['mpl_exponent']:.3f}  r2={p2['mpl_r2']:.3f}")

p3 = tests.p3_site_specialization(corpus)
print(f"\nP3 applicable: {p3['applicable']}")
if p3['applicable']:
    top = sorted(p3['residuals'].items(), key=lambda kv: -kv[1]['adj_resid'])[:5]
    for (site, atype), r in top:
        print(f"  {site} x {atype}: obs={r['observed']} exp={r['expected']:.1f} adj_resid={r['adj_resid']:.1f}")

p5 = tests.p5_positional_classes(corpus)
print(f"\nP5 class counts: {p5['class_counts']}")
print("  (by construction: AGENT marks should land initial-preferring, QUANTITY marks")
print("  final-preferring in 3-slot records, COMMODITY marks free/medial -- purely from fixed")
print("  slot order, with zero phonological or grammatical content behind any of it)")

p6 = tests.p6_cooccurrence_network(corpus)
print(f"\nP6: {len(p6['communities'])} communities, Q={p6['modularity_q']:.4f}")

p7 = tests.p7_segment_length(corpus)
print(f"\nP7: mean={p7['mean_length']:.3f} median={p7['median_length']} n={p7['n_segments']}")

community_classes = {s: cid for cid, members in p6['communities'].items() for s in members}
# every record is exactly 3 signs (the fixed AGENT/COMMODITY/QUANTITY slots), so lag 3 is
# structurally untestable (no segment is long enough) -- lags capped at 1-2 accordingly.
p9 = tests.p9_periodicity(corpus, community_classes, lags=(1, 2), n_perm=2000)
print(f"\nP9 (P6 community label, n_sequences={p9['n_sequences']}):")
for lag, r in p9['lag_results'].items():
    if r is None:
        print(f"  lag={lag}: not testable (no segment long enough)")
    else:
        print(f"  lag={lag}: z={r['z']:.2f}")

p10 = tests.p10_totaling_tablet_test(corpus)
print(f"\nP10 unrestricted: tested={p10['tested']} hit_rate={p10['hit_rate']*100:.1f}% null_rate={p10['null_rate']*100:.1f}%")
p10_total = tests.p10_totaling_tablet_test(corpus, marker_predicate=lambda signs: TOTAL_MARK in signs)
print(f"P10 restricted to TOTAL-marked records: tested={p10_total['tested']} "
      f"hit_rate={p10_total['hit_rate']*100:.1f}% null_rate={p10_total['null_rate']*100:.1f}%")
print("  (by construction: every TOTAL record's stated value IS the exact arithmetic sum of its")
print("  document's preceding quantities -- this should score at or near 100%, since it is pure")
print("  bookkeeping with no language involved, the same point Damerow & Englund's own method")
print("  makes about what P10 actually measures)")

p11 = tests.p11_conditional_entropy(corpus)
print(f"\nP11 entropy (bits): order0={p11['by_order'][0]['entropy']:.3f} "
      f"order1={p11['by_order'][1]['entropy']:.3f} order2={p11['by_order'][2]['entropy']:.3f}")
print(f"  within-segment shuffle null (order1): {p11['within_segment_shuffle_order1']['entropy']:.3f} "
      f"(weak null -- keeps each segment's own sign multiset, only randomizes order)")
print(f"  i.i.d. resample from marginal null (order1): {p11['iid_resample_order1']['entropy']:.3f} "
      f"(correct null -- shares the real corpus's sparsity bias)")
print(f"  real order1 vs i.i.d. null: {p11['by_order'][1]['entropy'] - p11['iid_resample_order1']['entropy']:+.3f} bits "
      f"(the genuine-structure signal; compare to order0-vs-order1 raw drop of "
      f"{p11['entropy_drop_0_to_1']:.3f}, which is inflated by estimator bias, not real structure alone)")
