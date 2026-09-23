"""Ur III Sumerian as a second deciphered-control corpus for the P1-P10 protocol.

Ur III (ca. 2100-2000 BCE) is the single largest administrative cuneiform
archive from the ancient world -- tens of thousands of economic tablets from
the bureaucracy of the Third Dynasty of Ur, almost entirely accounting
records. Sumerian has been read since the mid-19th century; this is a second,
independent retrospective check alongside Linear B, on a script with a
different structure (logo-syllabic cuneiform, not a syllabary) and, unlike
Linear B, in exactly the genre this protocol's P10 is modelled on -- Damerow
and Englund's own totaling-tablet method was developed ON administrative
cuneiform archives like this one.

Source: github.com/cdli-gh/data (the same official CDLI bulk-data dump this
project's Proto-Elamite adapter already uses), filtered to period == "Ur III
(ca. 2100-2000 BC)" and genre == "Administrative": 106,802 catalogued texts,
of which 74,678 (69.9%) have a transliteration in the ATF release -- the rest
are catalogued but not yet transliterated, a property of the source, reported
as a subset the same way this project's other partial corpora are.

Numeral system: Sumerian accounting used several parallel numeral systems --
a base counting system (disz=1, u=10, gesz2=60, szar2=3600, szargal=36000)
plus separate capacity, area, and weight systems with their own conversion
factors (barig, ban2, sila3 for grain; iku for area; and more). Only the base
counting system is decoded here; texts whose totaling line's value can't be
attributed to it are excluded from P10 rather than guessed at, the same
scoping decision this project made for Linear B's ideogram-coded quantities.
This likely undercounts P10's hit rate here more than for this project's
other corpora, since many Ur III accounts total in the grain/area/weight
systems this adapter doesn't decode.

To regenerate the source data:
    git clone --filter=blob:none https://github.com/cdli-gh/data
    (cdli_cat.csv and cdliatf_unblocked.atf are Git-LFS pointers; fetch the
    real files via https://media.githubusercontent.com/media/cdli-gh/data/master/<file>,
    the same approach this project's Proto-Elamite adapter documents.)
"""
import sys, os, re, csv
csv.field_size_limit(10**7)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from decipherment_protocol import Document, Corpus, tests

CDLI_DIR = os.environ.get('CDLI_DATA_DIR', '/home/user/cdli-data')
CAT = CDLI_DIR + '/cdli_cat_real.csv'
ATF = CDLI_DIR + '/cdliatf_real.atf'

BASE_VALUES = {'disz': 1, 'u': 10, 'gesz2': 60, 'szar2': 3600, 'szargal': 36000}
NUMERAL_RE = re.compile(r'^(\d+)\((' + '|'.join(BASE_VALUES) + r')\)$')

STRUCT_TAG = re.compile(r'^[@$#&]')
LINE_LABEL = re.compile(r"^\d+'?\.\s*")
SCORE_REF = re.compile(r'^>>')


def parse_word(token):
    token = token.strip('#!?*<>[]').strip()
    if not token or token in ('x', '...', 'n', 'N'):
        return None
    signs = []
    for part in re.split(r'[-.]', token):
        part = re.sub(r'\{([^}]*)\}', r'\1', part)  # determinative braces -> keep the sign inside
        part = part.strip('#!?*<>[]').strip()
        if not part or part in ('x', '...'):
            continue
        signs.append(part.lower())
    return signs or None


def parse_numeral(token):
    token = token.strip('#!?*')
    m = NUMERAL_RE.match(token)
    if m:
        return (float(m.group(1)) * BASE_VALUES[m.group(2)], 'count')
    return None


def parse_line(text):
    text = LINE_LABEL.sub('', text.strip())
    signs, nums = [], []
    for tok in text.split():
        if tok in ('$', '@'):
            continue
        n = parse_numeral(tok)
        if n:
            nums.append(n)
            continue
        ps = parse_word(tok)
        if ps:
            signs.extend(ps)
    return signs, nums


cat = {}
with open(CAT, encoding='utf-8', errors='replace') as f:
    for row in csv.DictReader(f):
        if row.get('period') == 'Ur III (ca. 2100-2000 BC)' and row.get('genre') == 'Administrative':
            pid = 'P' + str(row.get('id_text') or row.get('id')).zfill(6)
            cat[pid] = {'provenience': row.get('provenience')}
print(f"Ur III administrative catalog entries: {len(cat)}")

blocks = {}
cur_id, cur_lines = None, []
with open(ATF, encoding='utf-8', errors='replace') as f:
    for line in f:
        if line.startswith('&P'):
            if cur_id and cur_id in cat:
                blocks[cur_id] = cur_lines
            pid = line.split()[0][1:]
            cur_id = pid if pid in cat else None
            cur_lines = []
        elif cur_id:
            cur_lines.append(line.rstrip('\n'))
    if cur_id and cur_id in cat:
        blocks[cur_id] = cur_lines
print(f"transliterated blocks matched: {len(blocks)} ({len(blocks)/len(cat)*100:.1f}%)")

word_docs, line_docs = [], []
total_content_lines = 0
for pid, raw_lines in blocks.items():
    segments, seg_nums = [], []
    for raw in raw_lines:
        raw = raw.rstrip()
        if not raw or STRUCT_TAG.match(raw) or SCORE_REF.match(raw):
            continue
        if not re.match(r"^\d", raw):
            continue  # not a numbered text line (editorial notes, etc.)
        total_content_lines += 1
        signs, nums = parse_line(raw)
        if signs or nums:
            segments.append(signs)
            seg_nums.append(nums)
    if segments:
        meta = cat[pid]
        line_docs.append(Document(doc_id=pid, segments=segments, segment_numerals=seg_nums,
                                   site=meta.get('provenience') or None))
        word_segments = [s for s in segments if s]
        if word_segments:
            word_docs.append(Document(doc_id=pid, segments=word_segments,
                                       site=meta.get('provenience') or None))

print(f"documents with parsed content: {len(line_docs)}")

word_corpus = Corpus(name='Ur III Sumerian (words)', documents=word_docs)

p1 = tests.p1_census(word_corpus)
print(f"\nP1: occurrences={p1['total_occurrences']} inventory={p1['inventory_size']} "
      f"hapax={p1['hapax_legomena']} ({p1['hapax_share']*100:.1f}%)")

p2 = tests.p2_rank_frequency(word_corpus, cutoffs=(30, 60, 100, None))
print("\nP2 zipf:")
for cutoff, r in p2['zipf_by_cutoff'].items():
    print(f"  cutoff={cutoff}: s={r['exponent']:.3f} r2={r['r2']:.3f}")
print(f"MPL exponent = -{p2['mpl_exponent']:.3f}  r2={p2['mpl_r2']:.3f}")

p3 = tests.p3_site_specialization(word_corpus)
print(f"\nP3 applicable: {p3['applicable']}")
if p3['applicable']:
    top = sorted(p3['residuals'].items(), key=lambda kv: -kv[1]['adj_resid'])[:5]
    for site, r in top:
        print(f"  {site}: obs={r['observed']} exp={r['expected']:.1f} adj_resid={r['adj_resid']:.1f}")
else:
    print(f"  reason: {p3.get('reason')} -- this adapter only pulls a single genre "
          f"(Administrative) with no second artifact-type axis recorded, so unlike "
          f"Proto-Elamite's P3 inapplicability (a real 94.7%-one-site concentration), "
          f"this one is a scoping artifact of the adapter, not a corpus property")

p5 = tests.p5_positional_classes(word_corpus)
print(f"\nP5 class counts: {p5['class_counts']}")
known = {
    'sze3': 'terminative/allative case suffix "to/for" -- should be final-preferring',
    'ta': 'ablative case suffix "from" -- should be final-preferring',
    'ka': 'locative-genitive case suffix "in/of" -- should be final-preferring',
    'lugal': 'noun "king" -- should be free/content, not an ending',
}
rows_by_sign = {row['sign']: row for row in p5['rows']}
for sign, note in known.items():
    row = rows_by_sign.get(sign)
    if row:
        print(f"  known '{sign}' ({note}): class={row['class']} frac_final={row['frac_final']:.2f} n={row['n']}")

p6 = tests.p6_cooccurrence_network(word_corpus)
print(f"\nP6: {len(p6['communities'])} communities, Q={p6['modularity_q']:.4f}")

p7 = tests.p7_segment_length(word_corpus)
print(f"\nP7: mean={p7['mean_length']:.3f} median={p7['median_length']} n={p7['n_segments']}")

community_classes = {s: cid for cid, members in p6['communities'].items() for s in members}
# P9's permutation test is O(n_perm x total label positions); at ~1M word
# segments across this corpus, the full n_perm=2000 used for the other
# scripts is impractically slow in pure Python. A random 8,000-document
# subsample (still an order of magnitude larger than any of this package's
# other corpora) with n_perm=200 keeps the same test tractable while
# remaining a fair, representative sample -- documented here, not silently
# downsampled.
import random
random.seed(7)
p9_sample = random.sample(word_docs, min(8000, len(word_docs)))
p9_corpus = Corpus(name='Ur III Sumerian (P9 subsample)', documents=p9_sample)
print(f"\nP9 subsample: {len(p9_sample)} / {len(word_docs)} documents")
p9 = tests.p9_periodicity(p9_corpus, community_classes, lags=(1, 2, 3), n_perm=200)
print(f"P9 (P6 community label, n_sequences={p9['n_sequences']}):")
for lag, r in p9['lag_results'].items():
    print(f"  lag={lag}: z={r['z']:.2f}")

# ---------- P10 ----------
line_corpus = Corpus(name='Ur III Sumerian (lines)', documents=line_docs)

def has_sunigin(signs):
    return any(signs[i] == 'szu' and signs[i + 1].startswith('nigin') for i in range(len(signs) - 1))

p10 = tests.p10_totaling_tablet_test(line_corpus)
print(f"\nP10 unrestricted: tested={p10['tested']} hit_rate={p10['hit_rate']*100:.1f}% null_rate={p10['null_rate']*100:.1f}%")

p10_sunigin = tests.p10_totaling_tablet_test(line_corpus, marker_predicate=has_sunigin)
print(f"P10 restricted to szu-nigin ('grand total'): tested={p10_sunigin['tested']} "
      f"hit_rate={p10_sunigin['hit_rate']*100:.1f}% null_rate={p10_sunigin['null_rate']*100:.1f}%")
print("  NOTE: both rates beat their null by a real margin, but the absolute hit rate is")
print("  well under Linear A's 36.4% ku-ro-restricted figure -- likely because this adapter")
print("  decodes only Sumerian's base counting system (disz/u/gesz2/szar2/szargal), not the")
print("  separate grain, area, and weight systems many Ur III accounts total in, so many real")
print("  matches whose total is stated in those units are silently excluded rather than tested.")

# ---------- P11: conditional entropy (Rao et al. 2009-style; see README's
# "Raghavendra (2026) and the entropy test" section) ----------
p11 = tests.p11_conditional_entropy(word_corpus, max_order=6)
print(f"\nP11 entropy (bits), orders 0-6 (Rao et al. 2010's own block-entropy extension):")
for _k in range(7):
    _row = p11['by_order'][_k]
    print(f"  order{_k}: entropy={_row['entropy']:.3f}  n_contexts={_row['n_contexts']}  n_obs={_row['n_observations']}")
print(f"  within-segment shuffle null (order1): {p11['within_segment_shuffle_order1']['entropy']:.3f} "
      f"(weak null -- keeps each segment's own sign multiset, only randomizes order)")
print("  i.i.d. resample from marginal null, orders 0-6 (correct null -- shares the real corpus's sparsity bias):")
for _k in range(7):
    _row = p11['iid_resample_by_order'][_k]
    print(f"    order{_k}: entropy={_row['entropy']:.3f}  n_contexts={_row['n_contexts']}  n_obs={_row['n_observations']}")
print("  real vs i.i.d.-null gap by order (the genuine-structure signal; Rao et al. (2010) argue this keeps")
print("  scaling past order 1-2 for real language, against Sproat's (2010) order-independent counterexample --")
print("  see run_sproat_counterexample.py for this project's own test of that specific claim):")
for _k in sorted(p11['gap_by_order']):
    print(f"    order{_k}: {p11['gap_by_order'][_k]:+.3f} bits")

p11_boot = tests.p11_bootstrap_ci(word_corpus, orders=(1, 2), n_boot=15)
print(f"  bootstrap check (15 segment resamples w/ replacement -- see p11_bootstrap_ci's docstring:")
print("  sign_stability is the number to trust here, not the percentile CI, which is itself biased")
print("  by bootstrap-induced segment duplication):")
for _k in (1, 2):
    _r = p11_boot[_k]
    print(f"    order{_k}: point={_r['point_gap']:+.3f}  sign_stability={_r['sign_stability']:.3f}  (n_boot={_r['n_boot']})")
