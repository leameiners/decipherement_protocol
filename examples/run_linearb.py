"""Linear B as a deciphered control corpus for the P1-P10 protocol.

Linear B (Mycenaean Greek, ca. 1400-1200 BCE) was deciphered by Michael
Ventris in 1952. It is the direct sibling script of Linear A already used
elsewhere in this project: same sign inventory in spirit (many Linear A signs
carry a phonetic value inherited by Linear B), same administrative-ledger
genre, same word-divider convention, digitized from the same third-party
scraped-mirror tradition as this project's other corpora. Because the
language, the grammar, and hundreds of individual word meanings are known
with scholarly consensus, running the protocol on Linear B is not a fourth
open application -- it's a retrospective check: does the protocol's blind
statistical machinery recover signals that agree with what decipherment
actually revealed?

Source: github.com/InsiderPhD/Linear-B-Dataset, itself scraped from
minoan.deaditerranean.com's Linear B transliterations. Not a verified live
authoritative database, same caveat as this project's other three corpora.

Data-quality notes, found and handled the same way as this project's other
adapters (documented, not silently patched over):
  - 7,370 raw rows resolve to only 4,794 distinct tablet identifiers; 789 of
    796 duplicated identifiers are byte-identical repeats (a scrape artifact,
    not real re-editions), so rows are deduplicated by identifier, keeping
    the first occurrence.
  - The source mixes uppercase and lowercase transliteration for the same
    signs across records (e.g. "QE-RO2" in some tablets, "qe-ro2" in others)
    -- normalized to lowercase throughout, since case carries no phonemic
    distinction in Linear B transliteration.
  - "(so)" marks a paleographically uncertain but proposed reading -- kept,
    parentheses stripped, same convention as Indus's bracket-marked partial
    signs. "[...]" marks genuinely illegible/missing material -- dropped.
  - "*56"-style tokens are uncatalogued syllabograms without an agreed
    phonetic value, referenced by their conventional sign number -- kept as
    an opaque sign identifier, the same way Indus signs are referenced by
    number rather than by putative sound value.
"""
import sys, csv, re
import os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from decipherment_protocol import Document, Corpus, tests

# clone with: git clone https://github.com/InsiderPhD/Linear-B-Dataset
DATA = os.environ.get('LINEAR_B_DATA', '/home/user/insiderphd/linear-b-dataset/tablet-sets/tablets.csv')


def is_ideogram_code(token):
    """A bare all-uppercase, unhyphenated token (VIR, GRA, OVIS, TELA1...) is
    a logogram/commodity-ideogram transliteration, not a syllabic sign --
    excluded from the sign inventory the same way Linear A's ideogram and
    commodity signs are kept out of its P5/P6/P7 'real word token' subset.
    Some records in this source transliterate ordinary multi-syllable WORDS
    in uppercase too (e.g. "QE-RO2" for "qe-ro2"), so the hyphen is the
    signal that distinguishes a real word from a single ideogram code."""
    letters = re.sub(r'[^A-Za-z]', '', token)
    return bool(letters) and letters.isupper() and '-' not in token

def parse_word(token):
    if not token:
        return None
    token = token.strip().lower()
    if not token:
        return None
    signs = []
    for part in token.split('-'):
        part = part.strip()
        if not part:
            continue
        if '[' in part or ']' in part:
            continue  # illegible -- nothing to recover
        part = part.strip('()')  # uncertain-but-proposed reading; sign is otherwise known
        if part:
            signs.append(part)
    return signs or None


def dedupe(rows):
    seen = {}
    for r in rows:
        if r['identifier'] not in seen:
            seen[r['identifier']] = r
    return list(seen.values())


with open(DATA, encoding='utf-8') as f:
    raw_rows = list(csv.DictReader(f, delimiter=';'))
rows = dedupe(raw_rows)
print(f"raw rows: {len(raw_rows)}  distinct tablets after de-duplication: {len(rows)}")

# ---------- word-level corpus: P1, P2, P3, P5, P6, P7, P9 ----------
word_docs = []
total_word_tokens = 0
resolved_word_tokens = 0
for r in rows:
    words = [w for w in r['inscription'].split(',')]
    total_word_tokens += sum(1 for w in words if w.strip())
    segments = []
    for w in words:
        signs = parse_word(w)
        if signs:
            segments.append(signs)
            resolved_word_tokens += 1
    if segments:
        word_docs.append(Document(doc_id=r['identifier'], segments=segments,
                                   site=r['location'] or None, artifact_type=r['series'] or None))

word_corpus = Corpus(name='Linear B (word tokens)', documents=word_docs)
print(f"word-bearing tablets: {len(word_docs)} / {len(rows)}")
print(f"word tokens resolved: {resolved_word_tokens} / {total_word_tokens} "
      f"({resolved_word_tokens/total_word_tokens*100:.1f}%) -- the rest are wholly illegible ('[...]')")

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
    for (site, series), r in top:
        print(f"  {site} x {series}: obs={r['observed']} exp={r['expected']:.1f} adj_resid={r['adj_resid']:.1f}")

p5 = tests.p5_positional_classes(word_corpus)
print(f"\nP5 class counts: {p5['class_counts']}")
# ground-truth cross-check: known Linear B grammatical endings
known_endings = {'qe': 'enclitic "and" -- should be final-preferring', 'jo': 'common genitive-singular ending -- final-leaning',
                  'wa': 'common nominative/feminine ending', 'de': 'directional/allative particle -- should be final'}
rows_by_sign = {row['sign']: row for row in p5['rows']}
for sign, note in known_endings.items():
    row = rows_by_sign.get(sign)
    if row:
        print(f"  known ending '{sign}' ({note}): class={row['class']} frac_final={row['frac_final']:.2f} n={row['n']}")

p6 = tests.p6_cooccurrence_network(word_corpus)
print(f"\nP6: {len(p6['communities'])} communities, Q={p6['modularity_q']:.4f}")

p7 = tests.p7_segment_length(word_corpus)
print(f"\nP7: mean={p7['mean_length']:.3f} median={p7['median_length']} n={p7['n_segments']}")

community_classes = {s: cid for cid, members in p6['communities'].items() for s in members}
p9 = tests.p9_periodicity(word_corpus, community_classes, lags=(1, 2, 3), n_perm=2000)
print(f"\nP9 (P6 community label, n_sequences={p9['n_sequences']}):")
for lag, r in p9['lag_results'].items():
    print(f"  lag={lag}: z={r['z']:.2f}")

# ---------- line-level corpus: P10 ----------
# Many longer tablets in this source carry no ".1"/".2"/".a" line markers at
# all (they were apparently only kept for short 2-3 line records), so a
# marker-based split misses most of the genre entirely. Ledger tablets have a
# real structure regardless: an entry is a run of word signs followed by the
# number counted for it ("qa-ra-jo 1", "a-nu-wi-ko VIR 1", "to-so VIR 31"), so
# segmenting on numeral boundaries recovers each entry -- and the closing
# total -- directly from that structure instead of from markers the source
# doesn't reliably carry.
NUMERAL = re.compile(r'^\(?(\d+)\)?$')
# line-number markers (".1", ".a"), continuation slashes, and the editorial
# Latin annotations this source carries inline ("vacat" = blank space left on
# the tablet, "sup./inf. mut." = top/bottom broken off) -- none are sign
# content.
NON_SIGN = re.compile(r'^\.\w*$|^/$|^,$')
ANNOTATION_WORDS = {'vacat', 'mut', 'sup', 'inf', 'vest'}

def parse_original(text):
    segments = []
    cur_signs, cur_nums = [], []
    for tok in text.split():
        if is_ideogram_code(tok) or NON_SIGN.match(tok) or tok.strip('.').lower() in ANNOTATION_WORDS:
            continue
        m = NUMERAL.match(tok)
        if m:
            cur_nums.append((float(m.group(1)), 'unit'))
            segments.append((cur_signs, cur_nums))
            cur_signs, cur_nums = [], []
            continue
        ps = parse_word(tok)
        if ps:
            cur_signs.extend(ps)
    if cur_signs or cur_nums:
        segments.append((cur_signs, cur_nums))
    return segments

line_docs = []
for r in rows:
    if not r['original']:
        continue
    parsed = parse_original(r['original'])
    if not parsed:
        continue
    segments = [signs for signs, nums in parsed]
    seg_nums = [nums for signs, nums in parsed]
    line_docs.append(Document(doc_id=r['identifier'], segments=segments, segment_numerals=seg_nums,
                               site=r['location'] or None, artifact_type=r['series'] or None))
line_corpus = Corpus(name='Linear B (ATF-style lines)', documents=line_docs)
print(f"\nline-level documents: {len(line_docs)}")

p10 = tests.p10_totaling_tablet_test(line_corpus)
print(f"P10 unrestricted: tested={p10['tested']} hit_rate={p10['hit_rate']*100:.1f}% null_rate={p10['null_rate']*100:.1f}%")

def has_toso(signs):
    # parse_word splits "to-so" into its component syllables ['to', 'so'],
    # so the marker is a consecutive pair, not a single joined token.
    return any(signs[i] == 'to' and signs[i + 1] == 'so' for i in range(len(signs) - 1))

p10_toso = tests.p10_totaling_tablet_test(line_corpus, marker_predicate=has_toso)
print(f"P10 restricted to to-so ('so much/so many' -- Linear B's own totaling word): "
      f"tested={p10_toso['tested']} hit_rate={p10_toso['hit_rate']*100:.1f}% null_rate={p10_toso['null_rate']*100:.1f}%")
print("""  NOTE: this source's flattened line format (no reliable line breaks on
  longer tablets, no marker for 'tablet continues after this subtotal')
  makes the generic hit-rate an undercount, not a null result. Three
  examples checked by hand against the parsed segments directly:
    KN As 1517: 17 itemized VIR-1 entries, to-so VIR 17 -- exact match.
    KN As 1519: 8 parsed entries before to-so VIR 10 -- the tablet has a
      further unlabeled entry (\"ma-ri-ne-wo wo-i-ko-de\", no explicit
      count) folded into the total line by this parser, not a real miss.
    KN As 1516: 30 entries against to-so ... VIR 31 -- one entry's own
      count falls inside a damaged span ('[ ]'), so it parses as present
      but uncounted; short by exactly one, not a scribal error.
  The ledger-plus-stated-total structure P10 is built to find is visibly
  there in every one of these; what falls short of Linear A's and
  Proto-Elamite's precision is this adapter's numeral alignment on a
  source that does not mark line boundaries or breaks consistently, not
  the tablets themselves.""")

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
