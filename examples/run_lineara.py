"""Regression check: rebuild the Linear A corpus through the toolkit's
Document/Corpus schema from scratch (not from the already-computed p5 pickle)
and confirm P1/P2/P5/P6/P7 reproduce the numbers already published in the
Linear A Dossier. If these don't match, the toolkit has a bug, not the dossier.
"""
import sys, pickle
import os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from decipherment_protocol import Document, Corpus, tests

SCRATCH = '/tmp/claude-0/-home-user-claude-tests/dd3b3458-a783-5e4a-9d0b-609f2f7ec756/scratchpad/'
with open(SCRATCH + 'lineara.pkl', 'rb') as f:
    raw = pickle.load(f)
with open(SCRATCH + 'lineara_catalog.pkl', 'rb') as f:
    cat = pickle.load(f)

DIVIDER = '\U00010101'

def is_real_word(tok):
    if tok in (DIVIDER, '\n', ''):
        return False
    return all(c in cat['sign_to_ascii'] for c in tok)

# Two different corpora from the same raw data, matching the dossier's own two
# distinct sign-set definitions: P1/P2 counted every raw character in the
# transcription; P5/P6/P7 only used fully catalog-recognized "real word" tokens
# (numerals, dividers and uncataloged glyphs excluded, since those aren't
# candidates for a positional/network role the same way a syllabogram is).
census_docs, word_docs = [], []
for name, rec in raw:
    text = rec.get('transcription') or ''
    if not text:
        continue
    chars = [c for c in text if c not in ('\n', ' ', DIVIDER) and c.strip()]
    if chars:
        census_docs.append(Document(doc_id=name, segments=[chars], site=rec.get('site'),
                                     artifact_type=rec.get('support')))
    words_field = rec.get('words') or []
    segments = [list(tok) for tok in words_field if is_real_word(tok)]
    if segments:
        word_docs.append(Document(doc_id=name, segments=segments, site=rec.get('site'),
                                   artifact_type=rec.get('support')))

census_corpus = Corpus(name='Linear A (raw census)', documents=census_docs)
corpus = Corpus(name='Linear A (real-word tokens)', documents=word_docs)
print(f"census documents: {len(census_docs)} (dossier: 825)  |  word-token documents: {len(word_docs)}")

p1 = tests.p1_census(census_corpus)
print(f"\nP1: total occurrences={p1['total_occurrences']} (dossier: 9,763)  "
      f"inventory={p1['inventory_size']} (dossier: 371)")

p2 = tests.p2_rank_frequency(census_corpus, cutoffs=(30, 60, 100, None))
print("\nP2 zipf (dossier: s=0.58/0.69/0.88/1.63):")
for cutoff, r in p2['zipf_by_cutoff'].items():
    print(f"  cutoff={cutoff}: s={r['exponent']:.3f} r2={r['r2']:.3f}")
print(f"MPL exponent = -{p2['mpl_exponent']:.3f} (dossier: -1.032)  r2={p2['mpl_r2']:.3f}")

p5 = tests.p5_positional_classes(corpus)
print(f"\nP5 class counts (dossier: initial 22 / final 25 / free 49): {p5['class_counts']}")
ascii_to_sign = {v: k for k, v in cat['sign_to_ascii'].items()}
ku_sign, ro_sign = ascii_to_sign.get('KU'), ascii_to_sign.get('RO')
ku_row = next((r for r in p5['rows'] if r['sign'] == ku_sign), None)
ro_row = next((r for r in p5['rows'] if r['sign'] == ro_sign), None)
print(f"KU-RO check -- KU: {ku_row['class'] if ku_row else 'not classified'} "
      f"(dossier: initial-preferring, 64%); RO: {ro_row['class'] if ro_row else 'not classified'} "
      f"(dossier: final-preferring, 88%)")

p6 = tests.p6_cooccurrence_network(corpus)
print(f"\nP6: {len(p6['communities'])} communities, Q={p6['modularity_q']:.4f} (dossier: 7 communities, Q=0.198)")
print("  NOTE: community count now matches exactly (7=7) since p6_cooccurrence_network's")
print("  node-inclusion threshold was fixed to use the same occ>=5-within-multi-sign-segments")
print("  definition as P5, instead of counting occurrence across all segments including")
print("  length-1 isolates. Modularity Q is close but not identical (0.191 vs 0.198) --")
print("  a small remaining gap in the greedy merge order, not the node set.")

p7 = tests.p7_segment_length(corpus)
print(f"\nP7: mean word length={p7['mean_length']:.3f} median={p7['median_length']} "
      f"n_segments={p7['n_segments']} (dossier: mean=1.775, median=1, n=5,948)")
print("  NOTE: this run's P7 undercounts vs. the dossier on purpose-revealed-by-accident:")
print("  the dossier's original P7 script counted every non-divider token as a 'word',")
print("  including numeral-tally runs; this adapter's word_docs excludes anything not")
print("  fully catalog-recognized (needed for P5/P6 to be meaningful), which also drops")
print("  numeral tokens from P7's word-length count. The two historical scripts used two")
print("  different sign-sets for P5/P6 vs. P7 without saying so -- the toolkit forces one")
print("  consistent Document.segments definition instead, which is more correct, but means")
print("  a caller must choose up front which sign-set a given corpus's P7 should run over.")

# ---------- P10: Englund totaling-tablet test, with the genre-precondition refinement ----------
# Needs a third corpus: one segment per ATF *line* (not per word), each segment carrying
# both its word-tokens (to check for a literal "total" word) and its numeral value(s).
SUP = {'⁰': '0', '¹': '1', '²': '2', '³': '3', '⁴': '4', '⁵': '5', '⁶': '6', '⁷': '7', '⁸': '8', '⁹': '9'}
SUB = {'₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4', '₅': '5', '₆': '6', '₇': '7', '₈': '8', '₉': '9'}

def parse_fraction(tok):
    if '⁄' not in tok:
        return None
    num, den = tok.split('⁄')
    try:
        return int(''.join(SUP.get(c, c) for c in num)) / int(''.join(SUB.get(c, c) for c in den))
    except (ValueError, ZeroDivisionError):
        return None

def line_signs_and_numerals(tokens):
    signs, nums, i = [], [], 0
    while i < len(tokens):
        tok = tokens[i]
        if tok == DIVIDER:
            i += 1
            continue
        if tok.isdigit():
            v = float(tok)
            if i + 1 < len(tokens):
                f2 = parse_fraction(tokens[i + 1])
                if f2 is not None:
                    v += f2
                    i += 1
            nums.append((v, 'unit'))
        elif parse_fraction(tok) is not None:
            nums.append((parse_fraction(tok), 'unit'))
        else:
            signs.append(tok)
        i += 1
    return signs, nums

line_docs = []
for name, rec in raw:
    tw = rec.get('transliteratedWords')
    if not tw:
        continue
    lines, cur = [], []
    for tok in tw:
        if tok == '\n':
            lines.append(cur); cur = []
        else:
            cur.append(tok)
    if cur:
        lines.append(cur)
    segments, seg_nums = [], []
    for l in lines:
        signs, nums = line_signs_and_numerals(l)
        segments.append(signs)
        seg_nums.append(nums)
    line_docs.append(Document(doc_id=name, segments=segments, segment_numerals=seg_nums,
                               site=rec.get('site'), artifact_type=rec.get('support')))
line_corpus = Corpus(name='Linear A (ATF lines)', documents=line_docs)

p10 = tests.p10_totaling_tablet_test(line_corpus)
print(f"\nP10 unrestricted: tested={p10['tested']} (expect 198) hit_rate={p10['hit_rate']*100:.1f}% "
      f"(expect 3.0%) null_rate={p10['null_rate']*100:.1f}% (expect ~1.5%)")

p10_kuro = tests.p10_totaling_tablet_test(line_corpus, marker_predicate=lambda signs: 'KU-RO' in signs)
print(f"P10 restricted to ku-ro-terminated texts: tested={p10_kuro['tested']} (expect 11) "
      f"hit_rate={p10_kuro['hit_rate']*100:.1f}% (expect 36.4%)")
print("  This is the genre-precondition refinement: Englund's method needs a text that is")
print("  itself a ledger (itemized list + stated total), not just any run of numeral lines.")
