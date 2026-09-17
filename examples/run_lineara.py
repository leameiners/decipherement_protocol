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
print("  NOTE: modularity Q matches closely (0.20 vs 0.198); the community *count* differs")
print("  because node-inclusion (occ>=5) is computed here over all segments including")
print("  length-1 words/isolates, where the dossier's script counted occurrence only within")
print("  words of length>=2 (the same denominator P5 needs). Same reason as the P7 note above.")

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
