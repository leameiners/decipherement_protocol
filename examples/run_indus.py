"""Regression check against the Indus Script Dossier's published numbers.

Indus has no word-divider and no numeral sub-alphabet distinct from the main
inventory (numeral candidates are just low-numbered signs, 001-013), so each
Document here is the simplest possible shape: one segment holding the whole
text_code sign sequence, exactly matching how P1/P2/P6's original scripts
treated it. P5's segment-internal position, in this one-segment-per-text
setup, becomes "position within the whole text" -- the same thing the
dossier's original F/M/T positional-histogram script measured.

Unlike the Linear A and Proto-Elamite examples, this one does NOT reproduce
the dossier's numbers exactly (inventory 677 vs. 715, occurrences ~15,000 vs.
18,069) -- text_code carries ad hoc annotation conventions (bracket-enclosed
partial signs, slash-separated ambiguous alternate readings, a literal "000"
damage placeholder) that this adapter approximates rather than fully resolves,
and it's possible the original dossier numbers came partly from icit.js's own
independently-tallied per-sign totals rather than a pure text_code re-parse.
Left in on purpose, as an honestly-reported gap rather than a silently-dropped
example: it demonstrates the toolkit's math produces sane, right-shaped
numbers on Indus too (Zipf exponents and MPL land within ~0.03-0.04 of
published, site specialization detects the same Harappa/Mohenjo-daro
signature), without claiming a reproduction this adapter doesn't actually
deliver. Fixing the text_code parsing to match exactly is the concrete
next step before this counts as a real third validated corpus.
"""
import sys, json
sys.path.insert(0, '/tmp/claude-0/-home-user-claude-tests/dd3b3458-a783-5e4a-9d0b-609f2f7ec756/scratchpad/protocol_toolkit')
from decipherment_protocol import Document, Corpus, tests

with open('/tmp/texts.js', encoding='utf-8') as f:
    content = f.read()
content = content[content.index('{'):].rstrip().rstrip(';')
raw = json.loads(content)
cols = raw['columns']
idx = {c: i for i, c in enumerate(cols)}

def parse_text_code(tc):
    if not tc:
        return []
    out = []
    for tok in tc.strip('+').split('-'):
        tok = tok.strip('[]')
        if not tok:
            continue
        if '/' in tok:
            tok = tok.split('/')[0]  # ambiguous alternate reading -- take the first
        if tok and tok != '000' and tok.isdigit():
            out.append(tok)
    return out

docs = []
for row in raw['rows']:
    signs = parse_text_code(row[idx['text_code']])
    if not signs:
        continue
    docs.append(Document(doc_id=str(row[idx['id']]), segments=[signs],
                          site=row[idx['site']] or None, artifact_type=row[idx['type']] or None))

corpus = Corpus(name='Indus', documents=docs)
print(f"documents with a parseable text_code: {len(docs)} (dossier: 5,445 total records)")

p1 = tests.p1_census(corpus)
print(f"\nP1: occurrences={p1['total_occurrences']} (dossier: 18,069)  "
      f"inventory={p1['inventory_size']} (dossier: 715)  mean_len={p1['mean_text_length']:.2f} "
      f"(dossier: ~4.5)  max_len={p1['max_text_length']} (dossier: ~26)")

p2 = tests.p2_rank_frequency(corpus, cutoffs=(50, None))
print("\nP2 zipf (dossier: s=0.74 top-50, s=1.55 full):")
for cutoff, r in p2['zipf_by_cutoff'].items():
    print(f"  cutoff={cutoff}: s={r['exponent']:.3f} r2={r['r2']:.3f}")
print(f"MPL exponent = -{p2['mpl_exponent']:.3f} (dossier: -1.207 ours / -1.347 Fuls published)  "
      f"r2={p2['mpl_r2']:.3f}")

p3 = tests.p3_site_specialization(corpus)
print(f"\nP3 applicable: {p3['applicable']}")
if p3['applicable']:
    key_pairs = [('Harappa', 'TAB:B'), ('Harappa', 'TAB:I'), ('Mohenjo-daro', 'SEAL:S')]
    for site, atype in key_pairs:
        r = p3['residuals'].get((site, atype))
        if r:
            print(f"  {site} x {atype}: obs={r['observed']} exp={r['expected']:.1f} "
                  f"std_resid={r['std_resid']:.1f} adj_resid={r['adj_resid']:.1f}")
    print("  (dossier's Harappa x tablets figure pooled TAB:B+TAB:I+TAB:C as one 'tablets' "
          "category and Mohenjo-daro x seals pooled SEAL:S+SEAL:R+SEAL:C+SEAL:CY as one "
          "'seals' category; this run uses the raw un-pooled type codes, so exact figures "
          "won't match §16's pooled adj_resid=32.1/26.3 without the same pooling step -- "
          "shown here to demonstrate the function runs and returns sane numbers, not to "
          "reproduce that specific pooled figure.)")

p7 = tests.p7_segment_length(corpus)
print(f"\nP7 (segment = whole text, since there's no word-divider): mean={p7['mean_length']:.2f} "
      f"(dossier frames this as text length ~4.5, not word length -- Indus word length is a "
      f"separate inferred estimate the dossier gets from Fuls' connectivity formula, not from "
      f"this function)")
