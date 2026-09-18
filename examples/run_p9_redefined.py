"""Redefine P9 periodicity around P6's co-occurrence-network community label
instead of P5's within-segment position class.

P5's initial/final/free class is derived from *where in the segment* a sign
sits, so testing periodicity of P5 labels along the segment is partly
circular by construction (a lag-1 result is guaranteed to some degree by the
same positional information that produced the label). P6's community label
is derived from *which other signs a sign co-occurs with in the same
segment*, independent of position, so a periodicity test over the sequence
of a segment's community memberships is not circular in the same way.

Rebuilds all three corpora from the same cached raw/pickled data the other
regression scripts use, runs p6_cooccurrence_network to get community
labels, then runs p9_periodicity with those labels instead of p5's classes.
"""
import sys, json, pickle
import os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from decipherment_protocol import Document, Corpus, tests

SCRATCH = '/tmp/claude-0/-home-user-claude-tests/dd3b3458-a783-5e4a-9d0b-609f2f7ec756/scratchpad/'

def load_js_object(path):
    with open(path, encoding='utf-8') as f:
        content = f.read()
    content = content[content.index('{'):].rstrip().rstrip(';')
    return json.loads(content)


def run_for(name, corpus, lags=(1, 2, 3), n_perm=2000):
    p6 = tests.p6_cooccurrence_network(corpus)
    community_classes = {s: cid for cid, members in p6['communities'].items() for s in members}
    print(f"\n=== {name} ===")
    print(f"P6: {len(p6['communities'])} communities, Q={p6['modularity_q']:.4f}, "
          f"{len(community_classes)} labeled signs")
    p9 = tests.p9_periodicity(corpus, community_classes, lags=lags, n_perm=n_perm)
    print(f"P9 (community label, n_sequences={p9['n_sequences']}):")
    for lag, r in p9['lag_results'].items():
        print(f"  lag={lag}: z={r['z']:.2f}")
    return p9


# ---------- Indus ----------
icit = load_js_object('/tmp/icit.js')
raw = load_js_object('/tmp/texts.js')
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
            alts = [a.strip('[]') for a in tok.split('/')]
            legible = [a for a in alts if a != '000']
            tok = legible[0] if legible else alts[0]
        if tok and tok != '000' and tok.isdigit():
            out.append(tok)
    return out

indus_docs = []
for row in raw['rows']:
    signs = parse_text_code(row[idx['text_code']])
    if not signs:
        continue
    indus_docs.append(Document(doc_id=str(row[idx['id']]), segments=[signs],
                                site=row[idx['site']] or None, artifact_type=row[idx['type']] or None))
indus_corpus = Corpus(name='Indus (structural subset)', documents=indus_docs)
run_for('Indus', indus_corpus, lags=(1,))  # Indus segment = whole text, no internal lag>1 possible per P5 script

# ---------- Linear A ----------
with open(SCRATCH + 'lineara.pkl', 'rb') as f:
    la_raw = pickle.load(f)
with open(SCRATCH + 'lineara_catalog.pkl', 'rb') as f:
    la_cat = pickle.load(f)

DIVIDER = '\U00010101'

def is_real_word(tok):
    if tok in (DIVIDER, '\n', ''):
        return False
    return all(c in la_cat['sign_to_ascii'] for c in tok)

la_word_docs = []
for name, rec in la_raw:
    words_field = rec.get('words') or []
    segments = [list(tok) for tok in words_field if is_real_word(tok)]
    if segments:
        la_word_docs.append(Document(doc_id=name, segments=segments, site=rec.get('site'),
                                      artifact_type=rec.get('support')))
la_corpus = Corpus(name='Linear A (real-word tokens)', documents=la_word_docs)
run_for('Linear A', la_corpus, lags=(1, 2, 3))

# ---------- Proto-Elamite ----------
with open(SCRATCH + 'pe_parsed.pkl', 'rb') as f:
    D = pickle.load(f)
cat, parsed = D['cat'], D['parsed']
pe_docs = []
for pnum, lines in parsed.items():
    segments = [l['signs'] for l in lines]
    seg_nums = [l['nums'] for l in lines]
    meta = cat[pnum]
    pe_docs.append(Document(doc_id=pnum, segments=segments,
                             site=meta['provenience'], artifact_type=meta['object_type'],
                             segment_numerals=seg_nums))
pe_corpus = Corpus(name='Proto-Elamite', documents=pe_docs)
run_for('Proto-Elamite', pe_corpus, lags=(1, 2, 3))
