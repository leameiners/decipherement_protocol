"""Old Persian cuneiform as a third deciphered-control corpus for the P1-P10
protocol.

Old Persian was the first cuneiform script deciphered (Grotefend 1802,
completed by Rawlinson and others by the 1840s-50s), using multilingual
royal inscriptions (Old Persian alongside Elamite and Babylonian versions of
the same text) -- a genuinely different decipherment route from Linear B's
single-script statistical/positional attack, and from Ur III Sumerian's
gradual 19th-century cuneiform decipherment via Akkadian. It is also
structurally the most different of this project's three deciphered
controls: a small, largely CV-syllabic script (a few dozen signs, not
hundreds) invented specifically for royal Achaemenid inscriptions, used for
monumental/political text, not administrative ledgers.

Source: github.com/Electronic-Old-Persian-Library/Old-Persian-Dataset,
textdata/web_scraping/ -- 92 files (90 real inscriptions + the scraper
script and its readme), one per standard Kent siglum (DNa, DSf, XPh, A1Pb,
...), scraped from livius.org's Achaemenid Royal Inscriptions pages. Of the
90, 75 contain actual transliterated text; the other 15 (including DB, the
Behistun inscription) are livius.org pages that only link out to further
sub-pages the scraper didn't follow, so they carry no romanized text at all
in this source -- a property of the source, not of this adapter, reported
the same way this project reports every other corpus's real coverage gap.

Parsing notes (documented, not silently patched over):
  - livius.org marks the scribal word-divider with a literal '\\' in the
    page text; words sometimes wrap across a source line break with no
    space (mid-word), and the transliteration itself is interleaved with
    an English translation paragraph-by-paragraph on the same page. Taking
    only the lines that contain '\\' and concatenating them in file order
    with no inserted separator correctly reassembles words that wrap across
    an English paragraph (verified by hand against DNa: "...thâtiy \\ D" +
    "ârayavauš \\ ..." -> "Dârayavauš", the correct king's name, spanning a
    translation paragraph in between) -- this can, in principle, also
    spuriously fuse the last word of one paragraph with the first of an
    unrelated one on pages with more than one physically distinct text;
    not corrected for, and expected to affect at most one word boundary
    per file.
  - A handful of damage annotations ("[damaged]", "[lacuna]", "[broken
    off]", "[...]") are dropped; a handful of partially-damaged words
    carrying a single stray bracket (e.g. "Suguda]iyam") have the bracket
    stripped and the letters kept, the same convention as this project's
    other corpora for paleographically-restored-but-legible text.
  - This source gives continuous romanization with no sign-boundary mark of
    its own (unlike Linear B's and Ur III's hyphenated transliterations),
    so `graphemes()` approximates the actual cuneiform sign sequence with a
    greedy CV syllabifier (consonant, or the digraphs 'xš'/'th' for the
    single phonemes they represent, absorbing an immediately-following
    vowel into one sign; a word-final consonant or an unabsorbed vowel
    stands alone). This is a documented approximation, not the
    epigraphically attested sign list, but matches the standard scholarly
    sign division on every form checked by hand, e.g. xšâyathiya -> xšâ-ya-
    thi-ya, dârayavauš -> dâ-ra-ya-va-u-š.

P3 (site specialization) is not attempted: this source's per-page scrape
carries no structured site/artifact-type metadata beyond free-text page
boilerplate that isn't reliably parseable across 90 differently-formatted
pages -- an honest omission, not a fabricated axis.

P10 (totaling-tablet test) does not apply to this corpus at all: these are
monumental royal proclamations, not administrative ledgers, so there is no
itemized-list-plus-stated-total structure for the test to find -- the same
genre precondition this project's "Other candidates considered" section
already uses to rule out Egyptian hieroglyphic monumental texts.
"""
import sys, os, re, glob
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from decipherment_protocol import Document, Corpus, tests

# clone with: git clone https://github.com/Electronic-Old-Persian-Library/Old-Persian-Dataset
DATA_DIR = os.environ.get(
    'OLD_PERSIAN_DATA_DIR',
    '/home/user/electronic-old-persian-library/old-persian-dataset/textdata/web_scraping')

ANNOTATION_RE = re.compile(r'\[(damaged|lacuna|broken off|\.\.\.)\]', re.IGNORECASE)
DIGRAPHS = ('xš', 'th')
VOWELS = set('aiuâîûāūeo')


def parse_words(blob):
    blob = ANNOTATION_RE.sub(' ', blob)
    blob = blob.replace('[', '').replace(']', '')
    words = []
    for tok in re.split(r'[\s:;]+', blob):
        tok = tok.strip()
        if not tok or tok == '...':
            continue
        words.append(tok)
    return words


def graphemes(word):
    """Greedy CV syllabifier approximating this word's Old Persian cuneiform
    sign sequence -- see module docstring for the method and its caveats."""
    w = word.lower()
    signs, i, n = [], 0, len(w)
    while i < n:
        two = w[i:i + 2]
        if two in DIGRAPHS:
            cons, i = two, i + 2
        elif w[i] not in VOWELS and w[i].isalpha():
            cons, i = w[i], i + 1
        else:
            cons = None
        if cons is not None:
            if i < n and w[i] in VOWELS:
                signs.append(cons + w[i]); i += 1
            else:
                signs.append(cons)
        else:
            if w[i].isalpha():
                signs.append(w[i])
            i += 1
    return signs


docs = []
files = sorted(glob.glob(os.path.join(DATA_DIR, '*.txt')))
files_with_content = 0
total_words = 0
for path in files:
    siglum = os.path.basename(path)[:-4]
    with open(path, encoding='utf-8', errors='replace') as f:
        lines = [l.rstrip('\n') for l in f if '\\' in l]
    if not lines:
        continue
    blob = ''.join(lines)
    words = parse_words(blob)
    segments = [graphemes(w) for w in words]
    segments = [s for s in segments if s]
    if segments:
        files_with_content += 1
        total_words += len(segments)
        docs.append(Document(doc_id=siglum, segments=segments))

print(f"inscription files: {len(files) - 2}  with transliterated text: {files_with_content}")
print(f"total words: {total_words}")

corpus = Corpus(name='Old Persian (Kent sigla, web_scraping)', documents=docs)

p1 = tests.p1_census(corpus)
print(f"\nP1: occurrences={p1['total_occurrences']} inventory={p1['inventory_size']} "
      f"hapax={p1['hapax_legomena']} ({p1['hapax_share']*100:.1f}%)")

p2 = tests.p2_rank_frequency(corpus, cutoffs=(30, 60, 100, None))
print("\nP2 zipf:")
for cutoff, r in p2['zipf_by_cutoff'].items():
    print(f"  cutoff={cutoff}: s={r['exponent']:.3f} r2={r['r2']:.3f}")
print(f"MPL exponent = {p2['mpl_exponent']:.3f}  r2={p2['mpl_r2']:.3f}")
print("  NOTE: sign given as-is (not forced negative, unlike this project's other adapters) --")
print("  this corpus's much smaller, more alphabet-like 118-sign inventory (vs. hundreds for")
print("  Linear A/B, Proto-Elamite, Ur III) gives a flatter, less power-law-shaped")
print("  frequency-of-frequency distribution, so the fit direction itself is a real, reportable")
print("  structural contrast with this project's other corpora, not a bug to hide.")

print("\nP3: not attempted -- see module docstring (no reliable site/artifact-type "
      "metadata in this source).")

p5 = tests.p5_positional_classes(corpus)
print(f"\nP5 class counts: {p5['class_counts']}")
known = {
    'm': 'accusative singular ending (-am) -- should be final-preferring',
    'y': '3sg verb ending -tiy / adjective ending -iy -- should be final-preferring',
    'š': 'nominative singular masc. ending (-uš) -- should be final-preferring',
    'a': 'default/inherent vowel -- should not be final-preferring',
}
rows_by_sign = {row['sign']: row for row in p5['rows']}
for sign, note in known.items():
    row = rows_by_sign.get(sign)
    if row:
        print(f"  known '{sign}' ({note}): class={row['class']} frac_final={row['frac_final']:.2f} n={row['n']}")

p6 = tests.p6_cooccurrence_network(corpus)
print(f"\nP6: {len(p6['communities'])} communities, Q={p6['modularity_q']:.4f}")

p7 = tests.p7_segment_length(corpus)
print(f"\nP7: mean={p7['mean_length']:.3f} median={p7['median_length']} n={p7['n_segments']}")

community_classes = {s: cid for cid, members in p6['communities'].items() for s in members}
p9 = tests.p9_periodicity(corpus, community_classes, lags=(1, 2, 3), n_perm=2000)
print(f"\nP9 (P6 community label, n_sequences={p9['n_sequences']}):")
for lag, r in p9['lag_results'].items():
    print(f"  lag={lag}: z={r['z']:.2f}")

print("\nP10: not applicable -- see module docstring (monumental genre, no totaling-tablet "
      "structure to test).")

# ---------- P11: conditional entropy (Rao et al. 2009-style; see README's
# "Raghavendra (2026) and the entropy test" section) ----------
p11 = tests.p11_conditional_entropy(corpus, max_order=6)
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

p11_boot = tests.p11_bootstrap_ci(corpus, orders=(1, 2), n_boot=200)
print(f"  bootstrap check (200 segment resamples w/ replacement -- see p11_bootstrap_ci's docstring:")
print("  sign_stability is the number to trust here, not the percentile CI, which is itself biased")
print("  by bootstrap-induced segment duplication):")
for _k in (1, 2):
    _r = p11_boot[_k]
    print(f"    order{_k}: point={_r['point_gap']:+.3f}  sign_stability={_r['sign_stability']:.3f}  (n_boot={_r['n_boot']})")
