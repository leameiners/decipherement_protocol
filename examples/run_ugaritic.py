"""Ugaritic alphabetic cuneiform as a fourth deciphered-control corpus for the
P1-P10 protocol.

Ugaritic (ca. 1300-1190 BCE, Ras Shamra/Ugarit, Syria) was deciphered in
1930-32 by Bauer, Dhorme, and Virolleaud within months of the tablets'
first publication -- the fastest of this project's four deciphered
controls, and a genuinely different route again: a short (30-sign), fully
alphabetic cuneiform script (one wedge-shape per consonant, plus three
signs for aleph+a/i/u) cracked mostly from cryptanalytic letter-frequency
and positional reasoning on a still-undeciphered but structurally legible
script, the closest real-world precedent to this project's own method.
Ugaritic is also a West Semitic, triliteral-root language, structurally the
most different of the four controls from Linear B (Greek), Ur III
(Sumerian, a language isolate), and Old Persian (Indo-Iranian).

Source: github.com/alexsosn/cuc, auto_parsing/0.2.7/*.tsv -- one file per
KTU-numbered tablet (Cunchillos & Vita's "Concordance of Ugaritic
Texts" numbering), already fully morphologically parsed (surface form,
root/suffix parsing, DULAT dictionary lemma, part of speech, English
gloss) by the dataset's own maintainers. This particular release covers
279 tablets from KTU categories 1-3 only (139 KTU-1 literary/mythological
and ritual texts, 105 KTU-2 letters, 35 KTU-3 legal texts) -- it does not
include the KTU-4 economic/administrative category, a property of this
release, not of the Ugaritic corpus generally.

Data-quality notes (documented, not silently patched over):
  - Many word positions carry more than one competing morphological
    analysis as separate table rows sharing the same numeric word id (e.g.
    id 149475 in KTU 2.50 is parsed both as "to you" and as "lineage") --
    these are alternate readings of one physical occurrence, not two
    occurrences, so rows are deduplicated by id, keeping the first.
  - 401 word ids (1.1%) carry no surface form at all (destroyed/unreadable
    with nothing proposed); dropped, since there is nothing to recover.
  - This dataset marks an individually illegible sign within an otherwise
    legible word with a literal 'x' in the surface-form transliteration
    (e.g. "xgdlt" -> the word is "gdlt", a known word for a cattle
    offering, with its first sign unreadable). 'x' characters are dropped
    and the surviving legible signs of that word are kept, the same
    convention this project uses for Linear B's "[...]" and Old Persian's
    "[damaged]"; 978 words (3.6%) turn out to be wholly illegible ('x'
    only) once this is applied, and are dropped entirely. Across the whole
    corpus, 8.2% of individual sign characters are illegible.
  - Surface-form transliteration is otherwise already a clean, single
    Unicode character per sign (Ugaritic is genuinely alphabetic, unlike
    Old Persian's largely-CV syllabary), so no syllabifier is needed here
    -- each character is used directly as a P5/P6/P7/P9 "sign".

P3 (site specialization) is not attempted: essentially this whole corpus
comes from a single findspot (Ras Shamra), and this source's per-tablet
file carries no site field to build a second axis from regardless -- an
honest omission, not a fabricated one.

P10 (totaling-tablet test) is not run, despite this corpus itself
attesting a real Ugaritic word glossed literally "total (quantity or
price)" (kbd, e.g. KTU 3.13 id 152530) -- a genre precondition problem,
not a missing marker: inspecting KTU 3.13 by hand shows kbd closing each
individual line-item ("w ṯlṯm yn šbˤ kbd d ṯbṭ" = "and 37 [jars of] wine,
total, belonging to Tbt"), i.e. it marks that item's own stated quantity as
already-total, not a final tablet-wide sum over several preceding
itemized lines the way Linear B's to-so or Ur III's szu-nigin2 do. Running
p10_totaling_tablet_test's per-document last-segment-vs-sum-of-earlier-
segments model on that structure would silently misapply the test rather
than genuinely check it, so it is skipped and reported as such instead.
"""
import sys, os, csv, glob
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from decipherment_protocol import Document, Corpus, tests

# clone with: git clone https://github.com/alexsosn/cuc
DATA_DIR = os.environ.get('UGARITIC_DATA_DIR', '/home/user/alexsosn/cuc/auto_parsing/0.2.7')

files = sorted(glob.glob(os.path.join(DATA_DIR, '*.tsv')))

docs = []
total_raw_rows = 0
total_words = 0
total_signs = 0
illegible_signs = 0
dropped_fully_illegible = 0
empty_surface = 0
for path in files:
    siglum = os.path.basename(path)[:-4]
    with open(path, encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader, None)
        seen_ids = set()
        segments = []
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            total_raw_rows += 1
            wid = row[0]
            if wid in seen_ids:
                continue
            seen_ids.add(wid)
            sf = row[1].strip().lower() if len(row) > 1 else ''
            if not sf:
                empty_surface += 1
                continue
            signs = []
            for ch in sf:
                if ch == 'x':
                    illegible_signs += 1
                    continue
                signs.append(ch)
            total_signs += len(sf)
            if signs:
                total_words += 1
                segments.append(signs)
            else:
                dropped_fully_illegible += 1
    if segments:
        docs.append(Document(doc_id=siglum, segments=segments))

print(f"tablet files (KTU 1-3): {len(files)}")
print(f"raw word rows: {total_raw_rows}  distinct words after id-dedup: "
      f"{total_words + dropped_fully_illegible + empty_surface}")
print(f"words with no surface form (fully lost): {empty_surface}")
print(f"words wholly illegible ('x'-only): {dropped_fully_illegible}")
print(f"words with usable signs: {total_words}")
print(f"illegible sign chars dropped: {illegible_signs} / {total_signs} total sign chars "
      f"({illegible_signs/total_signs*100:.1f}%)")

corpus = Corpus(name='Ugaritic (KTU 1-3, auto_parsing)', documents=docs)

p1 = tests.p1_census(corpus)
print(f"\nP1: occurrences={p1['total_occurrences']} inventory={p1['inventory_size']} "
      f"hapax={p1['hapax_legomena']} ({p1['hapax_share']*100:.1f}%)")

p2 = tests.p2_rank_frequency(corpus, cutoffs=(30, 60, 100, None))
print("\nP2 zipf:")
for cutoff, r in p2['zipf_by_cutoff'].items():
    print(f"  cutoff={cutoff}: s={r['exponent']:.3f} r2={r['r2']:.3f}")
print(f"MPL exponent = {p2['mpl_exponent']:.3f}  r2={p2['mpl_r2']:.3f}")
print("  NOTE: as with Old Persian, this corpus's genuinely alphabetic, ~30-sign inventory")
print("  gives too few distinct frequency-of-frequency values for a meaningful MPL fit (r2 is")
print("  undefined here) -- a structural contrast with this project's syllabic/logo-syllabic")
print("  corpora (hundreds of signs), not a bug.")

print("\nP3: not attempted -- see module docstring (single findspot, no per-tablet site field).")

p5 = tests.p5_positional_classes(corpus)
print(f"\nP5 class counts: {p5['class_counts']}")
known = {
    'k': '2sg possessive/object suffix "your"',
    'y': '1sg possessive suffix "my" / nisba adjective ending',
    'n': 'energic/paragogic verbal ending, plural marker',
    'm': 'enclitic mem / masc. plural ending',
    'l': 'common root consonant / preposition "to" (contrast, not predicted final)',
    'b': 'common root consonant / preposition "in" (contrast, not predicted final)',
}
rows_by_sign = {row['sign']: row for row in p5['rows']}
for sign, note in known.items():
    row = rows_by_sign.get(sign)
    if row:
        print(f"  known '{sign}' ({note}): class={row['class']} frac_final={row['frac_final']:.2f} n={row['n']}")
print("""  NOTE: none of the four predicted suffix consonants cross this classifier's 0.6
  final-preferring threshold, unlike Linear B/Ur III/Old Persian's syllable-level
  signs -- but their frac_final values are still the four highest of the six
  (m=0.51, n=0.49, k=0.43, y=0.33) against the two root-consonant contrasts
  (l=0.32, b=0.17), a real directional signal, just a weaker one. The likely reason
  is structural, not a parsing error: Ugaritic's triliteral Semitic roots put a
  genuine root consonant in the word-final slot very often regardless of whether a
  suffix is also present, since almost any consonant can close a bare root -- unlike
  a syllable-final grammatical morpheme in the other three scripts, a final letter
  here is frequently just the third root radical. This is a legitimate limitation of
  a pure single-letter position statistic on a root-and-pattern language, worth
  reporting as such rather than discarding or re-cutting the threshold to force a
  cleaner-looking result.""")

p6 = tests.p6_cooccurrence_network(corpus)
print(f"\nP6: {len(p6['communities'])} communities, Q={p6['modularity_q']:.4f}")

p7 = tests.p7_segment_length(corpus)
print(f"\nP7: mean={p7['mean_length']:.3f} median={p7['median_length']} n={p7['n_segments']}")

community_classes = {s: cid for cid, members in p6['communities'].items() for s in members}
p9 = tests.p9_periodicity(corpus, community_classes, lags=(1, 2, 3), n_perm=2000)
print(f"\nP9 (P6 community label, n_sequences={p9['n_sequences']}):")
for lag, r in p9['lag_results'].items():
    print(f"  lag={lag}: z={r['z']:.2f}")

print("\nP10: not run -- see module docstring (kbd 'total' is attested, but its own use marks")
print("  per-line-item totals, not a tablet-wide closing sum this test's model needs).")

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
