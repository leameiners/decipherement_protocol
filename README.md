# decipherment_protocol

A reusable implementation of the P1–P10 statistical protocol developed across
three companion dossiers (Indus, Linear A, Proto-Elamite) for comparing
undeciphered-script corpora on equal footing. Pure Python 3, no dependencies
(no numpy/scipy/networkx) — that was an environment constraint originally,
kept on purpose so the toolkit runs anywhere.

## Why this exists

Each dossier was built by writing one-off scripts per script, per test. That
worked, but it meant every new corpus re-derived the same Zipf fit, the same
modularity algorithm, the same permutation test from scratch — and made it
easy for subtly different definitions to creep in between tests on the same
corpus (see "Known discrepancies" below, found *by* writing this toolkit).
This package is the fourth script's worth of that code, factored out once,
so a fifth script (or a revised P9, or someone else's corpus) doesn't have to
re-derive it again.

## The common schema

Every corpus becomes a `Corpus` of `Document`s (`decipherment_protocol.types`):

```python
from decipherment_protocol import Document, Corpus, tests

docs = [
    Document(doc_id="P008001", segments=[["M157"], ["M319", "M032"]],
             site="Susa", artifact_type="tablet"),
    ...
]
corpus = Corpus(name="my-script", documents=docs)
```

`segments` is the one field every P-test cares about: the text broken into
its natural sub-units, outermost first.

- No known word/line boundary (Indus) → one segment holding the whole text.
- A real word-divider (Linear A) → one segment per word.
- Numbered accounting lines (Proto-Elamite) → one segment per line.

P1/P2/P4/P6 pool signs across all segments; P5/P7/P9 work *within* a segment
(position, length, adjacency). Pick the segmentation that matches what the
corpus actually marks — don't invent word boundaries a script doesn't have.

Numerals that use a distinct sign-set from the main inventory (Proto-Elamite's
N-signs, Linear A's Aegean numeral block) go in `numerals` (flat) and, if P10's
totaling-tablet test is wanted, `segment_numerals` (per-segment, same order as
`segments`).

## Running the protocol

```python
from decipherment_protocol import tests

tests.p1_census(corpus)                        # texts, occurrences, inventory, hapax
tests.p2_rank_frequency(corpus)                 # Zipf + Modified Power Law
tests.p3_site_specialization(corpus)            # checks diversity first; may report not-applicable
tests.p4_allograph_variant_ratio(raw_by_base)   # generic helper only -- see below
tests.p5_positional_classes(corpus)             # initial/medial/final classifier
tests.p6_cooccurrence_network(corpus)           # CNM greedy modularity
tests.p7_segment_length(corpus)                 # word/line length stats
tests.p8_numeral_value_distribution(corpus, numeral_class="N01")
tests.p9_periodicity(corpus, classes=p5_result["classes"])
tests.p10_totaling_tablet_test(corpus)          # Englund-style external-validation helper
```

Each function returns a plain dict — no printing, no plotting. See
`examples/run_lineara.py` and `examples/run_protoelamite.py` for full working
examples that rebuild each corpus from its raw source data and check the
result against the numbers already published in that script's dossier.

## P3's precondition, made explicit

P3 needs real site/type variance to mean anything. `p3_site_specialization`
calls `stats.diversity_ok` first and returns `{"applicable": False, "reason": ...}`
instead of computing meaningless residuals on a near-monoculture table — this
is what happened on Proto-Elamite (94.7% of the accessible corpus is one site,
100% one artifact type). Check `result["applicable"]` before reading `result["residuals"]`.

## P4 and P10 are deliberately thin

What counts as allograph evidence differs by what's actually available for a
script: Indus had glyph-vector SVGs to run shape-similarity clustering on;
Linear A had none, but had a published referent-consolidation table instead;
Proto-Elamite had neither, only subscript-variant notation in the ATF text
itself. `p4_allograph_variant_ratio` only does the one piece that's always
available (how much does stripping a subscript/damage-flag convention compress
the raw token count) — shape clustering or a published-table lookup is the
caller's job, because it's genuinely a different kind of evidence each time.

Similarly, P10 is "check this script's own computed results against whatever
external ground truth is actually accessible for it" — which was an entire
2015 monograph for Indus, one securely-glossed word for Linear A, and a
partially-contested numeral-value literature for Proto-Elamite. The one piece
of P10 that *is* the same test every time — the Englund totaling-tablet
method — is implemented generically as `p10_totaling_tablet_test`; everything
else in a real P10 write-up is corpus-specific research, not a function call.

## P9's honest limitation, carried over

`stats.permutation_lag_test`'s docstring repeats the same caveat every dossier
that has used it has had to state: if the labels fed in were themselves
derived from position (as `p5_positional_classes`'s output is), the lag-1
result is partly circular by construction. This showed up as a different-looking
"clumps then alternates" or "alternates then clumps" pattern on all three
scripts tested so far — genuinely inconsistent findings, which is itself the
finding: **P9 is not yet a valid basis for a cross-script claim** until it's
redefined to use a class label that isn't itself position-based.

## Known discrepancies found while building this (worth reading before trusting a new result)

Validating this toolkit against the Linear A and Proto-Elamite dossiers' own
published numbers surfaced two things the original one-off scripts did
inconsistently, now fixed by forcing one `segments` definition per corpus:

1. **Linear A's original P7 script counted numeral-tally tokens as "words"**
   (5,948 words) while its P5/P6 scripts excluded them (using only fully
   catalog-recognized syllabic/ideographic tokens) — two different sign-sets
   for two tests on the same corpus, undocumented at the time. This toolkit's
   `examples/run_lineara.py` reproduces P1/P2/P5 exactly and explains the P7/P6
   gap in its own output rather than silently forcing a match.
2. **Network node-inclusion (`p6_cooccurrence_network`)** counts sign
   occurrence across *all* segments, including length-1 ones; the original
   per-corpus scripts computed the occurrence count used for the same
   min-frequency threshold only from segments of length ≥2 (the same
   denominator P5 needs). Modularity Q still lands within 0.005–0.04 of the
   original figure on both corpora tested; the exact community count can
   differ because the node set itself is slightly different.

Neither is a correctness bug in the toolkit's math (P1, P2, P5, P7's length
stats, P9, and P10 all reproduce their dossier's numbers exactly once fed an
equivalent `segments` definition) — they're both artifacts of the original
scripts' inconsistency, which building one shared implementation surfaced.
That's the actual case for a toolkit over one-off scripts: not "faster," but
"can't quietly disagree with itself between two tests on the same data."

## Layout

```
decipherment_protocol/
  types.py   Document, Corpus
  stats.py   power_law_fit, contingency_residuals, diversity_ok, cnm_modularity, permutation_lag_test
  tests.py   p1_census ... p10_totaling_tablet_test
examples/
  run_lineara.py        rebuilds the Linear A corpus from raw source, validates against the dossier
  run_protoelamite.py   same, for Proto-Elamite -- also exercises P3, P9, P10
  run_indus.py          same idea for Indus, but see the caveat below -- it doesn't fully land
```

## Validation status, honestly

- **Linear A**: P1, P2, and P5 reproduce the dossier's published numbers exactly.
  P6/P7 are close (Q within 0.005; length within 0.1) with the node-set/word-set
  discrepancy above fully diagnosed, not hand-waved.
- **Proto-Elamite**: P1, P2, P5, P7, P10 reproduce exactly; P3 correctly detects
  non-applicability; P6 is close (Q within 0.005); P9's z-scores are close
  (different permutation draw, as expected).
- **Indus**: does NOT reproduce exactly. `texts.js`'s `text_code` field carries
  ad hoc annotation conventions (bracket-enclosed partial signs, slash-separated
  ambiguous alternate readings) this adapter approximates rather than fully
  resolves -- inventory comes out 677 vs. the dossier's 715, occurrences ~15,000
  vs. 18,069. The Zipf/MPL exponents and the site-specialization signal are in
  the right range and point the same direction, which is enough to show the
  toolkit's math isn't broken on a third, differently-shaped corpus, but this
  is the one example still worth treating as approximate rather than verified.
  Resolving the `text_code` parsing precisely is the concrete next step before
  claiming all three scripts, not two of three, run through one shared,
  fully-validated implementation.
