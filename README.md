# decipherment_protocol

A dependency-free Python implementation of P1–P10: ten statistical tests for
comparing undeciphered-script corpora on equal footing (rank-frequency
fitting, site/artifact-type specialization, positional structural
classification, sign co-occurrence networks, segment-length measurement,
numeral-system analysis, periodicity, and external validation against
published scholarship). No numpy/scipy/networkx — every primitive (OLS power-law
fits, contingency-table residuals, greedy modularity maximization, permutation
testing) is implemented from scratch in `decipherment_protocol/stats.py`, so
the package runs anywhere a bare Python 3 interpreter does.

It was developed by running the identical ten tests on three real corpora —
the Indus Valley script, Linear A, and Proto-Elamite — three scripts whose
decipherment status differs in three genuinely different ways (Indus: sounds
unknown, structure statistically inferred; Linear A: ~100 signs carry a
scholarly-consensus sound value via shared Linear B shapes, but the language
they spell is unidentified; Proto-Elamite: numeral *values* were substantially
cracked in the 1980s via pure arithmetic cross-checking, while its ~600
ideographic signs and underlying language remain unread). Running one battery
of tests across three corpora with three different decipherment shapes is
what let several of these tests be validated as genuinely useful rather than
artifacts of one dataset — see "What generalized and what didn't" below.

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

`segments` is the one field every test cares about: the text broken into its
natural sub-units, outermost first.

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
tests.p1_census_from_freq(sign_freq)            # same, from a bare {sign: count} table -- no per-text data needed
tests.p2_rank_frequency_from_freq(sign_freq)    # same, when a script's own maintainers publish frequencies directly
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
`examples/run_lineara.py`, `examples/run_protoelamite.py`, and
`examples/run_indus.py` for full working examples that rebuild each corpus
from its raw source data and check the results for internal consistency.

## P3's precondition, made explicit

P3 needs real site/type variance to mean anything. `p3_site_specialization`
calls `stats.diversity_ok` first and returns `{"applicable": False, "reason": ...}`
instead of computing meaningless residuals on a near-monoculture table. This
is not a hypothetical: on the accessible Proto-Elamite corpus, 94.7% of texts
come from a single site (Susa) and 100% are the same artifact type (clay
tablets), so `p3_site_specialization` correctly refuses to compute residuals
rather than producing numbers with no interpretive content. Check
`result["applicable"]` before reading `result["residuals"]`.

## P4 and P10 are deliberately thin

What counts as allograph evidence differs by what's actually available for a
script. For Indus, glyph-vector SVGs exist and support pixel/shape-similarity
clustering. Linear A has none, but the LinearA Explorer project (built on
Godart & Olivier's *GORILA* corpus) publishes a referent-consolidation table
mapping distinct sign-shapes to a shared commodity meaning. Proto-Elamite has
neither, only ATF subscript-variant notation on individual signs.
`p4_allograph_variant_ratio` only does the one piece that's always available —
how much a subscript/damage-flag stripping convention compresses the raw
token count — leaving shape clustering or a published-table lookup as the
caller's job, because it is genuinely different evidence each time.

Similarly, P10 means "check this script's own computed results against
whatever external ground truth is actually accessible for it," and that
literature's size differs by roughly an order of magnitude across the three
scripts tested (Section "What generalized and what didn't"). The one piece of
P10 that *is* the same test every time — Englund's totaling-tablet
cross-check, the actual method that helped decipher proto-cuneiform and
Proto-Elamite numeral *values* (Englund 2004; Damerow & Englund 1989) — is
implemented generically as `p10_totaling_tablet_test`; everything else in a
real P10 analysis is corpus-specific research, not a function call.

## P9's honest limitation

`stats.permutation_lag_test`'s docstring states a caveat worth repeating here:
if the labels fed into it were themselves derived from position (as
`p5_positional_classes`'s output is), the lag-1 result is partly circular by
construction — nearby positions in a short segment correlate on a
position-derived label almost by definition. Running this test on all three
corpora produced a different-looking pattern each time (Indus: signs
alternate class at lag 1; Linear A: signs clump at lag 1, then alternate at
lags 2–3; Proto-Elamite: signs clump hard at lag 1, sit at chance at lag 2,
then alternate hard at lag 3) — three inconsistent shapes using an identical
method, which is itself the finding: **P9 is not yet a valid basis for a
cross-script claim** until it is redefined around a class label that isn't
itself position-based.

## What generalized and what didn't

Running the same ten tests on Indus, Linear A, and Proto-Elamite (built from
a third-party mirror of Wells & Fuls' *Interactive Corpus of Indus Texts*, the
LinearA Explorer's digitization of Godart & Olivier's *GORILA*, and the
official CDLI bulk-data dump respectively) produced a mix of results that
looked like genuine cross-script signal and results that turned out to be
script-specific or not yet trustworthy:

- **P2 (rank-frequency).** Fuls' Modified Power Law calibration places
  fusional/syllabic languages near an exponent of −1.1 and logo-syllabic
  Classic Maya near −1.4 (Fuls, in Wells 2015, Appendix III). Linear A (a real
  syllabary) landed at −1.03, near the syllabic pole as expected; Indus (mixed
  logo-syllabic) landed at −1.21, correctly between the two poles. Proto-Elamite
  — purely ideographic/numerical, no syllabic layer at all — was predicted to
  land beyond Indus, further from the syllabic pole; it instead landed at
  −1.10, closer to the syllabic pole than Indus. Two out of three following a
  known-typology prediction is not the same as the test generalizing — a third,
  genuinely out-of-sample script showed the pattern can break.
- **P3 (site specialization)** was extreme wherever it could be tested
  (Indus: Harappa/tablets; Linear A: Haghia Triada/nodules), with different
  sites and artifact types each time — consistent with a general property of
  Bronze Age administrative record-keeping rather than a civilization-specific
  pattern — and correctly detected as inapplicable on Proto-Elamite.
- **P5 (positional classification)**, built from pure sign-position statistics
  with zero semantic input, was independently confirmed against real external
  ground truth on every script tested: Indus sign 700 as the volumetric "V"
  marker (confirmed in Fuls 2024, matching a purely computational
  identification made without access to that paper); Linear A's KU/RO split
  matching the initial/final roles of the two syllables in *ku-ro*, the one
  Linear A word with a near-universally accepted meaning ("total," an
  identification attributed to M. Pope and widely cited in the field);
  Proto-Elamite's sign M288 (this corpus's single most frequent sign) coming
  out 94% final-preferring, matching its documented function as a container
  sign/totalizer.
- **P6 (co-occurrence networks)** found real, non-random structure on all
  three scripts (modularity Q = 0.137–0.198, weak-to-moderate by conventional
  thresholds), though of different character — Linear A's network cleanly
  isolated its entire numeral sub-system, while Indus's and Proto-Elamite's
  communities look more like genre overlap, except for one Indus community
  that blindly recovered the documented Harappa tablet/volumetric-accounting
  tradition.
- **P7 (segment length)** converged tightly across three unrelated
  measurement routes: Indus's two indirect estimates (a rank-frequency
  regression and a Fuls-style connectivity segmentation) bracket 1.81–2.16
  signs; Linear A's direct word-divider count gives 1.775; Proto-Elamite's
  direct accounting-line count gives 1.784. Three scripts, three measurement
  methods, one narrow sub-2-sign band.
- **P8 (numeral systems)** showed a real, checkable difference rather than a
  repeated finding: Indus's tally signs (001–007) break the cross-linguistic
  rule that number-word frequency decreases monotonically from 1 upward
  (Dehaene & Mehler 1992) — sign 002 is used nearly four times more than sign
  001, a pattern absent from every natural language surveyed in that study —
  while Proto-Elamite's base numeral sign (N01) obeys the same rule cleanly
  across values 1–9.
- **P10 (external validation)** is not comparable at face value across
  scripts: Indus has an entire monograph with technical appendices to check
  against (Wells 2015, cross-referencing Fuls' independent 2023–24 work),
  Linear A has essentially one securely-glossed word and some catalog
  metadata, and Proto-Elamite's numeral-value literature is itself internally
  contested (a mainstream sexagesimal-system value and a disputed
  reinterpretation for the same sign class disagree by 36×). The Englund
  totaling-tablet method specifically was underpowered on Indus's very short
  inscriptions (91% of texts carry 0–1 numeral-candidate signs) but produced a
  real, modest signal on Proto-Elamite's longer accounting texts (17.7% hit
  rate vs. a 13.1% permutation-shuffled null) — the corpus that method was
  historically built for.

## Corpora used to build and validate this package

| Script | Source | Scale |
| --- | --- | --- |
| Indus | A third-party GitHub mirror of Wells & Fuls' *Interactive Corpus of Indus Texts* (ICIT) | 5,445 records, 715 signs, 18,069 occurrences |
| Linear A | `github.com/mwenge/lineara.xyz` (LinearA Explorer), built from Godart & Olivier's *GORILA* | 1,722 catalogued (825 with transcription), 371 raw shapes, 9,763 occurrences |
| Proto-Elamite | `github.com/cdli-gh/data` (official CDLI bulk dump), filtered to Proto-Elamite | 1,558 tablets, 608 base signs, 14,711 occurrences |

None of these is a verified live authoritative database. The Indus mirror
disagrees with the primary published source in a documented, characterizable
way: `icit.js` (a per-sign aggregate frequency table shipped with the mirror)
sums to exactly 18,069 occurrences across exactly 715 signs, while
`texts.js`'s per-text `text_code` field — the only source with real per-text
sign sequences — resolves only 82.1% of that total even after handling its
bracket-enclosed partial-sign and slash-separated ambiguous-reading
conventions. P1/P2 in `examples/run_indus.py` use `icit.js` directly
(`p1_census_from_freq`/`p2_rank_frequency_from_freq`); P3/P5/P6/P7 still need
`texts.js`'s per-text sequences and are reported against that 82% subset's
actual, quantified size rather than treated as the whole corpus. To regenerate
these inputs from scratch:

- **Linear A**: clone `github.com/mwenge/lineara.xyz` and parse
  `LinearAInscriptions.js` — a JS `Map` literal, not a plain array, with two
  `Map`s in the same file (find the *first* closing `]);` after the opening
  `new Map([`, not the last) and ES6 `\u{XXXXX}` extended Unicode escapes
  (invalid strict JSON — convert with a regex before `json.loads`). `ideograms.js`
  and `network/commodities.js` give the published sign catalog and commodity-variant
  groups used for P4/P8/P10.
- **Proto-Elamite**: clone `github.com/cdli-gh/data`, then fetch the two
  Git-LFS-pointer files directly via
  `https://media.githubusercontent.com/media/cdli-gh/data/master/<file>`
  (works without a git-lfs client). Filter `cdli_cat.csv`'s `period` field for
  "Proto-Elamite" to get P-numbers, then pull the matching `&P######` ATF
  blocks out of `cdliatf_unblocked.atf`.
- **Indus**: `texts.js` and `icit.js` came from a third-party GitHub mirror of
  the ICIT database; see the References below for the primary published
  source (Wells 2015) this mirror is derived from.

## Layout

```
decipherment_protocol/
  types.py   Document, Corpus
  stats.py   power_law_fit, contingency_residuals, diversity_ok, cnm_modularity, permutation_lag_test
  tests.py   p1_census ... p10_totaling_tablet_test
examples/
  run_lineara.py        rebuilds the Linear A corpus from raw source
  run_protoelamite.py   same, for Proto-Elamite -- also exercises P3, P9, P10
  run_indus.py          same, for Indus -- exact match on P1/P2, subset-based on P3/P5/P6/P7
```

## Validation status

- **Linear A**: P1, P2, and P5 reproduce their expected published numbers
  exactly. P6/P7 are close (modularity Q within 0.005; length within 0.1)
  because the original per-corpus analysis used two subtly different sign-set
  definitions for P5/P6 versus P7 (whether numeral-tally tokens count as
  "words") that this package's single consistent `segments` definition
  doesn't reproduce bit-for-bit — documented in each example's own output
  rather than silently forced to match.
- **Proto-Elamite**: P1, P2, P5, P7, and P10 reproduce exactly; P3 correctly
  detects non-applicability; P6 is close (Q within 0.005); P9's z-scores are
  close (different permutation draw, as expected from a stochastic test).
- **Indus**: P1 and P2 reproduce exactly (18,069 occurrences, 715 signs, 216
  hapax, Zipf s=0.736/1.553, Modified Power Law exponent −1.207) once fed
  `icit.js`'s own frequency table rather than a re-parse of `texts.js`. P3,
  P5, P6, and P7 run on the 82.1%-resolvable `text_code` subset described
  above; P6's community count and modularity now land close to expected
  (6 communities, Q=0.143) since Indus's one-segment-per-text structure
  sidesteps the isolate-counting issue noted for Linear A above.

## References

- Wells, B.K. (2015). *The Archaeology and Epigraphy of Indus Writing*.
  Oxford: Archaeopress. Technical appendices (segmentation algorithm,
  Modified Power Law, Relative Information) by Andreas Fuls.
- Fuls, A. (2023). *A Catalog of Indus Signs*. Mathematica Epigraphica 4.
- Fuls, A. (2024). "Comparison of Linear Elamite and Indus Writing Systems."
  *Iranian Journal of Archaeological Studies* 14(1).
- Joshi, J.P. and Parpola, A., eds. (1987). *Corpus of Indus Seals and
  Inscriptions, 1: Collections in India*. Helsinki: Suomalainen Tiedeakatemia.
- Rao, R.P.N., Yadav, N., Vahia, M.N., Joglekar, H., Adhikari, R., and
  Mahadevan, I. (2009). "Entropic Evidence for Linguistic Structure in the
  Indus Script." *Science* 324(5931).
- Sproat, R. (2010). "Ancient symbols, computational linguistics, and the
  reviewing practices of the general science journals." *Computational
  Linguistics* 36(3). Sproat, R. (2014). "A statistical comparison of written
  language and nonlinguistic symbol systems." *Language* 90(2).
- Farmer, S., Sproat, R., and Witzel, M. (2004). "The Collapse of the
  Indus-Script Thesis: The Myth of a Literate Harappan Civilization."
  *Electronic Journal of Vedic Studies* 11(2).
- Godart, L. and Olivier, J.-P. *GORILA: Recueil des inscriptions en linéaire A*.
  Études Crétoises, 5 volumes.
- Damerow, P. and Englund, R.K. (1989). *The Proto-Elamite Texts from Tepe
  Yahya*. American School of Prehistoric Research Bulletin 39.
- Friberg, J. (1978–79). *The Third Millennium Roots of Babylonian
  Mathematics*, I–II. University of Gothenburg.
- Englund, R.K. (2004). "The State of Decipherment of Proto-Elamite." In S.
  Houston, ed., *The First Writing: Script Invention as History and Process*.
  Cambridge University Press.
- Dehaene, S. and Mehler, J. (1992). "Cross-linguistic regularities in the
  frequency of number words." *Cognition* 43(1).
