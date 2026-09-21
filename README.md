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

The package has also been run against two *deciphered* corpora as a
retrospective check — does a blind classifier's output agree with what
decipherment already established? See "Deciphered controls" below.

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
tests.p10_totaling_tablet_test(corpus, marker_predicate=lambda signs: "KU-RO" in signs)
                                                 # restricted to texts that actually close with a "total" marker
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

## The totaling-tablet test needs a genre precondition, not just numeral lines

`p10_totaling_tablet_test` takes an optional `marker_predicate(signs: list[str]) -> bool`,
restricting the test to documents whose *last* numeral-bearing segment's sign
list satisfies it — e.g. `lambda signs: 'KU-RO' in signs`. Without this,
the test conflates real ledgers (an itemized list followed by a stated total)
with ordinary multi-line inventories that never state a total at all, which
dilutes any real signal toward the null. Testing this directly:

- **Linear A, unrestricted** (every text with 3+ numeral-bearing lines,
  198 texts): 3.0% hit rate vs. a 1.5% null — close to chance, because most
  of those 198 texts simply aren't ledgers with a summary line.
- **Linear A, restricted to texts literally closing with a *ku-ro* ("total")
  line** (11 texts): 36.4% — a real jump, and three of the seven non-matches
  (HT9a, HT13, HT102) are independently tagged "Wrong Total" by the source
  corpus's own metadata, meaning the check finds real scribal errors, not noise.
- **Proto-Elamite, restricted to texts closing with sign M288** (the
  documented totalizer, 102 texts): 33.3% vs. a 22.8% null — only a modest
  improvement over the unrestricted rate (17.7% vs. 13.1%), nothing like
  Linear A's jump.
- **Indus, restricted to its one genre with an actual decoded numeral
  system** (V+# volumetric tablets — sign 700 followed by a long-linear-stroke
  numeral): nothing at all. `segment` here isn't a line within one text —
  Indus has no line structure — so the check groups *co-located artifacts*
  by fine-grained findspot instead (site/area/section/block/house/room) and
  asks whether any object's V+# value equals the sum of the others found with
  it. Only 14 findspot groups have enough decodable V+# objects to test:
  0% hit rate, 0% null. A more careful leave-one-out version across all 26
  groups with 2+ decodable objects finds just 4 apparent matches — every one
  two objects trivially sharing the same small value (4=4, 3=3) — against a
  shuffle-based null that produces *more* such coincidences on average (8.09)
  than the real data does.

The asymmetry across all three is the actual finding. M288's documented
function ("a container sign that may function as a unit marker or totalizer")
is weaker and more hedged than *ku-ro*'s dedicated meaning "total" — a sign
that sometimes plays a totalizing role is not the same restriction as a word
whose sole job is summation. And Indus's V+# genre fails for a more basic
reason than either: a V+# tablet records one container's fill count, never
several counts that a related object then sums into a stated total — there's
no ledger structure to restrict *to* in the first place, even in Indus's most
numerically legible genre. The genre precondition Englund's method needs is
about the marker's function and the text's structure, not just a plausible
sign showing up in the right place. See `examples/run_lineara.py`,
`examples/run_protoelamite.py`, and `examples/run_indus.py` for the exact
reproductions of all three results above.

## P9, redefined

`stats.permutation_lag_test`'s docstring states a caveat worth repeating here:
if the labels fed into it are themselves derived from position (as
`p5_positional_classes`'s output is), the lag-1 result is partly circular by
construction — nearby positions in a short segment correlate on a
position-derived label almost by definition. Running P9 on `p5`'s classes
produced a different-looking pattern on each of the three corpora (Indus:
signs alternate class at lag 1; Linear A: signs clump at lag 1, then
alternate at lags 2–3; Proto-Elamite: signs clump hard at lag 1, sit at
chance at lag 2, then alternate hard at lag 3) — three inconsistent shapes
using an identical method, which was itself the finding: not a valid basis
for a cross-script claim as originally defined.

`p9_periodicity` already accepted an arbitrary `sign -> label` mapping, so no
core change was needed to fix this, only which test supplies the label.
`examples/run_p9_redefined.py` reruns P9 using `p6_cooccurrence_network`'s
community assignment instead — a label built from which other signs a sign
co-occurs with across the whole corpus, not from its own position — which
removes the circularity. The result: every script clumps at lag 1 (Indus
z=+25.0, Linear A z=+2.2, Proto-Elamite z=+9.6), i.e. signs from the same
co-occurrence community tend to sit next to each other within a segment,
strongly on two scripts and weakly on the third. This is now a real, if
modest, generalizing result across all three corpora, not a discarded one.

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
sign sequences — resolves 82.9% of that total (14,978 occurrences). This is
not a parsing gap: the field's character set is closed (digits, `+`, `-`,
`[`, `]`, `/`), and every convention it uses is handled exactly — a leading
or trailing bracket marks partial legibility on an otherwise-legible sign; a
literal `000` is an explicit damage placeholder with nothing left to read;
`a/b` gives two alternate readings, and the parser takes whichever one isn't
itself `000`. The remaining 17.1% is confirmed missingness in the source:
902 of 5,445 records carry no `text_code` at all, and 1,442 further
occurrences within the remaining records are `000` damage placeholders.
P1/P2 in `examples/run_indus.py` use `icit.js` directly
(`p1_census_from_freq`/`p2_rank_frequency_from_freq`); P3/P5/P6/P7 still need
`texts.js`'s per-text sequences and are reported against that 82.1%-of-records
subset's actual, quantified size rather than treated as the whole corpus. To
regenerate these inputs from scratch:

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
  run_lineara.py         rebuilds the Linear A corpus from raw source
  run_protoelamite.py    same, for Proto-Elamite -- also exercises P3, P9, P10
  run_indus.py           same, for Indus -- exact match on P1/P2, subset-based on P3/P5/P6/P7
  run_p9_redefined.py    reruns P9 on all three corpora using P6 community labels, not P5's
  run_linearb.py          deciphered control: Linear B (Mycenaean Greek), full P1-P10
  run_urIII.py             deciphered control: Ur III Sumerian, full P1-P10
```

## Validation status

- **Linear A**: P1, P2, P5, and P6's community count reproduce their expected
  published numbers exactly (7 communities). P6's modularity Q is close but
  not identical (0.191 vs. 0.198 — a small remaining gap in the greedy merge
  order, not the node set, which is now defined identically to P5's). P7 is
  close (length within 0.1) because the original per-corpus analysis used two
  subtly different sign-set definitions for P5/P6 versus P7 (whether
  numeral-tally tokens count as "words") that this package's single
  consistent `segments` definition doesn't reproduce bit-for-bit — documented
  in the example's own output rather than silently forced to match.
- **Proto-Elamite**: P1, P2, P5, P6 (community count *and* modularity, both
  exact — 5 communities, Q=0.1433), P7, and P10 all reproduce exactly; P3
  correctly detects non-applicability; P9's z-scores are close (different
  permutation draw, as expected from a stochastic test).
- **Indus**: P1 and P2 reproduce exactly (18,069 occurrences, 715 signs, 216
  hapax, Zipf s=0.736/1.553, Modified Power Law exponent −1.207) once fed
  `icit.js`'s own frequency table rather than a re-parse of `texts.js`. P3,
  P5, P6, and P7 run on the 82.1%-resolvable `text_code` subset described
  above; P6 lands close to expected (6 communities, Q=0.141) since Indus's
  one-segment-per-text structure means there are effectively no length-1
  isolate segments to create a node-set mismatch in the first place.

P6's node-inclusion threshold (`p6_cooccurrence_network`) counts sign
occurrence only within segments of length ≥2 — the same definition
`p5_positional_classes` uses (`_occ_in_multi_sign_segments`), since a sign
occurring only in length-1 segments has no position information and can
never form a co-occurrence edge either. This is what fixed Linear A's and
Proto-Elamite's community counts to match exactly; the original version of
this package thresholded on raw occurrence across all segments, which
counted isolates that the source analyses this package validates against
had already excluded.

## A fourth script isn't a code problem, it's a data problem

The protocol is built to be corpus-agnostic -- `types.Document`/`Corpus` only
need a sign sequence per text, an optional sub-unit boundary, and site/type
metadata -- so applying it to a fourth script needs no new code, only a new
adapter like `run_indus.py`/`run_lineara.py`/`run_protoelamite.py`. What
actually blocks a fourth application in practice is that every realistic
undeciphered or partially-deciphered candidate we could identify -- Cypro-Minoan
(217 catalogued inscriptions), Linear Elamite (a corpus split across a Susa
sub-corpus of stone monuments and a Collection sub-corpus of silver vessels,
smaller still: 265 sign variants averaging only 5.6 occurrences each, against
905 variants averaging 19.9 for Indus), Rongorongo (roughly two dozen
surviving objects) -- is one to two orders of magnitude smaller than any of
the three corpora this package already validates against, and none has an
openly downloadable transliterated dataset reachable through this package's
usual approach (clone a mirror, fetch a bulk dump). A future adapter for any
of these should expect thin, noisy P2/P5/P6/P7/P9 results as a property of
the source, not a bug in the adapter.

## Deciphered controls

None of the four corpora above has a known answer to check the protocol's
classifiers against. `run_linearb.py` does: Linear B (Mycenaean Greek),
deciphered by Ventris in 1952, is a direct sibling of Linear A -- same
administrative-ledger genre, many shared sign shapes -- except its grammar
and hundreds of individual word meanings are now scholarly consensus. Source:
`github.com/InsiderPhD/Linear-B-Dataset`, 4,794 distinct tablets after
de-duplicating the raw scrape (789 of 796 duplicated identifiers were
byte-identical repeats).

P5's purely distributional classifier, given no semantic information, sorts
four well-documented Mycenaean grammatical particles: *-de* (allative,
"to/towards"), *-qe* (enclitic "and"), and *-jo* (genitive-singular ending)
all land final-preferring, matching their documented function exactly. The
fourth, *wa-*, lands free/medial -- also correctly, since *wa* is both a real
final inflectional element and one of the commonest syllables to *begin* a
Mycenaean name (*wa-na-ka*, "king"), so a mixed classification is the
linguistically accurate answer, not a miss. P10's genre-restricted totaling
test, run against Linear B's own summation particle *to-so* ("so much, so
many"), finds the same ledger-plus-stated-total structure Linear A's *ku-ro*
does, including one exact hand-verified match (tablet KN As 1517: seventeen
itemized entries, *to-so* VIR 17); the automated hit rate undercounts this
because the source's flattened line format doesn't reliably mark where a
tablet's lines break or whether it continues after a stated subtotal, not
because the structure isn't there. P1's blind census also recovers 90 signs,
close to Linear B's known ~87-sign core syllabary. P2's Modified Power Law
exponent (-0.26) is a genuine outlier against this package's other corpora
(-0.82 to -1.21): a closed phonetic syllabary transliterated at the syllable
level is a flatter, more alphabet-like distribution than the logogram-heavy
sign lists P2 was otherwise calibrated against, a scope difference rather
than a failed prediction.

`run_urIII.py` runs a second, independent check against Ur III Sumerian (ca.
2100-2000 BCE) -- the single largest administrative cuneiform archive that
survives, and the direct genre ancestor of this protocol's P10 (Damerow and
Englund developed the totaling-tablet method on archives exactly like this
one). Source: `github.com/cdli-gh/data`, the same official CDLI bulk dump
this package's Proto-Elamite adapter already uses, filtered to period "Ur III
(ca. 2100-2000 BC)" and genre "Administrative": 106,802 catalogued texts, of
which 74,678 (69.9%) carry a transliteration in the bulk release.

P1 counts 7,044 distinct signs across 4.62 million occurrences -- Sumerian's
transliterated "signs" at this level are closer to distinct word-readings
than a small closed syllabary, so both the inventory size and the hapax
share (43.8%) run well above this package's other corpora for that reason,
not because the tail is genuinely longer. P2's exponent (-0.82) is a third
distinct position on the same scale as Linear B's -0.26 and the original
three corpora's -1.03-to-1.35 range -- consistent with the exponent tracking
how word-like versus sign-like a script's transliterated units are, as much
as it tracks language typology. P3 does not apply, but for an adapter reason
rather than a corpus one: restricting to a single genre (Administrative)
leaves no second artifact-type axis to test site against. P5 again produces
the clearest check: two case suffixes with one grammatical function each
(*-sze3* "to/for", *-ta* "from") land final-preferring at 0.77 and 0.86 --
higher confidence than any single Linear B ending -- while the content noun
*lugal* ("king") is correctly classified free/medial rather than an ending.
A third sign, *-ka* ("in/of"), lands free/medial too, mirroring Linear B's
*wa-* complication: *-ka* is also the ordinary word for "mouth," and the
syllable recurs inside many unrelated word stems. P6 recovers 7 communities
at Q=0.291, the strongest modularity this package has found on any corpus.
P7's mean word length (4.4 signs) is longer than every other corpus here by
a wide margin, consistent with Sumerian's agglutinative morphology chaining
several affixes onto one word. P9, run on an 8,000-document subsample for
tractability (n_perm=200; the full corpus's ~90,000+ eligible segments make
the default n_perm=2000 impractically slow in pure Python), clumps sharply
at lag 1 (z=+115.9) and alternates beyond it, the same shape found
everywhere else. P10's totaling test restricted to *szu-nigin* ("grand
total") finds a real signal (7.1% hit rate against a 2.3% null on 98 texts)
smaller in absolute terms than Linear A's 36.4%, but Ur III's *unrestricted*
rate already beats its own null by a comparable margin (5.3% vs. 1.9%,
n=31,436) -- administrative ledgers are common with or without an explicit
total-word, unlike Linear A's mixed-genre corpus. The absolute rate is also
a likely undercount specifically: this adapter decodes only Sumerian's base
counting system, not the separate grain/area/weight systems many Ur III
accounts total in.

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
