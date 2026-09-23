"""P1-P10: the decipherment-statistics protocol itself, as functions over a Corpus.

Each function takes a decipherment_protocol.types.Corpus and returns a plain
dict of results -- no side effects, no printing, so callers can format results
however a given dossier needs. P4 and P10 are deliberately thin: what counts
as "allograph evidence" or "external validation" is corpus-specific by nature
(shape geometry vs. published referent tables; a whole monograph vs. one
glossed word), so those two only provide the generic helper a caller needs,
not a one-shot answer.
"""
import random
from collections import Counter, defaultdict
from . import stats


def p1_census(corpus):
    freq = corpus.sign_freq()
    lengths = [len(d.flat_signs) for d in corpus.documents]
    hapax = sum(1 for f in freq.values() if f == 1)
    lengths_sorted = sorted(lengths)
    n = len(lengths_sorted)
    return {
        'n_texts': len(corpus.documents),
        'total_occurrences': sum(freq.values()),
        'inventory_size': len(freq),
        'hapax_legomena': hapax,
        'hapax_share': hapax / len(freq) if freq else float('nan'),
        'mean_text_length': sum(lengths) / n if n else float('nan'),
        'median_text_length': lengths_sorted[n // 2] if n else float('nan'),
        'max_text_length': max(lengths) if lengths else 0,
        'sign_freq': dict(freq),
    }


def p1_census_from_freq(sign_freq: dict):
    """P1 restricted to what a bare frequency table can answer: total
    occurrences, inventory size, hapax rate. No text-length stats, since those
    need per-text sequences a frequency table alone doesn't carry.

    Use this instead of p1_census when a script's own maintainers publish an
    aggregate per-sign frequency table that is more authoritative than
    anything a caller could re-derive from parsing raw per-text sequences --
    e.g. Indus's icit.js, whose 715-sign/18,069-occurrence totals are the
    dossier's actual source for these two numbers, not a re-parse of
    texts.js's text_code field (see the Indus example and its README note on
    why the two disagree).
    """
    hapax = sum(1 for f in sign_freq.values() if f == 1)
    return {'total_occurrences': sum(sign_freq.values()), 'inventory_size': len(sign_freq),
            'hapax_legomena': hapax, 'hapax_share': hapax / len(sign_freq) if sign_freq else float('nan'),
            'sign_freq': dict(sign_freq)}


def p2_rank_frequency_from_freq(sign_freq: dict, cutoffs=(30, 60, 100, None)):
    """P2 on a bare frequency dict -- see p1_census_from_freq's docstring for
    when to prefer this over p2_rank_frequency(corpus, ...)."""
    ranked = sorted(sign_freq.values(), reverse=True)
    zipf = {}
    for c in cutoffs:
        label = c if c is not None else len(ranked)
        s, r2 = stats.power_law_fit(ranked, c)
        zipf[label] = {'exponent': s, 'r2': r2}
    g, r2_mpl = stats.modified_power_law(sign_freq)
    return {'zipf_by_cutoff': zipf, 'mpl_exponent': -g, 'mpl_r2': r2_mpl}


def p2_rank_frequency(corpus, cutoffs=(30, 60, 100, None)):
    return p2_rank_frequency_from_freq(corpus.sign_freq(), cutoffs)


def p3_site_specialization(corpus, min_dominant_share=0.85):
    table = Counter()
    for d in corpus.documents:
        if d.site and d.artifact_type:
            table[(d.site, d.artifact_type)] += 1
    ok, reason = stats.diversity_ok(dict(table), max_dominant_share=min_dominant_share)
    if not ok:
        return {'applicable': False, 'reason': reason, 'table': dict(table)}
    residuals = stats.contingency_residuals(dict(table))
    return {'applicable': True, 'table': dict(table), 'residuals': residuals}


def p4_allograph_variant_ratio(raw_tokens_by_base: dict[str, set]):
    """Generic P4 helper: how much does a base/raw sign-identity choice compress
    the inventory? Feed it {base_sign: {raw_token_1, raw_token_2, ...}} built
    however the caller strips subscripts/damage-flags for a given script's ATF
    or catalog convention. Shape-similarity clustering (Indus) or published
    referent-grouping (Linear A) both still need to be done by the caller --
    this only measures the cheap, always-available part of P4.
    """
    n_base = len(raw_tokens_by_base)
    n_raw = sum(len(v) for v in raw_tokens_by_base.values())
    multi = {b: v for b, v in raw_tokens_by_base.items() if len(v) > 1}
    return {'n_base_signs': n_base, 'n_raw_tokens': n_raw,
            'raw_to_base_ratio': n_raw / n_base if n_base else float('nan'),
            'n_base_with_multiple_variants': len(multi),
            'variant_share': len(multi) / n_base if n_base else float('nan')}


def _occ_in_multi_sign_segments(corpus):
    """Occurrence count restricted to segments of length >= 2.

    Shared by p5_positional_classes and p6_cooccurrence_network so both use
    the same node-eligibility threshold. A sign occurring only in length-1
    segments (isolates) has no position information at all (no "initial" or
    "final" is meaningful for it) and can never contribute a co-occurrence
    edge either, since there is nothing else in that segment to pair with --
    so excluding those occurrences from the frequency threshold both tests
    gate on is the correct behavior for each, not just a shared quirk.
    """
    occ = Counter()
    for d in corpus.documents:
        for seg in d.segments:
            if len(seg) < 2:
                continue
            for s in seg:
                occ[s] += 1
    return occ


def p5_positional_classes(corpus, min_n=5, threshold=0.6):
    pos_counts = defaultdict(lambda: {'initial': 0, 'medial': 0, 'final': 0})
    total_occ = Counter()
    for d in corpus.documents:
        for seg in d.segments:
            n = len(seg)
            if n < 2:
                continue
            for i, s in enumerate(seg):
                total_occ[s] += 1
                if i == 0:
                    pos_counts[s]['initial'] += 1
                elif i == n - 1:
                    pos_counts[s]['final'] += 1
                else:
                    pos_counts[s]['medial'] += 1
    classes, rows = {}, []
    for s, occ in total_occ.items():
        if occ < min_n:
            continue
        pc = pos_counts[s]
        fi, fm, ff = pc['initial'] / occ, pc['medial'] / occ, pc['final'] / occ
        cls = 'initial-preferring' if fi >= threshold else ('final-preferring' if ff >= threshold else 'free/medial')
        classes[s] = cls
        rows.append({'sign': s, 'n': occ, 'frac_initial': fi, 'frac_medial': fm, 'frac_final': ff, 'class': cls})
    rows.sort(key=lambda r: -r['n'])
    return {'classes': classes, 'rows': rows, 'class_counts': dict(Counter(r['class'] for r in rows))}


def fsw_formulaic_repetition(corpus, lengths=(3, 4, 5, 6)):
    """Farmer, Sproat & Witzel (2004)-style formulaic-repetition rate, as
    used in Nair (2026)'s synthetic-baseline scorecard (see this package's
    README's "Nair (2026)" discussion): at each phrase length L, the share
    of distinct L-sign phrases that recur in two or more different
    documents ("inscriptions"), among all distinct L-sign phrases found at
    all. A phrase's own repeat occurrences WITHIN one document count once
    toward that document's contribution -- only recurrence ACROSS
    documents counts toward "formulaic," matching the real-language
    intuition of a stock phrase reused across different texts, not a
    single text's own internal repetition (already covered by this
    package's P9/entropy tests). Phrases are drawn from each document's
    own flat sign sequence, not crossing document boundaries.

    Returns {L: {'total_types', 'repeated_types', 'repetition_rate'}}.
    """
    out = {}
    for L in lengths:
        phrase_docs = defaultdict(set)
        for i, d in enumerate(corpus.documents):
            signs = d.flat_signs
            seen_in_doc = {tuple(signs[j:j + L]) for j in range(len(signs) - L + 1)}
            for phrase in seen_in_doc:
                phrase_docs[phrase].add(i)
        total_types = len(phrase_docs)
        repeated_types = sum(1 for docs in phrase_docs.values() if len(docs) >= 2)
        out[L] = {'total_types': total_types, 'repeated_types': repeated_types,
                   'repetition_rate': repeated_types / total_types if total_types else float('nan')}
    return out


def fsw_positional_rigidity(corpus, top_n=10, min_n=5):
    """Farmer, Sproat & Witzel (2004)-style positional-rigidity score, as
    used in Nair (2026)'s synthetic-baseline scorecard: Cramer's V (see
    stats.cramers_v) between "is this sign" and "position within its
    segment" (initial/medial/final, the same convention p5_positional_
    classes uses), for each of the corpus's top-N most frequent signs,
    tested against a 2x3 table of {this sign, every other sign} x
    {initial, medial, final}. Higher V means that sign sticks to one
    position more rigidly than the corpus's signs do on average.

    Returns {'per_sign': {sign: V}, 'mean_v': ..., 'top_signs': [...]}
    -- mean_v is the single positional-rigidity score for the corpus,
    matching what Nair (2026) reports per corpus.
    """
    pos_counts = defaultdict(lambda: {'initial': 0, 'medial': 0, 'final': 0})
    total_occ = Counter()
    for d in corpus.documents:
        for seg in d.segments:
            n = len(seg)
            if n < 2:
                continue
            for i, s in enumerate(seg):
                total_occ[s] += 1
                if i == 0:
                    pos_counts[s]['initial'] += 1
                elif i == n - 1:
                    pos_counts[s]['final'] += 1
                else:
                    pos_counts[s]['medial'] += 1
    corpus_totals = {pos: sum(pc[pos] for pc in pos_counts.values()) for pos in ('initial', 'medial', 'final')}
    top_signs = [s for s, _ in sorted(total_occ.items(), key=lambda kv: -kv[1]) if total_occ[s] >= min_n][:top_n]
    per_sign = {}
    for s in top_signs:
        pc = pos_counts[s]
        table = {}
        for pos in ('initial', 'medial', 'final'):
            table[('sign', pos)] = pc[pos]
            table[('rest', pos)] = corpus_totals[pos] - pc[pos]
        per_sign[s] = stats.cramers_v(table)
    mean_v = sum(per_sign.values()) / len(per_sign) if per_sign else float('nan')
    return {'per_sign': per_sign, 'mean_v': mean_v, 'top_signs': top_signs}


def p6_cooccurrence_network(corpus, min_n=5):
    total_occ = _occ_in_multi_sign_segments(corpus)
    nodes = sorted([s for s, n in total_occ.items() if n >= min_n])
    node_idx = {s: i for i, s in enumerate(nodes)}
    edge_weight = defaultdict(float)
    for d in corpus.documents:
        for seg in d.segments:
            idxs = [node_idx[s] for s in seg if s in node_idx]
            for a in range(len(idxs)):
                for b in range(a + 1, len(idxs)):
                    i, j = idxs[a], idxs[b]
                    if i == j:
                        continue
                    if i > j:
                        i, j = j, i
                    edge_weight[(i, j)] += 1.0
    communities, q = stats.cnm_modularity(nodes, dict(edge_weight))
    named = {c: [nodes[i] for i in mem] for c, mem in communities.items()}
    return {'nodes': nodes, 'n_edges': len(edge_weight), 'communities': named, 'modularity_q': q}


def p7_segment_length(corpus):
    lengths = [len(seg) for d in corpus.documents for seg in d.segments if len(seg) > 0]
    lengths_sorted = sorted(lengths)
    n = len(lengths_sorted)
    return {'n_segments': n, 'mean_length': sum(lengths) / n if n else float('nan'),
            'median_length': lengths_sorted[n // 2] if n else float('nan'),
            'length_distribution': dict(Counter(lengths))}


def p8_numeral_value_distribution(corpus, numeral_class: str, value_range=range(1, 10)):
    """Generic P8 helper: usage-skew of a single numeral sign-class's multiplier
    values, for the Dehaene & Mehler monotonic-decrease cross-linguistic check.
    Caller supplies which numeral_class to test (e.g. Indus sign '002', Proto-
    Elamite 'N01') since sign-identification itself is corpus-specific.
    """
    counts = Counter()
    for d in corpus.documents:
        for val, ncls in d.numerals:
            if ncls == numeral_class:
                counts[val] += 1
    total = sum(counts.values())
    dist = {v: counts.get(v, 0) / total if total else float('nan') for v in value_range}
    is_monotonic = all(dist[v] >= dist[v + 1] - 1e-9 for v in list(value_range)[:-1])
    vals_present = [v for v in value_range if counts.get(v, 0) > 0]
    exponent, r2 = (stats.power_law_fit([counts[v] for v in vals_present]) if len(vals_present) >= 2
                    else (float('nan'), float('nan')))
    return {'counts': dict(counts), 'share_by_value': dist, 'monotonic_decrease': is_monotonic,
            'power_law_exponent': exponent, 'power_law_r2': r2}


def p9_periodicity(corpus, classes: dict, lags=(1, 2, 3), min_segment_len=3, n_perm=2000):
    label_seqs = []
    for d in corpus.documents:
        for seg in d.segments:
            if len(seg) < min_segment_len:
                continue
            labs = [classes.get(s) for s in seg]
            if any(l is None for l in labs):
                continue
            label_seqs.append(labs)
    return {'n_sequences': len(label_seqs),
            'lag_results': stats.permutation_lag_test(label_seqs, list(lags), n_perm=n_perm)}


def p11_conditional_entropy(corpus, max_order=6, seed=13, min_length=1):
    """P11: Rao et al. (2009)-style conditional entropy -- does entropy drop
    as more context (preceding signs) is added, the way it does in natural
    language, rather than staying flat (a maximally random sequence) or
    collapsing to near zero (a maximally rigid one)?

    Default max_order=6 matches Rao et al.'s (2010) reply to Sproat (2010):
    Sproat's original counterexample was a synthetic system (Zipfian
    frequencies, conditional independence beyond order 1) that matched
    Indus's order-0/1 entropy band, so a test that stops at order 2 -- as
    this test's own first version did -- cannot distinguish real script
    structure from Sproat's construction. Rao et al.'s reply argues real
    block entropy keeps SCALING past order 1-2 in a way Sproat's
    conditionally-independent construction does not; gap_by_order (below)
    is this package's way of checking that claim directly, corpus by
    corpus, against an explicit implementation of Sproat's counterexample
    (see examples/run_sproat_counterexample.py) rather than only citing it.

    Reports entropy at each order 0..max_order (see
    stats.conditional_entropy_by_order) against two nulls of increasing
    strength:
      - within_segment_shuffle_order1: signs reordered at random within
        their own segment. This keeps each segment's actual sign multiset
        intact, so it still carries real co-occurrence information and
        only randomizes order -- a weak null that can sit well below
        order-0 entropy on a corpus of short segments (a word's own small
        sign set constrains it regardless of order), and should be read as
        testing "does order matter beyond co-occurrence," not "is there
        structure at all."
      - iid_resample_by_order: every segment's signs redrawn i.i.d. from
        the corpus's own order-0 marginal (stats.iid_resample_entropy), at
        every order 0..max_order, keeping segment lengths but discarding
        everything else. Compare the real corpus's order-k entropy against
        THIS null at the SAME order, not against order-0: the plug-in
        conditional-entropy estimator used here is itself downward-biased
        at higher orders whenever the sign inventory is large relative to
        the corpus's occurrence count (most (context, next-sign) cells are
        seen only a handful of times, and a context seen once has zero
        empirical entropy by construction regardless of true randomness,
        an effect that gets sharply worse as order grows and segments run
        out of length to supply contexts at all) -- an i.i.d. resample
        shares the real corpus's occurrence count and segment-length
        distribution and so reproduces that same bias at every order,
        which is why it can sit well below order-0 entropy too. A real gap
        between the corpus's own order-k entropy and this null's AT ORDER
        k (gap_by_order[k]), not between the corpus and order-0, is the
        genuine-structure signal at that order (see
        stats.iid_resample_entropy's docstring). At high orders on
        short-segment corpora this gap can become uninformative rather
        than small: once almost every order-k context is unique in both
        the real corpus and its own null, both plug-in entropies collapse
        toward zero together and the gap says nothing either way -- report
        n_contexts/n_observations alongside the gap so a reader can tell
        a genuine null result from a data-starved one.
    Neither null failing is, on its own, evidence a corpus encodes language
    (see this package's README section on Raghavendra (2026) for why a
    purpose-built non-linguistic system can pass this and stronger checks
    too).

    min_length restricts to segments of at least this length before
    anything else runs (real, shuffle, and i.i.d.-null alike). Use it to
    length-match corpora before comparing their higher-order gaps: a
    corpus whose segments are mostly length 1-2 (e.g. Linear A, Proto-
    Elamite; median segment length 1 sign, see P7) contributes almost no
    order-2 observations at min_length=1, so its order-2 gap there is
    measured on a small, length-biased subset of its own data rather than
    a fair comparison against a corpus with longer segments throughout.
    """
    segments = [seg for d in corpus.documents for seg in d.segments if seg and len(seg) >= min_length]
    by_order = stats.conditional_entropy_by_order(segments, max_order=max_order)
    random.seed(seed)
    shuffled = []
    for seg in segments:
        s = seg[:]
        random.shuffle(s)
        shuffled.append(s)
    within_segment_shuffle_order1 = stats.conditional_entropy_by_order(shuffled, max_order=1)[1]
    order0_freq = defaultdict(int)
    for seg in segments:
        for s in seg:
            order0_freq[s] += 1
    iid_resample_by_order = stats.iid_resample_entropy(segments, dict(order0_freq), max_order=max_order, seed=seed)
    gap_by_order = {
        k: by_order[k]['entropy'] - iid_resample_by_order[k]['entropy']
        for k in range(1, max_order + 1)
        if by_order[k]['n_observations'] and iid_resample_by_order[k]['n_observations']
    }
    drop = (by_order[0]['entropy'] - by_order[1]['entropy']
            if max_order >= 1 and by_order[1]['n_observations'] else float('nan'))
    return {'by_order': by_order, 'within_segment_shuffle_order1': within_segment_shuffle_order1,
            'iid_resample_by_order': iid_resample_by_order, 'gap_by_order': gap_by_order,
            'entropy_drop_0_to_1': drop}


def p11_bootstrap_ci(corpus, orders=(1, 2), n_boot=200, seed=13, min_length=1):
    """Bootstrap check for P11's real-vs-i.i.d.-null entropy gap at each
    requested order (see p11_conditional_entropy and
    stats.bootstrap_entropy_gap's docstrings). Corpus sizes in this
    package range from a few thousand occurrences to several million, and
    the pure-Python plug-in estimator this package uses throughout has no
    fast path -- n_boot is the caller's to scale down for large corpora
    (see each run_*.py script's own choice, made and reported explicitly
    rather than silently thinned) rather than something this function
    picks for itself.

    IMPORTANT CAVEAT, found while building this: naive percentile
    bootstrap on this gap is itself asymmetrically biased, and the bias is
    not small. Resampling segments WITH REPLACEMENT necessarily creates
    exact-duplicate segments that a genuinely new same-size real corpus
    would not contain; the real (structured) side's conditional-entropy
    estimate is far more sensitive to duplicated, already-repeated
    sequences than the i.i.d. null's fresh-drawn side is -- confirmed
    empirically on Linear A order 1, where bootstrap resampling deflates
    the real side's mean entropy by ~0.36 bits but the null side's by only
    ~0.06. The bootstrap gap distribution therefore sits systematically
    below (more negative than) the true point estimate, and ci_lo/ci_hi
    below are reported for transparency only -- they should NOT be read
    as a classical 95% interval around point_gap; they will often not even
    contain it. sign_stability is the number this function's callers
    should actually use: the fraction of bootstrap replicates whose gap
    has the same sign as the real corpus's own point_gap. It sidesteps the
    duplication artifact entirely and answers the question that actually
    matters for a claim like "corpus X's order-2 gap is positive" -- is
    the SIGN robust to resampling, not the exact magnitude.

    min_length restricts to segments of at least this length before
    resampling -- see p11_conditional_entropy's docstring for why this
    matters for comparing corpora whose segments are mostly too short to
    supply a higher-order context at all (Linear A, Proto-Elamite).
    """
    segments = [seg for d in corpus.documents for seg in d.segments if seg and len(seg) >= min_length]
    order0_freq = defaultdict(int)
    for seg in segments:
        for s in seg:
            order0_freq[s] += 1
    freq = dict(order0_freq)
    real_by_order = stats.conditional_entropy_by_order(segments, max_order=max(orders))
    null_by_order = stats.iid_resample_entropy(segments, freq, max_order=max(orders), seed=seed)
    out = {}
    for k in orders:
        point_gap = real_by_order[k]['entropy'] - null_by_order[k]['entropy']
        gaps = stats.bootstrap_entropy_gap(segments, freq, k, n_boot=n_boot, seed=seed)
        lo_idx = max(0, int(round(0.025 * (len(gaps) - 1))))
        hi_idx = min(len(gaps) - 1, int(round(0.975 * (len(gaps) - 1))))
        if point_gap > 0:
            same_sign = sum(1 for g in gaps if g > 0)
        elif point_gap < 0:
            same_sign = sum(1 for g in gaps if g < 0)
        else:
            same_sign = sum(1 for g in gaps if g == 0)
        out[k] = {'point_gap': point_gap, 'sign_stability': same_sign / len(gaps) if gaps else float('nan'),
                   'ci_lo': gaps[lo_idx], 'ci_hi': gaps[hi_idx], 'median': gaps[len(gaps) // 2],
                   'n_boot': len(gaps)}
    return out


def p10_totaling_tablet_test(corpus, n_perm_or_shuffle=True, seed=1, marker_predicate=None):
    """Englund-style external-validation helper: does a text's final numeral
    line/segment total equal the sum of the numerals in the preceding ones?
    Needs Document.segment_numerals populated per-document, in the same order
    as Document.segments (see types.Document's docstring). For scripts that
    record numerals per-document rather than per-segment, group before
    calling. Returns the observed hit rate and a shuffled-order null for
    comparison -- the actual "does this beat chance" call is the caller's,
    since what counts as a fair null differs by corpus (see the Indus vs.
    Proto-Elamite writeups, where the same test was underpowered on one and
    informative on the other).

    marker_predicate: optional callable(signs: list[str]) -> bool, restricting
    the test to documents whose LAST numeral-bearing segment's sign list
    satisfies it -- e.g. `lambda signs: 'KU-RO' in signs` to test only texts
    that literally close with a word meaning "total", instead of testing
    every multi-numeral-line text regardless of whether its last line is
    actually a summary line. This is the genre precondition Englund's method
    needs (an itemized list *followed by a stated total*, not just any run of
    numeral-bearing lines) -- without it, the test conflates real ledgers with
    ordinary multi-line inventories that never state a total at all, diluting
    any real signal toward the null. Discovered by running this refinement on
    Linear A (unrestricted: 3.0% vs. 1.5% null on 198 texts; restricted to
    texts ending in a literal "total" word: 36.4% on 11 texts, three of the
    seven misses independently confirmed as scribal errors by the source
    corpus's own metadata) and on Proto-Elamite (restricting by the mere
    presence of a plausible totalizer sign, rather than a word whose sole
    function is "total", barely moved the rate: 33.3% vs. 22.8% null,
    close to the unrestricted 17.7% vs. 13.1% baseline) -- the precondition
    matters, but only when the marker is a dedicated summation word, not just
    a sign that sometimes plays that role.
    """
    import random
    random.seed(seed)

    def totals_by_class(numerals):
        d = defaultdict(int)
        for val, ncls in numerals:
            d[ncls] += val
        return d

    def numeral_segments(doc, order=None):
        pairs = list(zip(doc.segments, doc.segment_numerals))
        if order is not None:
            pairs = [pairs[i] for i in order]
        return [(signs, nums) for signs, nums in pairs if nums]

    def check(pairs):
        *body, (last_signs, last_nums) = pairs
        if marker_predicate is not None and not marker_predicate(last_signs):
            return None
        last_totals = totals_by_class(last_nums)
        body_totals = defaultdict(int)
        for _, seg_nums in body:
            for ncls, v in totals_by_class(seg_nums).items():
                body_totals[ncls] += v
        return any(body_totals.get(ncls) == tot and tot > 0 for ncls, tot in last_totals.items())

    hits, tested = 0, 0
    for d in corpus.documents:
        pairs = numeral_segments(d)
        if len(pairs) < 3:
            continue
        result = check(pairs)
        if result is None:
            continue
        tested += 1
        if result:
            hits += 1

    null_hits, null_tested = 0, 0
    if n_perm_or_shuffle:
        for d in corpus.documents:
            pairs = numeral_segments(d)
            if len(pairs) < 3:
                continue
            order = list(range(len(pairs)))
            random.shuffle(order)
            shuffled = [pairs[i] for i in order]
            result = check(shuffled)
            if result is None:
                continue
            null_tested += 1
            if result:
                null_hits += 1

    return {'tested': tested, 'hits': hits, 'hit_rate': hits / tested if tested else float('nan'),
            'null_hits': null_hits, 'null_tested': null_tested,
            'null_rate': null_hits / null_tested if null_tested else float('nan')}
