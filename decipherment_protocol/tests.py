"""P1-P10: the decipherment-statistics protocol itself, as functions over a Corpus.

Each function takes a decipherment_protocol.types.Corpus and returns a plain
dict of results -- no side effects, no printing, so callers can format results
however a given dossier needs. P4 and P10 are deliberately thin: what counts
as "allograph evidence" or "external validation" is corpus-specific by nature
(shape geometry vs. published referent tables; a whole monograph vs. one
glossed word), so those two only provide the generic helper a caller needs,
not a one-shot answer.
"""
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


def p2_rank_frequency(corpus, cutoffs=(30, 60, 100, None)):
    freq = corpus.sign_freq()
    ranked = sorted(freq.values(), reverse=True)
    zipf = {}
    for c in cutoffs:
        label = c if c is not None else len(ranked)
        s, r2 = stats.power_law_fit(ranked, c)
        zipf[label] = {'exponent': s, 'r2': r2}
    g, r2_mpl = stats.modified_power_law(freq)
    return {'zipf_by_cutoff': zipf, 'mpl_exponent': -g, 'mpl_r2': r2_mpl}


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


def p6_cooccurrence_network(corpus, min_n=5):
    total_occ = corpus.sign_freq()
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


def p10_totaling_tablet_test(corpus, n_perm_or_shuffle=True, seed=1):
    """Englund-style external-validation helper: does a text's final numeral
    line/segment total equal the sum of the numerals in the preceding ones?
    Needs Document.numerals populated per-document; for scripts that record
    numerals per-segment rather than per-document, group before calling.
    Returns the observed hit rate and a shuffled-order null for comparison --
    the actual "does this beat chance" call is the caller's, since what counts
    as a fair null differs by corpus (see the Indus vs. Proto-Elamite writeups,
    where the same test was underpowered on one and informative on the other).
    """
    import random
    random.seed(seed)

    def totals_by_class(numerals):
        d = defaultdict(int)
        for val, ncls in numerals:
            d[ncls] += val
        return d

    hits, tested = 0, 0
    for d in corpus.documents:
        segs_with_nums = [seg_nums for seg_nums in d.segment_numerals if seg_nums]
        if len(segs_with_nums) < 3:
            continue
        tested += 1
        *body, last = segs_with_nums
        last_totals = totals_by_class(last)
        body_totals = defaultdict(int)
        for seg_nums in body:
            for ncls, v in totals_by_class(seg_nums).items():
                body_totals[ncls] += v
        if any(body_totals.get(ncls) == tot and tot > 0 for ncls, tot in last_totals.items()):
            hits += 1

    null_hits = 0
    if n_perm_or_shuffle:
        for d in corpus.documents:
            segs_with_nums = [seg_nums for seg_nums in d.segment_numerals if seg_nums]
            if len(segs_with_nums) < 3:
                continue
            random.shuffle(segs_with_nums)
            *body, last = segs_with_nums
            last_totals = totals_by_class(last)
            body_totals = defaultdict(int)
            for seg_nums in body:
                for ncls, v in totals_by_class(seg_nums).items():
                    body_totals[ncls] += v
            if any(body_totals.get(ncls) == tot and tot > 0 for ncls, tot in last_totals.items()):
                null_hits += 1

    return {'tested': tested, 'hits': hits, 'hit_rate': hits / tested if tested else float('nan'),
            'null_hits': null_hits, 'null_rate': null_hits / tested if tested else float('nan')}
