"""Pure-Python statistical primitives shared across the P1-P10 tests.

No numpy/scipy/networkx dependency, by original necessity (none were available
in the environment this protocol was first built in) and now by design -- it
keeps the toolkit runnable anywhere with a bare Python 3 interpreter.
"""
import math
import random
from collections import defaultdict


def power_law_fit(values_desc: list[float], cutoff: int | None = None):
    """OLS fit of log(value) against log(rank) for a descending-sorted list.

    Used for both the raw Zipf rank-frequency fit (values = sign frequencies)
    and Fuls' Modified Power Law (values = counts of how many signs share each
    frequency value). Returns (exponent, r_squared) with exponent reported as
    a positive number (i.e. the fit is value ~ rank^-exponent).
    """
    n = len(values_desc) if cutoff is None else min(cutoff, len(values_desc))
    xs = [math.log(i + 1) for i in range(n)]
    ys = [math.log(v) for v in values_desc[:n]]
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    slope = num / den
    intercept = my - slope * mx
    ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - my) ** 2 for y in ys)
    r2 = 1 - ss_res / ss_tot if ss_tot else float('nan')
    return -slope, r2


def modified_power_law(sign_freq: dict):
    """Fuls-style MPL: rank by how many signs share each frequency value."""
    freq_of_freq = defaultdict(int)
    for f in sign_freq.values():
        freq_of_freq[f] += 1
    ranked = sorted(freq_of_freq.values(), reverse=True)
    return power_law_fit(ranked)


def contingency_residuals(table: dict[tuple, int]):
    """Standardized and margin-adjusted residuals for a site x type table.

    table: {(row_key, col_key): observed_count}
    Returns {(row_key, col_key): {'observed', 'expected', 'std_resid', 'adj_resid'}}.
    adj_resid is the Haberman margin-adjusted standardized residual:
        (obs-exp) / sqrt(exp * (1 - row_total/n) * (1 - col_total/n))
    which is what should be compared to a normal-distribution critical value;
    the plain std_resid is what most "corrected chi" style claims in the wild
    actually report, and the two are worth publishing side by side since they
    can differ by 30-50%.
    """
    row_tot, col_tot, n = defaultdict(int), defaultdict(int), 0
    for (r, c), obs in table.items():
        row_tot[r] += obs
        col_tot[c] += obs
        n += obs
    out = {}
    for (r, c), obs in table.items():
        exp = row_tot[r] * col_tot[c] / n
        std_resid = (obs - exp) / math.sqrt(exp) if exp > 0 else float('nan')
        margin = (1 - row_tot[r] / n) * (1 - col_tot[c] / n)
        adj_resid = (obs - exp) / math.sqrt(exp * margin) if exp > 0 and margin > 0 else float('nan')
        out[(r, c)] = {'observed': obs, 'expected': exp, 'std_resid': std_resid, 'adj_resid': adj_resid}
    return out


def cramers_v(table: dict[tuple, int]):
    """Cramer's V, the chi-square-derived association-strength measure in
    [0, 1] (0 = no association, 1 = perfect association) between the row
    and column categories of an r x c contingency table of raw counts.

    table: {(row_key, col_key): observed_count} -- every (row, col)
    combination should have an explicit entry (including 0) so row/column
    counts are correct even when a cell is empty.

    Used here for Nair (2026)/Farmer-Sproat-Witzel (2004)'s "positional
    rigidity" metric (see tests.fsw_positional_rigidity): a 2x3 table of
    {one sign, all other signs} x {initial, medial, final position}. No
    small-sample (bias) correction is applied -- this is the same
    asymptotic V the FSW literature itself reports, not a corrected
    variant, so this package's own numbers are comparable to theirs.
    """
    row_tot, col_tot, n = defaultdict(int), defaultdict(int), 0
    for (r, c), obs in table.items():
        row_tot[r] += obs
        col_tot[c] += obs
        n += obs
    if n == 0:
        return float('nan')
    chi2 = 0.0
    for (r, c), obs in table.items():
        exp = row_tot[r] * col_tot[c] / n
        if exp > 0:
            chi2 += (obs - exp) ** 2 / exp
    denom = n * (min(len(row_tot), len(col_tot)) - 1)
    return math.sqrt(chi2 / denom) if denom > 0 else float('nan')


def diversity_ok(table: dict[tuple, int], min_rows=2, min_cols=2, max_dominant_share=0.85):
    """P3 precondition check: is there enough site/type variance to test at all?

    Returns (ok: bool, reason: str). Call this before contingency_residuals --
    running the test anyway on a near-monoculture table (one dominant site or
    type) produces numbers with no interpretive content (see the Proto-Elamite
    dossier's P3 section, where 94.7% of the corpus was a single site).
    """
    rows = {r for r, c in table}
    cols = {c for r, c in table}
    if len(rows) < min_rows or len(cols) < min_cols:
        return False, f"only {len(rows)} row categories / {len(cols)} col categories present"
    n = sum(table.values())
    row_tot = defaultdict(int)
    for (r, c), obs in table.items():
        row_tot[r] += obs
    dominant = max(row_tot.values()) / n
    if dominant > max_dominant_share:
        return False, f"one row category holds {dominant*100:.1f}% of the corpus (>{max_dominant_share*100:.0f}% threshold)"
    return True, "ok"


def cnm_modularity(nodes: list, edge_weight: dict[tuple[int, int], float]):
    """Greedy Clauset-Newman-Moore modularity maximization.

    edge_weight: {(i, j): w} with i < j, indices into `nodes`.
    Returns (communities: dict[repr_idx -> list[node_idx]], Q: float) at the
    highest-modularity point reached during the greedy merge sequence (not
    necessarily the fully-merged end state -- Q can rise then fall).
    """
    N = len(nodes)
    deg = [0.0] * N
    for (i, j), w in edge_weight.items():
        deg[i] += w
        deg[j] += w
    m2 = sum(deg)
    if m2 == 0:
        return {i: [i] for i in range(N)}, 0.0

    comm = list(range(N))
    comm_deg = deg[:]
    members = {i: [i] for i in range(N)}
    Q_running = -sum((d / m2) ** 2 for d in deg)
    best_Q, best_partition = Q_running, {i: [i] for i in range(N)}

    improved = True
    while improved:
        improved = False
        pair_weight = defaultdict(float)
        for (i, j), w in edge_weight.items():
            ci, cj = comm[i], comm[j]
            if ci == cj:
                continue
            key = (ci, cj) if ci < cj else (cj, ci)
            pair_weight[key] += w
        best_dq, best_pair = 0.0, None
        for (ci, cj), w_ij in pair_weight.items():
            dq = 2 * (w_ij / m2 - (comm_deg[ci] * comm_deg[cj]) / (m2 * m2))
            if dq > best_dq:
                best_dq, best_pair = dq, (ci, cj)
        if best_pair and best_dq > 1e-12:
            ci, cj = best_pair
            for n in members[cj]:
                comm[n] = ci
            members[ci].extend(members[cj])
            members[cj] = []
            comm_deg[ci] += comm_deg[cj]
            comm_deg[cj] = 0
            Q_running += best_dq
            improved = True
            if Q_running > best_Q:
                best_Q = Q_running
                best_partition = {c: list(mem) for c, mem in members.items() if mem}
    return {c: mem for c, mem in best_partition.items() if mem}, best_Q


def shannon_entropy(freq: dict) -> float:
    """Shannon entropy, in bits, of a frequency distribution."""
    total = sum(freq.values())
    if total == 0:
        return float('nan')
    h = 0.0
    for c in freq.values():
        if c == 0:
            continue
        p = c / total
        h -= p * math.log2(p)
    return h


def conditional_entropy_by_order(segments: list[list[str]], max_order: int = 2):
    """Rao et al. (2009)-style conditional entropy H(X_n | X_{n-k}..X_{n-1})
    at orders 0..max_order, over sign sequences within segments (protocol P11).

    Order 0 is plain sign-frequency entropy; order k>=1 is the entropy of a
    sign given the k signs immediately before it, estimated directly from
    observed (context, next-sign) counts -- the same maximum-likelihood
    estimate Rao et al. use, with no smoothing. Contexts never cross a
    segment boundary, the same convention this package's other
    position-aware tests (P5, P7, P9) use. Each order's result reports
    n_contexts and n_observations alongside the entropy value, since a
    higher-order estimate on a small corpus can be thin enough that its
    number is more a measurement of data sparsity than of the script.
    """
    order0_freq = defaultdict(int)
    for seg in segments:
        for s in seg:
            order0_freq[s] += 1
    results = {0: {'entropy': shannon_entropy(order0_freq), 'n_contexts': 1 if order0_freq else 0,
                   'n_observations': sum(order0_freq.values())}}
    for k in range(1, max_order + 1):
        context_next = defaultdict(lambda: defaultdict(int))
        for seg in segments:
            for i in range(k, len(seg)):
                context_next[tuple(seg[i - k:i])][seg[i]] += 1
        total_obs = sum(sum(d.values()) for d in context_next.values())
        if total_obs == 0:
            results[k] = {'entropy': float('nan'), 'n_contexts': 0, 'n_observations': 0}
            continue
        h = 0.0
        for nxts in context_next.values():
            p_ctx = sum(nxts.values()) / total_obs
            h += p_ctx * shannon_entropy(nxts)
        results[k] = {'entropy': h, 'n_contexts': len(context_next), 'n_observations': total_obs}
    return results


def iid_resample_entropy(segments: list[list[str]], freq: dict, max_order: int = 1, seed: int = 13):
    """Redraws every segment's signs i.i.d. from `freq` (typically the
    corpus's own order-0 marginal), keeping each segment's length but
    discarding everything else -- the correct null for isolating real
    sequential/co-occurrence structure from two things that inflate a raw
    order0-vs-order1 entropy drop even when none is present.

    The first is what shuffling signs within their own segment controls
    for: an in-segment shuffle keeps each segment's actual sign multiset
    intact (so it still carries real co-occurrence information) and only
    randomizes order, so it isolates "does order matter beyond
    co-occurrence" -- a weaker null than this one.

    The second, which this function exists to isolate, is the plug-in
    conditional-entropy estimator's own small-sample bias: with an
    inventory of hundreds of signs, the number of possible (context,
    next-sign) cells is large relative to any real corpus's occurrence
    count, so most contexts are seen only a handful of times -- a context
    seen once has, by construction, zero empirical conditional entropy
    regardless of the script's true randomness. This estimator is
    downward-biased on ANY sequence at this order and sample size,
    including a genuinely i.i.d. one, which is exactly why an i.i.d.
    resample sharing the real corpus's occurrence count and segment-length
    distribution is the correct baseline: it reproduces that same bias, so
    it does NOT come out close to order-0 entropy in general -- the real
    corpus's order-1 entropy should be compared against THIS null, not
    against order-0, and a real gap between the two (rather than between
    real and order-0) is the genuine-structure signal.
    """
    signs = list(freq.keys())
    weights = list(freq.values())
    random.seed(seed)
    resampled = [random.choices(signs, weights=weights, k=len(seg)) for seg in segments if seg]
    return conditional_entropy_by_order(resampled, max_order=max_order)


def bootstrap_entropy_gap(segments: list[list[str]], freq: dict, order: int, n_boot: int = 200, seed: int = 13):
    """Nonparametric bootstrap for P11's real-vs-i.i.d.-null conditional-
    entropy gap at a single order (see tests.p11_conditional_entropy's
    docstring for what the gap means and why the i.i.d. null, not order-0
    entropy, is the right comparison in the first place).

    Each of n_boot replicates resamples len(segments) segments WITH
    REPLACEMENT from the real corpus, computes that resample's own
    order-`order` conditional entropy, and compares it against a fresh
    i.i.d. resample drawn from the same fixed marginal `freq` but sharing
    THAT replicate's own resampled segment-length multiset -- the same
    pairing principle iid_resample_entropy itself uses, applied per
    replicate so the null's finite-sample bias tracks the real resample's
    own sparsity at every draw, not just the original corpus's.

    Returns the n_boot gap values, sorted ascending. Percentiles of this
    list (e.g. indices at 2.5%/97.5% for a 95% CI) are the caller's to pull
    out, since which interval matters depends on the specific claim being
    tested. A local random.Random instance is used for the outer resampling
    draw so this function's own randomness doesn't get clobbered by
    iid_resample_entropy's internal global re-seeding on each call.
    """
    rng = random.Random(seed)
    n = len(segments)
    gaps = []
    for _ in range(n_boot):
        resampled = [segments[rng.randrange(n)] for _ in range(n)]
        real_h = conditional_entropy_by_order(resampled, max_order=order)[order]['entropy']
        null_h = iid_resample_entropy(resampled, freq, max_order=order,
                                       seed=rng.randrange(1_000_000))[order]['entropy']
        gaps.append(real_h - null_h)
    gaps.sort()
    return gaps


def permutation_lag_test(label_seqs: list[list], lags: list[int], n_perm: int = 2000, seed: int = 42):
    """Permutation test for same-label recurrence at each lag (protocol P9).

    label_seqs: one list of labels per segment (already filtered to segments
        with every position labeled -- see decipherment_protocol.tests.p9).
    For each lag k, shuffles labels *within* each segment n_perm times to build
    a null distribution, then reports the observed same-label rate against it.

    Known limitation, carried over honestly from every dossier that has used
    this test so far: if the labels themselves were derived from position
    (as P5's initial/medial/final classes are), the lag-1 result is partly
    circular by construction. Do not treat this test's lag-1 output as
    independent evidence without a position-free label source.
    """
    random.seed(seed)

    def rate(seqs, lag):
        same, total = 0, 0
        for labs in seqs:
            n = len(labs)
            for i in range(n - lag):
                total += 1
                if labs[i] == labs[i + lag]:
                    same += 1
        return (same / total, total) if total else (float('nan'), 0)

    results = {}
    for lag in lags:
        obs_rate, obs_total = rate(label_seqs, lag)
        if obs_total == 0:
            results[lag] = None
            continue
        perm_rates = []
        for _ in range(n_perm):
            shuffled = []
            for labs in label_seqs:
                s = labs[:]
                random.shuffle(s)
                shuffled.append(s)
            r, _ = rate(shuffled, lag)
            perm_rates.append(r)
        mean_p = sum(perm_rates) / n_perm
        sd_p = math.sqrt(sum((r - mean_p) ** 2 for r in perm_rates) / (n_perm - 1))
        z = (obs_rate - mean_p) / sd_p if sd_p > 0 else float('nan')
        p = math.erfc(abs(z) / math.sqrt(2))
        results[lag] = {'observed_rate': obs_rate, 'n': obs_total, 'null_mean': mean_p,
                         'null_sd': sd_p, 'z': z, 'p': p}
    return results
