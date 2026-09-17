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
