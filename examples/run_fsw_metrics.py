"""Nair (2026)'s own four Farmer-Sproat-Witzel (2004) metrics -- text
brevity, formulaic repetition, hapax legomenon rate, and positional
rigidity -- run directly across this package's corpora, rather than only
cited (see README's "Nair (2026)" section for the comparison this extends).
Nair's own scorecard applies these four metrics to Indus alone, calibrated
against two purpose-built synthetic baselines; this script instead reuses
this package's own eight corpora (seven real, one synthetic-control), which
already carry an independent deciphered/undeciphered split P11 doesn't (see
tests.p11_conditional_entropy's docstring and the paper's P11 Results),
making a direct check possible: do these four metrics track that same
split, independently of P11's own entropy-based one?

Metric definitions (see tests.fsw_formulaic_repetition and
tests.fsw_positional_rigidity's docstrings for the two novel ones):
  - text brevity: mean inscription length in signs (this package's own
    per-document flat_signs length, i.e. p1_census's mean_text_length).
  - hapax legomenon rate: this package's own p1_census hapax_share.
  - formulaic repetition: at phrase lengths L=3..6, the share of distinct
    L-sign phrases recurring in two or more different inscriptions.
  - positional rigidity: mean Cramer's V (stats.cramers_v) of the
    corpus's ten most frequent signs' start/middle/end distribution
    against the corpus-wide positional marginal.

Each corpus is rebuilt by re-running the corpus-construction portion of its
own dedicated run_*.py script (via AST, stopping right after the line that
assigns its corpus variable) rather than reimplementing each script's own,
sometimes intricate, source parsing here a second time -- this guarantees
the exact same corpus every other test in this package runs against, with
one source of truth for how each corpus is built.
"""
import sys, os, ast
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from decipherment_protocol import tests

EXAMPLES = os.path.dirname(os.path.abspath(__file__)) + os.sep


def _prefix_through_assignment(src, varname):
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == varname for t in node.targets):
            end_line = node.end_lineno
            lines = src.splitlines(keepends=True)
            return ''.join(lines[:end_line])
    raise ValueError(f"no assignment to {varname!r} found in this script")


def _load_corpus(fname, varname):
    src = open(EXAMPLES + fname, encoding='utf-8').read()
    prefix = _prefix_through_assignment(src, varname)
    ns = {'__file__': EXAMPLES + fname, '__name__': '__main__'}
    exec(prefix, ns)
    return ns[varname]


# (script, corpus variable name in that script, decipherment status)
JOBS = [
    ('run_indus.py', 'corpus', 'undeciphered'),
    ('run_lineara.py', 'corpus', 'undeciphered'),
    ('run_protoelamite.py', 'corpus', 'undeciphered'),
    ('run_linearb.py', 'word_corpus', 'deciphered'),
    ('run_urIII.py', 'word_corpus', 'deciphered'),
    ('run_oldpersian.py', 'corpus', 'deciphered'),
    ('run_ugaritic.py', 'corpus', 'deciphered'),
    ('run_synthetic_control.py', 'corpus', 'non-linguistic (positive control)'),
    ('run_sproat_counterexample.py', 'corpus', 'non-linguistic (negative control)'),
]

if __name__ == '__main__':
    results = {}
    for fname, varname, status in JOBS:
        corpus = _load_corpus(fname, varname)
        print(f"\n=== {fname} ({status}) ===")

        p1 = tests.p1_census(corpus)
        lengths = [len(d.flat_signs) for d in corpus.documents if d.flat_signs]
        mean_text_len = sum(lengths) / len(lengths) if lengths else float('nan')
        print(f"  text brevity: mean inscription length = {mean_text_len:.2f} signs (n={len(lengths)} inscriptions)")
        print(f"  hapax legomenon rate: {p1['hapax_share']*100:.2f}% "
              f"({p1['hapax_legomena']}/{p1['inventory_size']} sign types)")

        fr = tests.fsw_formulaic_repetition(corpus)
        print("  formulaic repetition rate (share of distinct L-sign phrases recurring in >=2 inscriptions):")
        for L, r in fr.items():
            print(f"    L={L}: {r['repeated_types']}/{r['total_types']} = {r['repetition_rate']*100:.2f}%")

        pr = tests.fsw_positional_rigidity(corpus)
        print(f"  positional rigidity: mean Cramer's V (top-10 signs) = {pr['mean_v']:.3f}")

        results[fname] = {
            'status': status, 'mean_text_len': mean_text_len, 'hapax_rate': p1['hapax_share'],
            'formulaic': {L: r['repetition_rate'] for L, r in fr.items()}, 'rigidity': pr['mean_v'],
        }

    print("\n\n=== SUMMARY TABLE ===")
    print(f"{'Corpus':<28} {'Status':<32} {'Brevity':>8} {'Hapax%':>8} {'FR@3':>7} {'FR@6':>7} {'Rigidity':>9}")
    for fname, r in results.items():
        print(f"{fname:<28} {r['status']:<32} {r['mean_text_len']:>8.2f} {r['hapax_rate']*100:>7.2f}% "
              f"{r['formulaic'][3]*100:>6.2f}% {r['formulaic'][6]*100:>6.2f}% {r['rigidity']:>9.3f}")
