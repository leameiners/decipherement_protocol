"""Common input schema for the P1-P10 decipherment-statistics protocol.

Every script this protocol has been run on (Indus, Linear A, Proto-Elamite) has a
different native structure -- flat sign strings with no word breaks, a real
word-divider glyph, or numbered accounting lines with a separate numeral
sub-alphabet. A Document normalizes all three down to one shape so the same
P1-P10 functions run unmodified on any of them.
"""
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Document:
    """One inscription / tablet / text.

    segments: the text broken into its natural sub-units, outermost first.
        - No known word/line boundary (Indus): segments = [[all signs in the text]]
          -- one segment containing everything.
        - A real word-divider (Linear A): segments = [[word1 signs], [word2 signs], ...]
        - Numbered accounting lines (Proto-Elamite): segments = [[line1 signs], [line2 signs], ...]
      Every P-test that needs an internal-position concept (P5, P7, P9) operates
      *within* a segment; P1/P2/P4/P6 pool across segments.
    numerals: flat stream of (value, numeral_class) tuples across the whole
      document, for scripts whose number system uses a distinct sign set from
      the main inventory (Proto-Elamite's N-signs; Linear A's Aegean numeral
      block). Leave as [] for a script whose numerals are ordinary inventory
      signs (Indus signs 001-013), which just live in `segments` instead.
    segment_numerals: the same numeral stream, but split per-segment in the
      same order as `segments` -- one list of (value, class) tuples per
      segment. Only needed by p10_totaling_tablet_test, which checks whether a
      later segment's total equals the sum of the earlier ones; leave as []
      if that test isn't being run on this corpus.
    """
    doc_id: str
    segments: list[list[str]]
    site: Optional[str] = None
    artifact_type: Optional[str] = None
    numerals: list[tuple[int, str]] = field(default_factory=list)
    segment_numerals: list[list[tuple[int, str]]] = field(default_factory=list)

    @property
    def flat_signs(self) -> list[str]:
        return [s for seg in self.segments for s in seg]


@dataclass
class Corpus:
    name: str
    documents: list[Document]

    def sign_freq(self):
        from collections import Counter
        c = Counter()
        for d in self.documents:
            c.update(d.flat_signs)
        return c
