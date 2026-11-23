# Mathar's inverse-quintic recurrence

`note.pdf` proves the six-term recurrence for OEIS A392976 at every integer
index n >= 5 and constructs a nonzero homogeneous differential operator of
order five. It also proves the generating-function and binomial-sum formulas
for A392976 and the indicated formulas in A392977, A001002 and A370839.
A separate trace-zero argument proves that A392976 also has a homogeneous
differential equation of order at most four; an explicit operator of that
smaller order is not required for this existence proof.
The inverse family of `t((1-t^2)^r-t)` has integral coefficients and an
algebraic, D-finite quotient by x for every integer r. Its positive integral
powers have an explicit generalized-binomial formula; coefficients are
strictly positive for r >= 0.

The cubic specializations are reproved with attribution to Barry (2011),
Schuetz and Whieldon (2014), Birmajer, Gil and Weiner (2015), and Liszewska
and Mlotkowski (2019). Further written corollaries cover A001002's alternate
sum, recurrence, derivative representations, reciprocal convolutions,
logarithmic derivative and hypergeometric formula. The denominator inside
the arcsine in Barry's Example 7 is corrected from 6 to 16, as in OEIS.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively use a standard LaTeX installation:

```sh
pdflatex note.tex
pdflatex note.tex
```

The bibliography is included in `note.tex`; no external bibliography file is
required. The document uses `geometry`, `amsmath`, `amssymb`, `amsthm`,
`booktabs`, `array`, `fontenc`, `lmodern` and `hyperref`.

## Exact arithmetic

Python 3.8 or later with the standard library suffices:

```sh
python3 verify.py --bound 40 --output verification.json
```

The checker verifies a zero integer polynomial of degree at most 62 and the
nonzero leading coefficient of the order-five differential operator. It
compares recursive formal reversion and binomial sums for indices 0 through
40 in each of A001002, A392976, A392977 and A370839. The six-term recurrence
is also checked at indices 5 through 40. Its universal proof is the formal
change-of-variable argument in the note, supported by the polynomial identity.
It also verifies the cubic's inhomogeneous operator identity `L_B B = -x`,
which gives its recurrence at every n >= 2. Generalized-binomial formulas
and recursive reversion agree at n = 0..22 for integer r = -3..5 and powers
s = 1..4. Cubic reciprocal, logarithmic-derivative and hypergeometric checks
extend through n = 40. These finite comparisons supplement the written proofs.

The pinned declaration
`D5.S1.Recurrence.Algebraic.A392976InverseRecurrence.result` in
`A392976InverseRecurrence.lean` proves inverse compositions, inverse uniqueness,
the quotient-series identity, D-finiteness and the recurrence for every n >= 5.
Its exact scientific source is linked in the PDF and `metadata.json`.
The order-at-most-four bound, parameterized formulas and cubic corollaries
have written proofs in the note; the pinned `result` does not assert those
extensions or the initial five coefficient values.

## Attribution

Creators: Haobo Ma (ChronoAI Pte Ltd, Singapore), followed by Wenlin Zhang
(National University of Singapore; ChronoAI Pte Ltd, Singapore).
License for this note and its accompanying original source: CC BY 4.0.
The sequence is due to Seiichi Manyama; the recurrence conjecture is due to
R. J. Mathar. Published neighboring formulas retain their source attribution.
The linked OEIS entries are credited to the OEIS Foundation and contributors
and are available under CC BY-SA 4.0. The source archive does not reproduce
the OEIS entry exports or the Lean development.
