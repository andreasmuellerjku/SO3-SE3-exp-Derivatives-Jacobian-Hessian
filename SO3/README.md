# SO3Core

Minimal Wolfram Language utilities for working with the Lie algebra $\mathfrak{so}(3)$ and its identification with vectors in $\mathbb{R}^3$.

- Author: Andreas Mueller
- Created: 2 Jan 2005
- Updated: 1 Apr 2026

## Features

- `R3Toso3[x]`: map a 3‑vector $x$ to the skew‑symmetric matrix $\tilde{x}$ such that $\tilde{x}\,y = x \times y$.
- `so3ToR3[w]`: inverse of the above map from a $3\times 3$ skew‑symmetric matrix to its axis vector.
- `SO3EPS`: default numerical threshold (`1.*^-12`) for zero‑like values.

## Install

Place `SO3Core.m` on your `$Path` (e.g., `~/Library/Mathematica/Applications/`, `~/.Mathematica/Applications/`, or `%AppData%\Mathematica\Applications\`), then:

```wl
<< SO3Core`
```

## Usage

```wl
<< SO3Core`

x = {1, 2, 3};
X = R3Toso3[x]            (* {{0,-3,2},{3,0,-1},{-2,1,0}} *)
so3ToR3[X]                (* {1, 2, 3} *)

y = {4, -1, 5};
X.y == Cross[x, y]        (* True *)
```

## Notes

- Inputs are not validated; ensure vectors are length‑3 and matrices are skew‑symmetric.
- Works with numeric or symbolic data. Use `Chop[..., SO3EPS]` to suppress floating‑point noise.
