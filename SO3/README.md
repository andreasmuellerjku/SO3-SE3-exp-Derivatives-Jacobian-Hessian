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

# SO3Exp

Wolfram Language package for fast, numerically stable computations on the rotation group SO(3) using canonical coordinates (rotation vectors). Provides exponential/logarithm maps and the differential of the exponential (dexp), its inverse, and first/second Fréchet derivatives, implemented via sinc-based closed forms.

- Author: Andreas Mueller
- Created: 2 Jan 2005 • Updated: 1 Apr 2026
- Depends on: SO3Core`

## Features

- `SO3Exp[x]`: rotation matrix $R = \exp(\tilde{x})$ using a branch-free, sinc-based Rodrigues form.
- `SO3Log[R, eps]`: rotation vector $x$ with $R=\exp(\tilde{x})$; small-angle threshold `eps` (default `SO3EPS`).
- `SO3dexp[x, eps]`: right-trivialized differential of dexp at $x$.
- `SO3dexpInv[x, eps]`: inverse $(\mathrm{dexp}_x)^{-1}$.
- `SO3Ddexp[x, y, eps]`: Fréchet derivative $D_x(\mathrm{dexp})(y)$ of dexp at $x$ along $y$.
- `SO3D2dexp[x, y, u, v, eps]`: D derivative $D_X(\mathrm{Ddexp})(U)$ of $\mathrm{Ddexp}(X):=(D_x\mathrm{dexp})(y)$ at $X=(x,y)$.
- `SO3DdexpInv[x, y, eps]`: derivative of $\mathrm{dexp}^{-1}$ at $x$.
- `SO3D2dexpInv[x, y, u, v, eps]`: second derivative of $\mathrm{dexp}^{-1}$ at $x$.

All functions accept symbolic or numeric inputs; small-angle handling governed by `eps` (default `SO3EPS = 1.*^-12`).

## Install

Place `SO3Core.m` and `SO3Exp.m` on your `$Path` (e.g., `~/Library/Mathematica/Applications/`, `~/.Mathematica/Applications/`, or `%AppData%\Mathematica\Applications\`), then:

```wl
<< SO3Core`
<< SO3Exp`
```

## Usage

```wl
<< SO3Core`
<< SO3Exp`

x = {0.1, -0.2, 0.05};

(* Exponential and logarithm *)
R = SO3Exp[x];
SO3Log[R]                (* ≈ x *)

(* dexp and its inverse *)
J  = SO3dexp[x];
Jinv = SO3dexpInv[x];
Chop[J . Jinv - IdentityMatrix[3], SO3EPS]  (* 0 *)

(* First derivative of dexp along y *)
y = {0.3, 0.2, -0.1};
D1 = SO3Ddexp[x, y];

(* Second derivative examples *)
u = {0.1, 0.0, 0.2};
v = {-0.05, 0.02, 0.03};
D2  = SO3D2dexp[x, y, u, v];
D2i = SO3D2dexpInv[x, y, u, v];
```

Mathematical form (rotation vector $x$, norm $n=\|x\|$):
$$
\exp(\tilde{x}) = I + \mathrm{sinc}(n)\,\tilde{x} + \tfrac12\,\mathrm{sinc}^2\!\left(\tfrac{n}{2}\right)\tilde{x}^2,
\quad \mathrm{sinc}(t)=\frac{\sin t}{t}.
$$

## Notes

- Inputs are not shape‑validated; ensure vectors are length‑3 and matrices are valid rotations when using `SO3Log`.
- For very small angles, functions use `eps` to avoid division by tiny numbers; adjust if needed for your tolerance.

## References

- A. Mueller, Review of the exponential and Cayley map on SE(3) as relevant for Lie group integration..., Proc. Royal Soc. A, 2021. doi: 10.1098/rspa.2021.0303
- Preprint: https://arxiv.org/abs/2303.07928
