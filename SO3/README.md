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
- Created: 2 Jan 2005
- Updated: 1 Apr 2026
- Depends on: SO3Core`

## Features

- `SO3Exp[x]`: rotation matrix $R = \exp(\tilde{x})$ using a branch-free, sinc-based Rodrigues form.
- `SO3Log[R, eps]`: rotation vector $x$ with $R=\exp(\tilde{x})$; small-angle threshold `eps` (default `SO3EPS`).
- `SO3dexp[x, eps]`: right-trivialized differential of dexp at $x$.
- `SO3dexpInv[x, eps]`: right-trivialized differential of $(\mathrm{dexp}_x)^{-1}$ at $x$.
- `SO3Ddexp[x, y, eps]`: Fréchet derivative $D_x(\mathrm{dexp})(y)$ of dexp at $x$ along $y$.
- `SO3DdexpInv[x, y, eps]`: Fréchet derivative $D_x(\mathrm{dexp}^{-1})(y)$ of dexp at $x$ along $y$.
- `SO3D2dexp[x, y, u, v, eps]`: derivative $D_X(\mathrm{Ddexp})(U)$ of $\mathrm{Ddexp}(X):=(D_x\mathrm{dexp})(y)$ at $X=(x,y)$ along $U=(u,v)$.
- `SO3D2dexpInv[x, y, u, v, eps]`: derivative $D_X(\mathrm{Ddexp}^{-1})(U)$ of $\mathrm{Ddexp^{-1}}(X):=(D_x\mathrm{dexp}^{-1})(y)$ at $X=(x,y)$ along $U=(u,v)$.

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
x == SO3Log[R]

(* dexp and its inverse *)
J  = SO3dexp[x];
Jinv = SO3dexpInv[x];
Chop[J . Jinv - IdentityMatrix[3], SO3EPS]  (* 0 *)

(* First derivative of dexp along y *)
y = {0.3, 0.2, -0.1};
D1 = SO3Ddexp[x, y];
D1inv = SO3DdexpInv[x, y];

(* Second derivatives *)
u = {0.1, 0.0, 0.2};
v = {-0.05, 0.02, 0.03};
D2  = SO3D2dexp[x, y, u, v];
D2inv = SO3D2dexpInv[x, y, u, v];
```

## Notes

- Inputs are not shape‑validated; ensure vectors are length‑3 and matrices are valid rotations when using `SO3Log`.
- For very small angles, functions use `eps` to avoid division by tiny numbers; adjust if needed for your tolerance.

## References

- A. Mueller: Review of the exponential and Cayley map on SE(3) as relevant for Lie group integration of the generalized Poisson equation and flexible multibody systems,  Proc. Royal Soc. A, September 2021; 477 (2253): 20210303.
https://doi.org/10.1098/rspa.2021.0303
Corrected preprint: https://arxiv.org/abs/2303.07928
