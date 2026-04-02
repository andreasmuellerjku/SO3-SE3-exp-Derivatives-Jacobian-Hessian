# SE3Core

Minimal Wolfram Language utilities for working with the Lie algebra se(3) and the Lie group $\mathrm{SE}(3)$ (rigid‑body motions). Provides screw vector to $4\times 4$ matrix maps, the adjoint (Lie bracket) in 6‑vector form, homogeneous transform assembly, and group inverse.

- Author: Andreas Mueller
- Created: 2 Jan 2005
- Updated: 1 Apr 2026
- Depends on: SO3Core`

## Features

- `TwistTose3[X]`: hat map for a twist $X=(w,v)\in\mathbb{R}^6$,
- `se3ToTwist[Xhat]`: map from a $4\times 4$ $\mathfrak{se}(3)$ matrix to the twist $(w,v)$.
- `SE3adMatrix[X]`: $6\times 6$ adjoint matrix $\mathrm{ad}_X$ such that $\mathrm{ad}_X Y = [X,Y]$,
- `RotPosToSE3[R,p]`: assemble homogeneous transformation matrix $H\in SE(3)$ from Rotation matrix R and position vector p.
- `SE3Inverse[H]`: inverse $H^{-1}$ of $H\in SE(3).$
- `SE3EPS`: default numerical threshold (`1.*^-12`) for zero‑like values.

## Install

Place `SO3Core.m` and `SE3Core.m` on your `$Path` (e.g., `~/Library/Mathematica/Applications/`, `~/.Mathematica/Applications/`, or `%AppData%\Mathematica\Applications\`), then:

```wl
<< SO3Core`
<< SE3Core`
```

## Usage

```wl
<< SO3Core`
<< SE3Core`

(* Twist <-> se(3) *)
X = {0.1, -0.2, 0.05,  1.0, 0.3, -0.4};   (* {w, v} *)
Xhat = TwistTose3[X];
X == se3ToTwist[Xhat]

(* Adjoint (Lie bracket in 6D) *)
A = SE3adMatrix[X];
Y = {0.2, 0.1, 0.0,  -0.3, 0.2, 0.1};
A.Y

(* Homogeneous transforms *)
R = SO3Exp[{0.1, 0.0, 0.2}];
p = {1., 2., 3.};
H = RotPosToSE3[R, p];
Hinverse = SE3Inverse[H];
Chop[Hinverse.H - IdentityMatrix[4], SE3EPS]  (* 0 *)
```

# SE3PExp

Wolfram Language package for numerically stable computations on the rigid‑body group SE(3) using canonical coordinates and the semi‑direct product block structure SE(3) = SO(3) ⋉ ℝ³. Implements exponential/logarithm maps for twists, the differential of the exponential (dexp), its inverse, and first derivatives, by leveraging SO(3) building blocks.

- Author: Andreas Mueller
- Created: 2 Jan 2005
- Updated: 1 Apr 2026
- Depends on: SO3Core, SO3Exp, SE3Core

## Features

- `SE3PExp[X, eps]`: homogeneous transform $H=\exp(\widehat{X})$ for screw coordinates $X=\{x,y\}\in\mathbb{R}^6$; small‑angle threshold `eps` (default `SE3EPS`).
- `SE3PLog[H, eps]`: screw coordinates $X$ with $H=\exp(\widehat{X})$; uses block structure with SO(3) log and dexp inverse.
- `SE3Pdexp[X, eps]`: right-trivialized differential dexp of $exp$ at $X$ computed using the $3\times 3$ block partitioning.
- `SE3PdexpInv[X, eps]`: inverse $\mathrm{dexp}^{-1}_X$ of $exp_X$.
- `SE3PDdexp[X, U, eps]`: Frechet derivative $D_X(\mathrm{dexp})(U)$ of $\mathrm{dexp}$ at $X$ along $U$.
- `SE3PDdexpInv[X, U, eps]`: Frechet derivative $D_X(\mathrm{dexp}^{-1})(U)$ of of $\mathrm{dexp}^{-1}$ at $X$.

## Install

Place `SO3Core.m`, `SO3Exp.m`, `SE3Core.m`, and `SE3PExp.m` on your `$Path` (e.g., `~/Library/Mathematica/Applications/`, `~/.Mathematica/Applications/`, or `%AppData%\Mathematica\Applications\`), then:

```wl
<< SO3Core`
<< SO3Exp`
<< SE3Core`
<< SE3PExp`
```

## Usage

```wl
<< SO3Core`
<< SO3Exp`
<< SE3Core`
<< SE3PExp`

X = {0.1, -0.2, 0.05,  1.0, 0.3, -0.4};  (* twist {x, y} *)

(* Exponential and logarithm *)
H = SE3PExp[X];
Xrec = se3ToTwist @ SE3PLog[H];          (* ≈ X *)

(* dexp and its inverse *)
J   = SE3Pdexp[X];
Jinv = SE3PdexpInv[X];
Chop[J . Jinv - IdentityMatrix[6], SE3EPS]  (* 0 *)

(* Derivatives *)
U = {0.02, 0.01, 0.0,  -0.05, 0.02, 0.01};
D1  = SE3PDdexp[X, U];
D1i = SE3PDdexpInv[X, U];
```

For $x\in\mathbb{R}^3$ with $n=\|x\|$, the rotation and translation blocks are computed via SO(3) sinc‑based formulas; for small $n$, the `eps` threshold avoids numerical issues.

## Notes

- Inputs are not shape‑validated; ensure twists are length‑6 and transforms are $4\times 4$ SE(3) elements.
- Works with numeric or symbolic data. Use `Chop[..., SE3EPS]` to suppress floating‑point noise.

## References

- A. Mueller, Review of the exponential and Cayley map on SE(3) as relevant for Lie group integration…, Proc. Royal Soc. A, 2021. doi: 10.1098/rspa.2021.0303
- Preprint: https://arxiv.org/abs/2303.07928
