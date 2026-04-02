# SE3Core

Minimal Wolfram Language utilities for working with the Lie algebra se(3) and the Lie group $\mathrm{SE}(3)$ (rigid‑body motions). Provides screw vector to $4\times 4$ matrix maps, the adjoint (Lie bracket) in 6‑vector form, homogeneous transform assembly, and group inverse.

- Author: Andreas Mueller
- Created: 2 Jan 2005
- Updated: 1 Apr 2026
- Depends on: SO3Core

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
<< SE3PExp`
```

## Usage

```wl
<< SE3PExp`

X = {0.1, -0.2, 0.05,  1.0, 0.3, -0.4};  (* twist {x, y} *)

(* Exponential and logarithm *)
H = SE3PExp[X];
Xrec = SE3PLog[H];

(* dexp and its inverse *)
J   = SE3Pdexp[X];
Jinv = SE3PdexpInv[X];
Chop[J . Jinv - IdentityMatrix[6], SE3EPS]

(* Derivatives *)
U = {0.02, 0.01, 0.0,  -0.05, 0.02, 0.01};
D1  = SE3PDdexp[X, U];
D1inv = SE3PDdexpInv[X, U];
```
## Notes

- Inputs are not shape‑validated; ensure twists are length‑6 and transforms are $4\times 4$ SE(3) elements.
- Works with numeric or symbolic data. Use `Chop[..., SE3EPS]` to suppress floating‑point noise.

## References

- A. Mueller: Review of the exponential and Cayley map on SE(3) as relevant for Lie group integration of the generalized Poisson equation and flexible multibody systems,  Proc. Royal Soc. A, September 2021; 477 (2253): 20210303 https://doi.org/10.1098/rspa.2021.0303
- Corrected preprint: https://arxiv.org/abs/2303.07928


# SE3Exp (Adjoint Representation)

Wolfram Language package for computing on SE(3) using canonical coordinates without block partitioning, via the 6×6 adjoint representation. Provides exponential/logarithm on SE(3), the differential of the exponential (dexp) and its inverse in closed form, first and second derivatives, Jacobians, Hessians, and higher-order approximations thereof.

- Author: Andreas Mueller
- Created: 9 Jun 2024
- Updated: 1 Apr 2026
- Depends on: SO3Core, SO3Exp, SE3Core

## Features

- `SE3Exp[X, eps]`: homogeneous transform $H=\exp(\widehat{X})$ using a stable polynomial in $\widehat{X}$.
- `SE3dexp[X, eps]`: right-trivialized differential of $dexp$ on $SE(3)$
- `SE3dexpInv[X, eps]`: inverse of $dexp$ on $SE(3)$
- `SE3Ddexp[X, U, eps]`: Fréchet derivative $D_X(\mathrm{dexp})(U)$ of $dexp$ at $X$ along $U$.
- `SE3DdexpInv[X, U, eps]`: Fréchet derivative $D_X(\mathrm{dexp}^{-1})(U)$ of $dexp^{-1}$ at $X$ along $U$.

Jacobians and Hessians:
- `SE3dexpJac[X, Z, eps]`, `SE3dexpInvJac[X, Z, eps]`, `SE3dexpJacTranspose[X, Z, eps]`
- `SE3dexpHessian[X, Q, Z, eps]`, `SE3dexpInvHessian[X, Q, Z]`
- `SE3dexpHessianSwitchOrder2[X, Q, Z, eps]`: switch to exact Hessian when $\|x\|>\text{eps}$, else use order‑2 series.

Series/truncated approximations (useful for very small $\|x\|$):
- `SE3dexpOrder{k}[X]`, `SE3dexpInvOrder{k}[X]` for $k=1..4$
- `SE3DdexpOrder{0|1|2|3}[X, U]`, `SE3DdexpInvOrder{0|1|2|3}[X, U]`
- `SE3D2dexpOrder{0|1|2}[X, U, S]`, `SE3D2dexpInvOrder{0|2}[X, U, S]`
- `SE3dexpJacOrder{0|1|2}[X, Z]`, `SE3dexpHessianOrder{0|1|2}[X, Q, Z]`

Utilities:
- `SE3adBarMatrix[X]`: skew‑symmetric 6×6 matrix such that $\mathrm{adBar}_X\,U = \mathrm{ad}_U^\top X$.

Small‑angle handling uses `eps` (default `SE3EPS` from SE3Core`).

## Install

Place `SO3Core.m`, `SO3Exp.m`, `SE3Core.m`, and `SE3Exp.m` (this file) on your `$Path`, then:

```wl
<< SO3Core`
<< SO3Exp`
<< SE3Core`
<< SE3Exp`
```

## Usage

```wl
<< SO3Core`
<< SO3Exp`
<< SE3Core`
<< SE3Exp`

X = {0.1, -0.2, 0.05,  1.0, 0.3, -0.4};   (* twist {x,y} *)

(* Exponential on SE(3) *)
H = SE3Exp[X];

(* dexp and its inverse (6×6) *)
J   = SE3dexp[X];
Jinv = SE3dexpInv[X];
Chop[J . Jinv - IdentityMatrix[6], SE3EPS]    (* 0 *)

(* First derivative along U *)
U = {0.02, 0.01, 0.0, -0.05, 0.02, 0.01};
D1  = SE3Ddexp[X, U];
D1i = SE3DdexpInv[X, U];

(* Jacobian and Hessian examples *)
Z   = {0.1, 0.0, -0.1, 0.2, -0.1, 0.0};
Q   = RandomReal[{-1,1}, 6];
Jxz = SE3dexpJac[X, Z];
Hzz = SE3dexpHessian[X, Q, Z];

(* Small-angle series as a fallback/check *)
Js2 = SE3dexpOrder2[X];
```

Notes:
- Inputs are not shape‑validated; ensure twists are length‑6 and transforms are $4\times 4$.
- Works with numeric or symbolic data. For very small $\|x\|$, prefer the provided series or use the `eps` threshold; use `Chop[..., SE3EPS]` to suppress numerical noise.
