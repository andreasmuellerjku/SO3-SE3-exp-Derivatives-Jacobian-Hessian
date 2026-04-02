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
- `SE3D2dexp[X, U, S, eps]`: Second Fréchet derivative $D^2_X(\mathrm{dexp})(U)(S)$ of $dexp$ at $X$.
- `SE3D2dexpInv[X, U, S, eps]`: Second Fréchet derivative $D^2_X(\mathrm{dexp}^{-1})(U)(S)$ of $dexp^{-1}$ at $X$.

Jacobians and Hessians:
- `SE3dexpJac[X, Z, eps]`: Jacobian of $dexp_X Z$ w.r.t. $X$.
- `SE3dexpInvJac[X, Z, eps]`: Jacobian of $dexp^{-1}_XZ$ w.r.t. $X$.
- `SE3dexpJacTranspose[X, Z, eps]`: Jacobian of $dexp^T_XZ$ w.r.t. $X$.
- `SE3dexpHessian[X, Q, Z, eps]`: Hessian of $Q^T dexp_X Z$ w.r.t. $X$.
- `SE3dexpInvHessian[X, Q, Z]`: Hessian of $Q^T dexp^{-1}_X Z$ w.r.t. $X$.

Series/truncated approximations (useful for small $\|x\|$):
- `SE3dexpOrder1[X]`,`SE3dexpOrder2[X]`,`SE3dexpOrder3[X]`,`SE3dexpOrder4[X]`: Approximation of $dexp$ of order $1,\ldots,4$.
- `SE3dexpInvOrder1[X]`,`SE3dexpInvOrder2[X]`,`SE3dexpInvOrder3[X]`,`SE3dexpInvOrder4[X]`: Approximation of $dexp^-1$ of order $1,\ldots,4$.
- `SE3DdexpOrder0[X, U]`,`SE3DdexpOrder1[X, U]`,`SE3DdexpOrder2[X, U]`,`SE3DdexpOrder3[X, U]`: Approximation of $D_X(dexp)(U)$ of order $0,\ldots,3$.
- `SE3DdexpInvOrder0[X, U]`,`SE3DdexpInvOrder1[X, U]`,`SE3DdexpInvOrder2[X, U]`,`SE3DdexpInvOrder3[X, U]`: Approximation of $D_X(dexp^{-1})(U)$ of order $0,\ldots,3$.
- `SE3D2dexpOrder0[X, U, S]`,`SE3D2dexpOrder1[X, U, S]`,`SE3D2dexpOrder2[X, U, S]`: Approximation of $D^2_X(dexp)(U)(S)$ of order $0,1,2$.
- `SE3D2dexpInvOrder0[X, U, S]`,`SE3D2dexpInvOrder2[X, U, S]`: Approximation of $D^2_X(dexp^{-1})(U)(S)$ of order $0,2$.
- `SE3dexpJacOrder0[X, Z]`,`SE3dexpJacOrder1[X, Z]`,`SE3dexpJacOrder2[X, Z]`: Approximation of the Jacobian of $dexp_XZ$ of order $0,1,2$.
- `SE3dexpHessianOrder0[X, Q, Z]`,`SE3dexpHessianOrder1[X, Q, Z]`,`SE3dexpHessianOrder2[X, Q, Z]`: Approximation of the Hessian of $Q^T dexp_XZ$ of order $0,1,2$.
- `SE3dexpHessianSwitchOrder2[X, Q, Z, eps]`: switch to exact Hessian when $\|x\|>\text{eps}$, else use order‑2 series.

Utilities:
- `SE3adBarMatrix[X]`: skew‑symmetric 6×6 matrix such that $\mathrm{adBar}_X\,U = \mathrm{ad}_U^\top X$.

Small‑angle handling uses `eps` (default `SE3EPS` from SE3Core`).

## Install

Place `SO3Core.m`, `SO3Exp.m`, `SE3Core.m`, and `SE3Exp.m` (this file) on your `$Path`, then:

```wl
<< SE3Exp`
```

## Usage

```wl
<< SE3Exp`

X = {0.1, -0.2, 0.05,  1.0, 0.3, -0.4};   (* twist {x,y} *)

(* Exponential on SE(3) *)
H = SE3Exp[X];

(* dexp and its inverse (6×6) *)
J   = SE3dexp[X];
Jinv = SE3dexpInv[X];
Chop[J . Jinv - IdentityMatrix[6], SE3EPS]

(* First derivative along U *)
U = {0.02, 0.01, 0.0, -0.05, 0.02, 0.01};
D1  = SE3Ddexp[X, U];
D1iinv = SE3DdexpInv[X, U];

(* Second derivative along U, S *)
U = {0.02, 0.01, 0.0, -0.05, 0.02, 0.01};
S = {1.7, -2.9, -9.2, 7.6, 6.7, 2.4};
D2  = SE3Ddexp[X, U, S];
D2iinv = SE3DdexpInv[X, U, S];

(* Jacobian and Hessian examples *)
Z   = {0.1, 0.0, -0.1, 0.2, -0.1, 0.0};
Q   = RandomReal[{-1,1}, 6];
Jxz = SE3dexpJac[X, Z];
Hzz = SE3dexpHessian[X, Q, Z];

(* Approximations *)
Js1 = SE3dexpOrder1[X];
Js2 = SE3dexpOrder2[X];
Js3 = SE3dexpOrder3[X];
Js4 = SE3dexpOrder4[X];

D10 = SE3DdexpOrder0[X, U, S];
D11 = SE3DdexpOrder1[X, U, S];
D12 = SE3DdexpOrder2[X, U, S];
D13 = SE3DdexpOrder3[X, U, S];

D20 = SE3D2dexpOrder0[X, U, S];
D21 = SE3D2dexpOrder1[X, U, S];
D22 = SE3D2dexpOrder2[X, U, S];

```

## Notes

- Inputs are not shape‑validated.
- Works with numeric or symbolic data. Use `Chop[..., SE3EPS]` to suppress floating‑point noise.
