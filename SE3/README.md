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

## Notes

- Inputs are not shape‑validated; ensure twists are length‑6, rotations are $3\times 3$, and transforms are $4\times 4$.
- Works with numeric or symbolic data. Use `Chop[..., SE3EPS]` to suppress floating‑point noise.
