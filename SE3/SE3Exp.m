(* ::Package:: *)

(* 
Author: Andreas Mueller
Created: 9 June 2024
Last updated: 1 April 2026
Purpose: Computing on SE(3) in terms of canonical coordinates without block partitioning using 6x6 adjoint representation
References: 
[1] A. Mueller: 
[2] A. Mueller: Review of the exponential and Cayley map on SE(3) as relevant for Lie group integration of the generalized Poisson equation and flexible multibody systems,  Proc. Royal Soc. A, September 2021; 477 (2253): 20210303. 
https://doi.org/10.1098/rspa.2021.0303
https://arxiv.org/abs/2303.07928
*)

BeginPackage["SE3Exp`","SE3Core`","SO3Core`","SO3Exp`"];

Unprotect[SE3adBarMatrix,SE3Exp,a1,a2,a3,a4,SE3dexp,b1,b2,b4,SE3dexpInv,a1bar,a2bar,a3bar,a4bar,Da1,Da2,Da3,Da4,SE3Ddexp,b2bar,b4bar,Db2,Db4,SE3DdexpInv,SE3dexpJac,SE3dexpInvJac,SE3dexpJacTranspose,DP1,DP2,DP3,DP4,D2P2,D2P3,D2P4,a1circ,a2circ,a3circ,a4circ,D2a1,D2a2,D2a3,D2a4,SE3D2dexp,b2circ,b4circ,D2b2,D2b4,SE3D2dexpInv,H11,H12,H13,H14,H22,H23,H24,H31,H32,H33,H34,SE3dexpHessian,H12Inv,H14Inv,H22Inv,H24Inv,H32Inv,H34Inv,SE3dexpInvHessian,SE3dexpOrder1,SE3dexpOrder2,SE3dexpOrder3,SE3dexpOrder4,SE3dexpInvOrder1,SE3dexpInvOrder2,SE3dexpInvOrder3,SE3dexpInvOrder4,SE3DdexpOrder0,SE3DdexpOrder1,SE3DdexpOrder2,SE3DdexpOrder3,SE3DdexpInvOrder0,SE3DdexpInvOrder1,SE3DdexpInvOrder2,SE3DdexpInvOrder3,SE3D2dexpOrder0,SE3D2dexpOrder1,SE3D2dexpOrder2,SE3D2dexpInvOrder0,SE3D2dexpInvOrder2,SE3dexpJacOrder0,SE3dexpJacOrder1,SE3dexpJacOrder2,SE3dexpHessianOrder0,SE3dexpHessianOrder1,SE3dexpHessianOrder2,SE3dexpHessianSwitchOrder2]; 

(* ==============================================
   SE(3) using the 6x6 ad matrix representation
   exp, log, dexp, its inverse, and derivatives
   ============================================== *)

SE3Exp::usage = "SE3Exp[X, eps] returns returns the SE(3) matrix exp(X) for twist X, with small-angle threshold eps.";

SE3dexp::usage = "SE3dexp[X, eps:-1] returns a polynomial approximation of dexp in the adjoint representation (6×6).";
SE3dexpInv::usage = "SE3dexpInv[X, eps:-1] returns a polynomial approximation of dexp^{-1} in the adjoint representation (6×6).";

SE3Ddexp::usage = "SE3Ddexp[X, U, eps:-1] returns the bilinear form for D(dexp) as a polynomial in ad_X.";
SE3DdexpInv::usage = "SE3DdexpInv[X, U, eps:-1] returns the bilinear form for D(dexp^{-1}) as a polynomial in ad_X.";

(* Jacobians with respect to X applied to Z *)
SE3dexpJac::usage = "SE3dexpJac[X, Z, eps:-1] returns the Jacobian (w.r.t X) of dexp[X]·Z.";
SE3dexpInvJac::usage = "SE3dexpInvJac[X, Z, eps:-1] returns the Jacobian (w.r.t X) of dexp[X]^{-1}·Z.";
SE3dexpJacTranspose::usage = "SE3dexpJacTranspose[X, Z, eps:-1] returns the Jacobian of dexp[X]^T·Z.";

(* Hessians (Q^T · dexp · Z and Q^T · dexp^{-1} · Z) *)
SE3dexpHessian::usage = "SE3dexpHessian[X, Q, Z, eps:-1] returns the full Hessian of Q^T·dexp[X]·Z (6×6).";
SE3dexpInvHessian::usage = "SE3dexpInvHessian[X, Q, Z] returns the full Hessian of Q^T·dexp[X]^{-1}·Z (6×6).";

(* Truncated series approximations *)
SE3dexpOrder1::usage = "SE3dexpOrder1[X] first-order series for dexp (6×6).";
SE3dexpOrder2::usage = "SE3dexpOrder2[X] second-order series for dexp (6×6).";
SE3dexpOrder3::usage = "SE3dexpOrder3[X] third-order series for dexp (6×6).";
SE3dexpOrder4::usage = "SE3dexpOrder4[X] fourth-order series for dexp (6×6).";

SE3dexpInvOrder1::usage = "SE3dexpInvOrder1[X] first-order series for dexp^{-1} (6×6).";
SE3dexpInvOrder2::usage = "SE3dexpInvOrder2[X] second-order series for dexp^{-1} (6×6).";
SE3dexpInvOrder3::usage = "SE3dexpInvOrder3[X] third-order series for dexp^{-1} (6×6).";
SE3dexpInvOrder4::usage = "SE3dexpInvOrder4[X] fourth-order series for dexp^{-1} (6×6).";

SE3DdexpOrder0::usage = "SE3DdexpOrder0[X, U] zeroth-order series for D(dexp).";
SE3DdexpOrder1::usage = "SE3DdexpOrder1[X, U] first-order series for D(dexp).";
SE3DdexpOrder2::usage = "SE3DdexpOrder2[X, U] second-order series for D(dexp).";
SE3DdexpOrder3::usage = "SE3DdexpOrder3[X, U] third-order series for D(dexp).";

SE3DdexpInvOrder0::usage = "SE3DdexpInvOrder0[X, U] zeroth-order series for D(dexp^{-1}).";
SE3DdexpInvOrder1::usage = "SE3DdexpInvOrder1[X, U] first-order series for D(dexp^{-1}).";
SE3DdexpInvOrder2::usage = "SE3DdexpInvOrder2[X, U] second-order series for D(dexp^{-1}).";
SE3DdexpInvOrder3::usage = "SE3DdexpInvOrder3[X, U] third-order series for D(dexp^{-1}).";

SE3D2dexpOrder0::usage = "SE3D2dexpOrder0[X, U, S] zeroth-order series for D^2(dexp).";
SE3D2dexpOrder1::usage = "SE3D2dexpOrder1[X, U, S] first-order series for D^2(dexp).";
SE3D2dexpOrder2::usage = "SE3D2dexpOrder2[X, U, S] second-order series for D^2(dexp).";

SE3D2dexpInvOrder0::usage = "SE3D2dexpInvOrder0[X, U, S] zeroth-order series for D^2(dexp^{-1}).";
SE3D2dexpInvOrder2::usage = "SE3D2dexpInvOrder2[X, U, S] second-order series for D^2(dexp^{-1}).";

SE3dexpJacOrder0::usage = "SE3dexpJacOrder0[X, Z] zeroth-order series for Jacobian of dexp·Z.";
SE3dexpJacOrder1::usage = "SE3dexpJacOrder1[X, Z] first-order series for Jacobian of dexp·Z.";
SE3dexpJacOrder2::usage = "SE3dexpJacOrder2[X, Z] second-order series for Jacobian of dexp·Z.";

SE3dexpHessianOrder0::usage = "SE3dexpHessianOrder0[X, Q, Z] zeroth-order series for Hessian of Q^T·dexp·Z.";
SE3dexpHessianOrder1::usage = "SE3dexpHessianOrder1[X, Q, Z] first-order series for Hessian of Q^T·dexp·Z.";
SE3dexpHessianOrder2::usage = "SE3dexpHessianOrder2[X, Q, Z] second-order series for Hessian of Q^T·dexp·Z.";
SE3dexpHessianSwitchOrder2::usage = "SE3dexpHessianSwitchOrder2[X, Q, Z, eps:-1] switches to exact Hessian when ||x|| > eps otherwise uses order-2 series.";



(* -------- SyntaxInformation for public API -------- *)

SyntaxInformation[SE3Exp] = {"ArgumentsPattern" -> {_}};

SyntaxInformation[SE3dexp] = {"ArgumentsPattern" -> {_, _.}};
SyntaxInformation[SE3dexpInv] = {"ArgumentsPattern" -> {_, _.}};

SyntaxInformation[SE3Ddexp] = {"ArgumentsPattern" -> {_, _, _.}};
SyntaxInformation[SE3DdexpInv] = {"ArgumentsPattern" -> {_, _, _.}};

SyntaxInformation[SE3dexpJac] = {"ArgumentsPattern" -> {_, _, _.}};
SyntaxInformation[SE3dexpInvJac] = {"ArgumentsPattern" -> {_, _, _.}};
SyntaxInformation[SE3dexpJacTranspose] = {"ArgumentsPattern" -> {_, _, _.}};

SyntaxInformation[SE3D2dexp] = {"ArgumentsPattern" -> {_, _, _, _.}};
SyntaxInformation[SE3D2dexpInv] = {"ArgumentsPattern" -> {_, _, _, _.}};

SyntaxInformation[SE3dexpHessian] = {"ArgumentsPattern" -> {_, _, _, _.}};

SyntaxInformation[SE3dexpInvHessian] = {"ArgumentsPattern" -> {_, _, _}};

SyntaxInformation[SE3dexpOrder1] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[SE3dexpOrder2] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[SE3dexpOrder3] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[SE3dexpOrder4] = {"ArgumentsPattern" -> {_}};

SyntaxInformation[SE3dexpInvOrder1] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[SE3dexpInvOrder2] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[SE3dexpInvOrder3] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[SE3dexpInvOrder4] = {"ArgumentsPattern" -> {_}};

SyntaxInformation[SE3DdexpOrder0] = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[SE3DdexpOrder1] = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[SE3DdexpOrder2] = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[SE3DdexpOrder3] = {"ArgumentsPattern" -> {_, _}};

SyntaxInformation[SE3DdexpInvOrder0] = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[SE3DdexpInvOrder1] = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[SE3DdexpInvOrder2] = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[SE3DdexpInvOrder3] = {"ArgumentsPattern" -> {_, _}};

SyntaxInformation[SE3D2dexpOrder0] = {"ArgumentsPattern" -> {_, _, _}};
SyntaxInformation[SE3D2dexpOrder1] = {"ArgumentsPattern" -> {_, _, _}};
SyntaxInformation[SE3D2dexpOrder2] = {"ArgumentsPattern" -> {_, _, _}};

SyntaxInformation[SE3D2dexpInvOrder0] = {"ArgumentsPattern" -> {_, _, _}};
SyntaxInformation[SE3D2dexpInvOrder2] = {"ArgumentsPattern" -> {_, _, _}};

SyntaxInformation[SE3dexpJacOrder0] = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[SE3dexpJacOrder1] = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[SE3dexpJacOrder2] = {"ArgumentsPattern" -> {_, _}};

SyntaxInformation[SE3dexpHessianOrder0] = {"ArgumentsPattern" -> {_, _, _}};
SyntaxInformation[SE3dexpHessianOrder1] = {"ArgumentsPattern" -> {_, _, _}};
SyntaxInformation[SE3dexpHessianOrder2] = {"ArgumentsPattern" -> {_, _, _}};
SyntaxInformation[SE3dexpHessianSwitchOrder2] = {"ArgumentsPattern" -> {_, _, _, _.}};


(* Begin private section of the package *)
Begin["`Private`"];

(* Input: X = {x, y} in R^6
Output: skew-symmetric 6×6 matrix adBar_X such that adBar_X.U = (ad_U)^T.X
*)
SE3adBarMatrix[X_] := Module[{x = Take[X, {1, 3}], y = Take[X, {4, 6}]}, 
       Return[ArrayFlatten[{{R3Toso3[x], R3Toso3[y]}, {R3Toso3[y], 0*IdentityMatrix[3]}}]]
];
Protect[SE3adBarMatrix]; 

(* SE3Exp - exponential map on SE(3)
Input: X = {x, y} in R^6, eps > 0
Output: H = exp(X^) in SE(3)
Reference: [1] eq. (17) *)
SE3Exp[X_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha, Delta, x, Xhat}, 
   x = X[[1 ;; 3]]; 
   n2 = x . x; 
   n = Sqrt[n2]; 
   If[n < eps, Return[RotPosToSE3[IdentityMatrix[3], X[[4 ;; 6]]]]];
   Alpha = Sinc[n]; 
   Beta = Sinc[n/2]^2; 
   Delta = (1 - Alpha)/n2;
   Xhat = TwistTose3[X]; 
   Return[IdentityMatrix[4] + Xhat + (Beta/2)*Xhat . Xhat + Delta*Xhat . Xhat . Xhat];
]; 
Protect[SE3Exp];

(* P - powers of ad_X
Input: X in R^6, power i
Output: P[X,i] = (ad_X)^i *)
Unprotect[P];
P[X_, i_] := MatrixPower[SE3adMatrix[X], i];
Protect[P]; 


(* SE3dexp - dexp on se(3)
Input: X = {x, y} in R^6, eps > 0
Output: 6×6 matrix dexp_X 
Reference: [1] eq. (21)  *)
SE3dexp[X_, eps_:SE3EPS] := Module[{n, n2, x},
   x = Take[X, {1, 3}];
   n2 = x . x; 
   n = Sqrt[n2];
   If[n < eps, Return[IdentityMatrix[6] + (1/2)*SE3adMatrix[X] + (1/6)*MatrixPower[SE3adMatrix[X], 2] + (1/24)*MatrixPower[SE3adMatrix[X], 3]]]; 
   Return[IdentityMatrix[6] + a1[x, eps]*P[X, 1] + a2[x, eps]*P[X, 2] + a3[x, eps]*P[X, 3] + a4[x, eps]*P[X, 4]];
];
Protect[SE3dexp];

(* Helper functions *)

a1[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha, Delta}, 
   n2 = x . x;
   n = Sqrt[n2]; 
   If[n < eps, Return[1/2]];
   Alpha = Sinc[n]; 
   Beta = Sinc[n/2]^2; 
   Delta = (1 - Alpha)/n2; 
   Return[Beta - Alpha/2];
]; 
Protect[a1];

a2[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha, Delta},
   n2 = x . x;
   n = Sqrt[n2]; 
   If[n < eps, Return[1/6]];
   Alpha = Sinc[n];
   Beta = Sinc[n/2]^2;
   Delta = (1 - Alpha)/n2;
   Return[(1/2)*(5*Delta - Beta/2)];
]; 
Protect[a2];

a3[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha, Delta},
   n2 = x . x; 
   n = Sqrt[n2]; 
   If[n < eps, Return[1/24]];
   Alpha = Sinc[n];
   Beta = Sinc[n/2]^2;
   Delta = (1 - Alpha)/n2;
   Return[(1/(2*n2))*(Beta - Alpha)];
]; 
Protect[a3];

a4[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha, Delta}, 
   n2 = x . x; 
   n = Sqrt[n2]; 
   If[n < eps, Return[1/120]]; 
   Alpha = Sinc[n]; 
   Beta = Sinc[n/2]^2; 
   Delta = (1 - Alpha)/n2; 
   Return[(1/(2*n2))*(3*Delta - Beta/2)];
]; 
Protect[a4];


(* SE3dexpInv - dexp on se(3)
Input: X = {x, y} in R^6, eps > 0
Output: 6×6 matrix dexp_X^-1 
Reference: [1] eq. (22)  *)
SE3dexpInv[X_, eps_:SE3EPS] := Module[{n, n2, Beta, x},
   x = Take[X, {1, 3}];
   n2 = x . x;
   n = Sqrt[n2]; 
   If[n < eps, Return[IdentityMatrix[6] - (1/2)*SE3adMatrix[X] + (1/6)*(1/2!)*MatrixPower[SE3adMatrix[X], 2]]]; 
   Return[IdentityMatrix[6] + b1[x]*P[X, 1] + b2[x, eps]*P[X, 2] + b4[x, eps]*P[X, 4]];
];
Protect[SE3dexpInv];

(* Helper functions *)

b1[x_] := -1/2
Protect[b1];

b2[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha},
   n2 = x . x;
   n = Sqrt[n2]; 
   If[n < eps, Return[1/12]]; 
   Alpha = Sinc[n]; 
   Beta = Sinc[n/2]^2;
   Return[(1/n2)*(2 - (1 + 3*Alpha)/(2*Beta))];
]; 
Protect[b2];

b4[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha},
   n2 = x . x; 
   n = Sqrt[n2]; 
   If[n < eps, Return[-1/720]];
   Alpha = Sinc[n]; 
   Beta = Sinc[n/2]^2;
   Return[(1/n2^2)*(1 - (1 + Alpha)/(2*Beta))];
];
Protect[b4];


(* SE3Ddexp - derivative of dexp on se(3)
Input: X, U in R^6, eps > 0
Output: 6×6 matrix D(dexp)_X(U)
Reference: [1] eq. (28) *)
SE3Ddexp[X_, U_, eps_:SE3EPS] := Module[{x, u, adU = SE3adMatrix[U], adX = SE3adMatrix[X], P1, P2}, 
   x = Take[X, {1, 3}]; 
   u = Take[U, {1, 3}]; 
   P1 = SE3adMatrix[X]; 
   P2 = SE3adMatrix[X] . P1; 
   If[x . x < eps, Return[(1/2)*adU + (1/6)*(adU . P1 + P1 . adU) + (1/24)*(P2 . adU + adU . P2 + P1 . adU . P1)]]; 
   Return[a1[x, eps]*P[U, 1] + a2[x, eps]*(P[X, 1] . P[U, 1] + P[U, 1] . P[X, 1]) + 
   a3[x, eps]*(P[X, 2] . P[U, 1] + P[U, 1] . P[X, 2] + P[X, 1] . P[U, 1] . P[X, 1]) + 
   a4[x, eps]*(P[X, 3] . P[U, 1] + P[X, 1] . P[U, 1] . P[X, 2] + P[X, 2] . P[U, 1] . P[X, 1] + P[U, 1] . P[X, 3]) + 
   Da1[X, U, eps]*P[X, 1] + Da2[X, U, eps]*P[X, 2] + Da3[X, U, eps]*P[X, 3] + Da4[X, U, eps]*P[X, 4]];
];
Protect[SE3Ddexp];

(* Helper functions *)

a1bar[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha, Delta}, 
   n2 = x . x; n = Sqrt[n2]; 
   If[n < eps, Return[0]]; 
   Alpha = Sinc[n]; 
   Beta = Sinc[n/2]^2; 
   Delta = (1 - Alpha)/n2; 
   Return[(1/4)*Beta + (1/n2)*(2*(Alpha - Beta) + Alpha/2 - 1/2)];
]; 
Protect[a1bar];

a2bar[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha, Delta}, 
   n2 = x . x; 
   n = Sqrt[n2];
   If[n < eps, Return[0]]; 
   Alpha = Sinc[n]; 
   Beta = Sinc[n/2]^2; 
   Delta = (1 - Alpha)/n2; 
   Return[(1/n2)*((7/4)*Beta - (1/2)*Alpha - (15/2)*Delta)];
]; 
Protect[a2bar];

a3bar[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha, Delta}, 
   n2 = x . x; 
   n = Sqrt[n2]; 
   If[n < eps, Return[-180^(-1)]]; 
   Alpha = Sinc[n]; 
   Beta = Sinc[n/2]^2; 
   Delta = (1 - Alpha)/n2; 
   Return[(1/n2^2)*((5/2)*Alpha - 1/2 + (n2/4 - 2)*Beta)];
]; 
Protect[a3bar];

a4bar[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha, Delta}, 
   n2 = x . x; 
   n = Sqrt[n2]; 
   If[n < eps, Return[-1260^(-1)]];
   Alpha = Sinc[n]; 
   Beta = Sinc[n/2]^2; 
   Delta = (1 - Alpha)/n2; 
   Return[(1/n2^3)*((15/2 - n2/2)*Alpha - 15/2 + (7/4)*n2*Beta)];
]; 
Protect[a4bar];

Da1[X_, U_, eps_:SE3EPS] := Module[{x, u}, 
   x = Take[X, {1, 3}]; 
   u = Take[U, {1, 3}]; 
   Return[a1bar[x, eps]*x . u];
]; 
Protect[Da1];

Da2[X_, U_, eps_:SE3EPS] := Module[{x, u}, 
   x = Take[X, {1, 3}]; 
   u = Take[U, {1, 3}]; 
   Return[a2bar[x, eps]*x . u];
]; 
Protect[Da2];

Da3[X_, U_, eps_:SE3EPS] := Module[{x, u}, 
   x = Take[X, {1, 3}]; 
   u = Take[U, {1, 3}]; 
   Return[a3bar[x, eps]*x . u];
]; 
Protect[Da3];

Da4[X_, U_, eps_:SE3EPS] := Module[{x, u}, 
   x = Take[X, {1, 3}]; 
   u = Take[U, {1, 3}]; 
   Return[a4bar[x, eps]*x . u];
]; 
Protect[Da4];



(* SE3DdexpInv - derivative of dexp^-1 on se(3)
Input: X, U in R^6, eps > 0
Output: 6×6 matrix D(dexp)^-1_X(U)
Reference: [1] eq. (29) *)
SE3DdexpInv[X_, U_, eps_:SE3EPS] := Module[{x, u},
   x = Take[X, {1, 3}];
   u = Take[U, {1, 3}];
   If[x . x < eps, Return[(-1/2)*SE3adMatrix[U]]];
   Return[b1[x]*P[U, 1] + b2[x, eps]*(P[X, 1] . P[U, 1] + P[U, 1] . P[X, 1])
   + b4[x, eps]*(P[X, 3] . P[U, 1] + P[X, 1] . P[U, 1] . P[X, 2] + P[X, 2] . P[U, 1] . P[X, 1] + P[U, 1] . P[X, 3]) 
   + Db2[X, U, eps]*P[X, 2] + Db4[X, U, eps]*P[X, 4]];
];
Protect[SE3DdexpInv];

(* Helper functions *)

b2bar[X_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha, Delta, x, Gamma}, 
   x = Take[X, {1, 3}]; 
   n2 = x . x; 
   n = Sqrt[n2]; 
   If[n < eps, Return[0]];
   Alpha = Sinc[n];
   Beta = Sinc[n/2]^2;
   Delta = (1 - Alpha)/n2;
   Gamma = Alpha/Beta; 
   Return[(1/n2)*(3/4) + (1/n2^2)*(Gamma*(3/2 + 1/Beta + 3*Gamma) - 3/(2*Beta) - 4)];
]; 
Protect[b2bar];

b4bar[X_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha, Delta, x, Gamma}, 
   x = Take[X, {1, 3}]; 
   n2 = x . x; 
   n = Sqrt[n2]; 
   If[n < eps, Return[-1/7560]]; 
   Alpha = Sinc[n]; 
   Beta = Sinc[n/2]^2; 
   Delta = (1 - Alpha)/n2; 
   Gamma = Alpha/Beta; 
   Return[1/(4*n2^2) + (1/n2^3)*(1/(2*Beta) + (3/2)*Gamma + Gamma/Beta + Gamma^2 - 4)];
]; 
Protect[b4bar];

Db2[X_, U_, eps_:SE3EPS] := Module[{x, u}, 
   x = Take[X, {1, 3}]; 
   u = Take[U, {1, 3}]; 
   Return[b2bar[x, eps]*x . u];
]; 
Protect[Db2];

Db4[X_, U_, eps_:SE3EPS] := Module[{x, u}, 
   x = Take[X, {1, 3}]; 
   u = Take[U, {1, 3}]; 
   Return[b4bar[x, eps]*x . u]]; 
Protect[Db4];



(* SE3dexpJac - Jacobian of dexp[X]·Z w.r.t X
Input: X, Z in R^6, eps > 0
Output: 6×6 matrix J = d(dexp[X] Z)/dX 
Reference: [1]  eq. ??? *)
SE3dexpJac[X_, Z_, eps_:SE3EPS] := Module[{x, Z1, Z2, Z3, Z4}, 
   x = Take[X, {1, 3}]; 
   Z1 = P[X, 1] . Z; 
   Z2 = P[X, 2] . Z; 
   Z3 = P[X, 3] . Z; Z4 = P[X, 4] . Z; If[x . x < eps, Return[(-2^(-1))*SE3adMatrix[Z]]]; 
   Return[(-a1[x, eps])*P[Z, 1] + a2[x, eps]*(P[Z, 1] . P[X, 1] - 2*P[X, 1] . P[Z, 1]) 
   - a3[x, eps]*(P[X, 2] . P[Z, 1] + P[X, 1] . P[Z1, 1] + P[Z2, 1]) 
   - a4[x, eps]*(P[X, 3] . P[Z, 1] + P[X, 1] . P[Z2, 1] + P[X, 2] . P[Z1, 1] + P[Z3, 1]) 
   + a1bar[x, eps]*Transpose[{Z1}] . {Join[x, {0, 0, 0}]} 
   + a2bar[x, eps]*Transpose[{Z2}] . {Join[x, {0, 0, 0}]} 
   + a3bar[x, eps]*Transpose[{Z3}] . {Join[x, {0, 0, 0}]} 
   + a4bar[x, eps]*Transpose[{Z4}] . {Join[x, {0, 0, 0}]}];
];
Protect[SE3dexpJac];


(* SE3dexpInvJac - Jacobian w.r.t X of dexp[X]^{-1}·Z
Input: X, Z in R^6, eps > 0
Output: 6×6 matrix J = d(dexp[X]^-1 Z)/dX
Reference: [1]  eq. ??? *)
SE3dexpInvJac[X_, Z_, eps_:SE3EPS] := Module[{x, u, Z1, Z2, Z3, Z4}, 
   x = Take[X, {1, 3}]; 
   Z1 = P[X, 1] . Z; 
   Z2 = P[X, 2] . Z; 
   Z3 = P[X, 3] . Z; 
   Z4 = P[X, 4] . Z; 
   If[x . x < eps, Return[(1/2)*SE3adMatrix[Z] - (1/12)*(SE3adMatrix[Z1] + P[X, 1] . Z)]]; 
   Return[(-b1[x])*P[Z, 1] + b2[x, eps]*(P[Z, 1] . P[X, 1] - 2*P[X, 1] . P[Z, 1]) 
   - b4[x, eps]*(P[X, 3] . P[Z, 1] + P[X, 1] . P[Z2, 1] + P[X, 2] . P[Z1, 1] + P[Z3, 1]) 
   + b2bar[x, eps]*Transpose[{Z2}] . {Join[x, {0, 0, 0}]} 
   + b4bar[x, eps]*Transpose[{Z4}] . {Join[x, {0, 0, 0}]}];
];
Protect[SE3dexpInvJac];


(* SE3dexpJacTranspose - Jacobian of dexp[X]^T·Z
Input: X, Z in R^6, eps > 0
Output: 6×6 matrix J = d(dexp[X]^-T Z)/dX
Reference: [1]  eq. ??? *)
SE3dexpJacTranspose[X_, Z_, eps_:SE3EPS] := Module[{x, u, Z1, Z2, Z3, Z4, P1, P2, P3}, 
   x = Take[X, {1, 3}];
   Z1 = Transpose[P[X, 1]] . Z;
   Z2 = Transpose[P[X, 2]] . Z;
   Z3 = Transpose[P[X, 3]] . Z;
   Z4 = Transpose[P[X, 4]] . Z; 
   P1 = P[X, 1];
   P2 = P[X, 2];
   P3 = P[X, 3]; 
   If[x . x < eps, Return[(1/2)*SE3adBarMatrix[Z]]]; 
   Return[a1[x, eps]*SE3adBarMatrix[Z] + a2[x, eps]*(SE3adBarMatrix[Z1] + Transpose[P1] . SE3adBarMatrix[Z]) 
   + a3[x, eps]*(SE3adBarMatrix[Z2] + Transpose[P2] . SE3adBarMatrix[Z] + Transpose[P1] . SE3adBarMatrix[Z1]) 
   + a4[x, eps]*(SE3adBarMatrix[Z3] + Transpose[P2] . SE3adBarMatrix[Z1] + Transpose[P1] . SE3adBarMatrix[Z2] + Transpose[P3] . SE3adBarMatrix[Z]) 
   + a1bar[x, eps]*Transpose[{Z1}] . {Join[x, {0, 0, 0}]}
   + a2bar[x, eps]*Transpose[{Z2}] . {Join[x, {0, 0, 0}]} 
   + a3bar[x, eps]*Transpose[{Z3}] . {Join[x, {0, 0, 0}]} 
   + a4bar[x, eps]*Transpose[{Z4}] . {Join[x, {0, 0, 0}]}];
];
Protect[SE3dexpJacTranspose];



(* SE3D2dexp - second directional derivative of dexp on se(3)
Input: X, U, S in R^6, eps > 0
Output: 6×6 matrix (D^2dexp_X)(U,S)
Reference: [1]  eq. ??? *)
SE3D2dexp[X_, U_, S_, eps_:SE3EPS] := Module[{x, u, s, adU = SE3adMatrix[U], adS = SE3adMatrix[S], adX = SE3adMatrix[X], P1, P2},
   x = Take[X, {1, 3}];
   u = Take[U, {1, 3}];
   s = Take[S, {1, 3}];
   P1 = SE3adMatrix[X];
   P2 = SE3adMatrix[X] . P1;
   If[x . x < eps, Return[(1/6)*(adU . adS + adS . adU) 
   + (1/24)*(P1 . adS + adS . P1) . adU 
   + (1/24)*(adU . adS + adS . adU) . P1 
   + (1/24)*(P1 . adU + adU . P1) . adS 
   + (1/120)*adS . (P2 . adU + adU . P2 + P1 . adU . P1) 
   + (1/120)*adU . (P2 . adS + adS . P2 + P1 . adS . P1) 
   + (1/120)*P2 . (adS . adU + adU . adS) 
   + (1/120)*P1 . (adS . (P1 . adU + adU . P1) 
   +  adU . (P1 . adS + adS . P1))]];
   Return[D2a1[X, U, S, eps]*P[X, 1] + D2a2[X, U, S, eps]*P[X, 2] + D2a3[X, U, S, eps]*P[X, 3] + D2a4[X, U, S]*P[X, 4] 
   + a2[x, eps]*D2P2[X, U, S] + a3[x, eps]*D2P3[X, U, S] + a4[x, eps]*D2P4[X, U, S] 
   + Da1[X, U, eps]*DP1[X, S] + Da2[X, U, eps]*DP2[X, S] + Da3[X, U, eps]*DP3[X, S]
   + Da4[X, U, eps]*DP4[X, S] + Da1[X, S, eps]*DP1[X, U] + Da2[X, S, eps]*DP2[X, U] 
   + Da3[X, S, eps]*DP3[X, U] + Da4[X, S, eps]*DP4[X, U]];
];
Protect[SE3D2dexp];

(* Helper functions *)

DP1[X_, U_] := SE3adMatrix[U]
Protect[DP1];

DP2[X_, U_] := SE3adMatrix[X] . SE3adMatrix[U] + SE3adMatrix[U] . SE3adMatrix[X]
Protect[DP2];

DP3[X_, U_] := SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[U] + 
   SE3adMatrix[U] . SE3adMatrix[X] . SE3adMatrix[X] + 
   SE3adMatrix[X] . SE3adMatrix[U] . SE3adMatrix[X]
Protect[DP3];

DP4[X_, U_] := SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[U] + 
   SE3adMatrix[U] . SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[X] + 
   SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[U] . SE3adMatrix[X] + 
   SE3adMatrix[X] . SE3adMatrix[U] . SE3adMatrix[X] . SE3adMatrix[X]
Protect[DP4];

D2P2[X_, U_, S_] := SE3adMatrix[S] . SE3adMatrix[U] + SE3adMatrix[U] . SE3adMatrix[S]
Protect[D2P2];

D2P3[X_, U_, S_] := (SE3adMatrix[S] . SE3adMatrix[X] + SE3adMatrix[X] . SE3adMatrix[S]) . SE3adMatrix[U] +
   (SE3adMatrix[S] . SE3adMatrix[U] + SE3adMatrix[U] . SE3adMatrix[S]) . SE3adMatrix[X] +
   (SE3adMatrix[U] . SE3adMatrix[X] + SE3adMatrix[X] . SE3adMatrix[U]) . SE3adMatrix[S]
Protect[D2P3];

D2P4[X_, U_, S_] := SE3adMatrix[S] . (SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[U] +
   SE3adMatrix[U] . SE3adMatrix[X] . SE3adMatrix[X] +
   SE3adMatrix[X] . SE3adMatrix[U] . SE3adMatrix[X]) +
   SE3adMatrix[U] . (SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[S] +
   SE3adMatrix[S] . SE3adMatrix[X] . SE3adMatrix[X] +
   SE3adMatrix[X] . SE3adMatrix[S] . SE3adMatrix[X]) +
   SE3adMatrix[X] . SE3adMatrix[X] . (SE3adMatrix[S] . SE3adMatrix[U] + SE3adMatrix[U] . SE3adMatrix[S]) +
   SE3adMatrix[X] . (SE3adMatrix[S] . (SE3adMatrix[U] . SE3adMatrix[X] + SE3adMatrix[X] . SE3adMatrix[U]) +
   SE3adMatrix[U] . (SE3adMatrix[S] . SE3adMatrix[X] + SE3adMatrix[X] . SE3adMatrix[S]))
Protect[D2P4];

a1circ[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha}, 
    n2 = x . x;
    n = Sqrt[n2];
    If[n < eps, Return[-1/90]];
    Alpha = Sinc[n];
    Beta = Sinc[n/2]^2;
    Return[(1/n2)*(Alpha/2 - (7/4)*Beta) + (1/n2^2)*(7/2 - (23/2)*Alpha + 8*Beta)];
]; 
Protect[a1circ];

a2circ[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha}, 
   n2 = x . x;
   n = Sqrt[n2];
   If[n < eps, Return[-1/630]]; 
   Alpha = Sinc[n];
   Beta = Sinc[n/2]^2;
   Return[(1/n2)*(1/4)*Beta + (1/n2^2)*(5*Alpha - (43/4)*Beta - 1/2) + (1/n2^3)*(75/2 - (75/2)*Alpha)];
]; 
Protect[a2circ];

a3circ[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha}, 
   n2 = x . x;
   n = Sqrt[n2];
   If[n < eps, Return[1/1680]]; 
   Alpha = Sinc[n];
   Beta = Sinc[n/2]^2;
   Return[(1/n2^2)*((1/2)*Alpha - (9/4)*Beta) + (1/n2^3)*(9/2 - (33/2)*Alpha + 12*Beta)];
]; 
Protect[a3circ];

a4circ[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha, Delta}, 
   n2 = x . x;
   n = Sqrt[n2];
   If[n < eps, Return[1/15120]]; 
   Alpha = Sinc[n];
   Beta = Sinc[n/2]^2;
   Delta = (1 - Alpha)/n2; 
   Return[(1/n2^2)*(1/4)*Beta + (1/n2^3)*(6*Alpha - 1/2 - (57/4)*Beta) + (1/n2^4)*(105/2 - (105/2)*Alpha)];
]; 
Protect[a4circ];

D2a1[X_, U_, S_, eps_:SE3EPS] := Module[{x, u, s},
   x = Take[X, {1, 3}];
   u = Take[U, {1, 3}];
   s = Take[S, {1, 3}]; 
   Return[a1bar[x, eps]*s . u + a1circ[x, eps]*x . u*x . s];
]; 
Protect[D2a1];

D2a2[X_, U_, S_, eps_:SE3EPS] := Module[{x, u, s}, 
   x = Take[X, {1, 3}];
   u = Take[U, {1, 3}]; 
   s = Take[S, {1, 3}]; 
   Return[a2bar[x, eps]*s . u + a2circ[x, eps]*x . u*x . s];
]; 
Protect[D2a2];

D2a3[X_, U_, S_, eps_:SE3EPS] := Module[{x, u, s}, 
   x = Take[X, {1, 3}]; 
   u = Take[U, {1, 3}]; 
   s = Take[S, {1, 3}]; 
   Return[a3bar[x, eps]*s . u + a3circ[x, eps]*x . u*x . s];
]; 
Protect[D2a3];

D2a4[X_, U_, S_, eps_:SE3EPS] := Module[{x, u, s},
   x = Take[X, {1, 3}];
   u = Take[U, {1, 3}];
   s = Take[S, {1, 3}]; 
   Return[a4bar[x, eps]*s . u + a4circ[x, eps]*x . u*x . s];
]; 
Protect[D2a4];



(* SE3D2dexpInv - second directional derivative of dexp on se(3)
Input: X, U, S in R^6, eps > 0
Output: 6×6 matrix (D^2dexp^-1_X)(U,S)
Reference: [1]  eq. ??? *)
SE3D2dexpInv[X_, U_, S_, eps_:SE3EPS] := Module[{x, u, s, adU = SE3adMatrix[U], adS = SE3adMatrix[S], 
    adX = SE3adMatrix[X]},
    x = Take[X, {1, 3}];
    u = Take[U, {1, 3}];
    s = Take[S, {1, 3}]; 
    If[x . x < eps, Return[(1/12)*(adU . adS + adS . adU)]];
    Return[D2b2[X, U, S, eps]*P[X, 2] + D2b4[X, U, S, eps]*P[X, 4] 
    + b2[x, eps]*D2P2[X, U, S] + b4[x, eps]*D2P4[X, U, S] 
    + Db2[X, U, eps]*DP2[X, S] + Db4[X, U, eps]*DP4[X, S] 
    + Db2[X, S, eps]*DP2[X, U] + Db4[X, S, eps]*DP4[X, U]];
];
Protect[SE3D2dexpInv];

(* Helper functions *)

b2circ[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha, Gamma}, 
   n2 = x . x; 
   n = Sqrt[n2];
   If[n < eps, Return[-1/3780]];
   Alpha = Sinc[n];
   Beta = Sinc[n/2]^2;
   Gamma = Alpha/Beta; 
   Return[(1/n2^3)*(16 + 1/Beta^2 + (8/Beta - 9/2)*Gamma - (9 + 4/Beta)*Gamma^2 - 12*Gamma^3 + 9/(2*Beta)) - (1/n2^2)*(9/4 + 3*Gamma + 1/(2*Beta))];
]; 
Protect[b2circ];

b4circ[x_, eps_:SE3EPS] := Module[{n, n2, Beta, Alpha, Gamma},
   n2 = x . x;
   n = Sqrt[n2];
   If[n < eps, Return[-1/50400]];
   Alpha = Sinc[n];
   Beta = Sinc[n/2]^2;
   Gamma = Alpha/Beta; 
   Return[(1/n2^4)*(24 - 2/Beta - (15/2 + 4/Beta)*Gamma - (11/2 + 3/Beta)*Gamma^2 - 2*Gamma^3) 
   - (1/n2^3)*(11/8 + 1/(4*Beta) + (1/2)*Gamma)];
]; 
Protect[b4circ];

D2b2[X_, U_, S_, eps_:SE3EPS] := Module[{x, u, s}, 
   x = Take[X, {1, 3}];
   u = Take[U, {1, 3}];
   s = Take[S, {1, 3}];
   Return[b2bar[x, eps]*s . u + b2circ[x, eps]*x . u*x . s];
]; 
Protect[D2b2];

D2b4[X_, U_, S_, eps_:SE3EPS] := Module[{x, u, s},
   x = Take[X, {1, 3}];
   u = Take[U, {1, 3}];
   s = Take[S, {1, 3}]; 
   Return[b4bar[x, eps]*s . u + b4circ[x, eps]*x . u*x . s];
]; 
Protect[D2b4];



(* SE3dexpHessian - Hessian of Q^T·dexp[X]·Z
Input: X, Q, Z in R^6
Output: 6×6 Hessian
Reference: [1]  eq. ??? *)
SE3dexpHessian[X_, Q_, Z_, eps_:SE3EPS] := H11[X, Q, Z, eps] + H12[X, Q, Z, eps] + H13[X, Q, Z, eps] + H14[X, Q, Z, eps] + 
   H22[X, Q, Z, eps] + Transpose[H22[X, Q, Z, eps]] + H23[X, Q, Z, eps] + Transpose[H23[X, Q, Z, eps]] + 
   H24[X, Q, Z, eps] + Transpose[H24[X, Q, Z, eps]] + H31[X, Q, Z, eps] + Transpose[H31[X, Q, Z, eps]] + 
   H32[X, Q, Z, eps] + Transpose[H32[X, Q, Z, eps]] + H33[X, Q, Z, eps] + Transpose[H33[X, Q, Z, eps]] + 
   H34[X, Q, Z, eps] + Transpose[H34[X, Q, Z, eps]]
Protect[SE3dexpHessian];

(* Helper functions *)

H11[X_, Q_, Z_, eps_:SE3EPS] := Module[{x, zero3}, 
   x = Take[X, {1, 3}];
   zero3 = 0*IdentityMatrix[3]; 
   Return[Q . P[X, 1] . Z*ArrayFlatten[{{a1bar[x, eps]*IdentityMatrix[3] + a1circ[x, eps]*Transpose[{x}] . {x}, zero3}, {zero3, zero3}}]];
];
Protect[H11];

H12[X_, Q_, Z_, eps_:SE3EPS] := Module[{x, zero3}, 
   x = Take[X, {1, 3}];
   zero3 = 0*IdentityMatrix[3]; 
   Return[Q . P[X, 2] . Z*ArrayFlatten[{{a2bar[x, eps]*IdentityMatrix[3] + a2circ[x, eps]*Transpose[{x}] . {x}, zero3}, {zero3, zero3}}]];
];
Protect[H12];

H13[X_, Q_, Z_, eps_:SE3EPS] := Module[{x, zero3},
   x = Take[X, {1, 3}];
   zero3 = 0*IdentityMatrix[3]; 
   Return[Q . P[X, 3] . Z*ArrayFlatten[{{a3bar[x, eps]*IdentityMatrix[3] + a3circ[x, eps]*Transpose[{x}] . {x}, zero3}, {zero3, zero3}}]];
];
Protect[H13];

H14[X_, Q_, Z_, eps_:SE3EPS] := Module[{x, zero3},
   x = Take[X, {1, 3}];
   zero3 = 0*IdentityMatrix[3];
   Return[Q . P[X, 4] . Z*ArrayFlatten[{{a4bar[x, eps]*IdentityMatrix[3] + a4circ[x, eps]*Transpose[{x}] . {x}, zero3}, {zero3, zero3}}]];
];
Protect[H14];

H22[X_, Q_, Z_, eps_:SE3EPS] := a2[Take[X, {1, 3}], eps]*SE3adBarMatrix[Q] . SE3adMatrix[Z]; 
Protect[H22];

H23[X_, Q_, Z_, eps_:SE3EPS] := Module[{x, Z1, Z2, Q1}, 
   x = Take[X, {1, 3}]; Z1 = P[X, 1] . Z; Z2 = P[X, 2] . Z; 
   Q1 = Transpose[P[X, 1]] . Q; 
   Return[a3[x, eps]*(SE3adBarMatrix[Q1] . SE3adMatrix[Z] + SE3adBarMatrix[Q] . SE3adMatrix[X] . SE3adMatrix[Z] + SE3adBarMatrix[Q] . SE3adMatrix[Z1])];
];
Protect[H23];

H24[X_, Q_, Z_, eps_:SE3EPS] := Module[{x, u, Z1, Z2, Q1, Q2},
   x = Take[X, {1, 3}];
   Z1 = P[X, 1] . Z;
   Z2 = P[X, 2] . Z; 
   Q1 = Transpose[P[X, 1]] . Q;
   Q2 = Transpose[P[X, 2]] . Q; 
   Return[a4[x, eps]*(SE3adBarMatrix[Q2] . SE3adMatrix[Z] + SE3adBarMatrix[Q1] . (SE3adMatrix[Z1] + SE3adMatrix[X] . SE3adMatrix[Z]) 
   + SE3adBarMatrix[Q] . (SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[Z] + SE3adMatrix[X] . SE3adMatrix[Z1] + SE3adMatrix[Z2]))];
];
Protect[H24];

H31[X_, Q_, Z_, eps_:SE3EPS] := Module[{x}, 
   x = Take[X, {1, 3}];
   Return[(-a1bar[x, eps])*(Transpose[{Join[x, {0, 0, 0}]}] . {Q}) . SE3adMatrix[Z]];
];
Protect[H31];

H32[X_, Q_, Z_, eps_:SE3EPS] := Module[{x}, 
   x = Take[X, {1, 3}]; 
   Return[a2bar[x, eps]*(Transpose[{Join[x, {0, 0, 0}]}] . {Q}) . (SE3adMatrix[Z] . SE3adMatrix[X] - 2*SE3adMatrix[X] . SE3adMatrix[Z])];
];
Protect[H32];

H33[X_, Q_, Z_, eps_:SE3EPS] := Module[{x, Z1, Z2}, 
   x = Take[X, {1, 3}];
   Z1 = P[X, 1] . Z;
   Z2 = P[X, 2] . Z;
   Return[(-a3bar[x, eps])*(Transpose[{Join[x, {0, 0, 0}]}] . {Q}) . (SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[Z] + 
   SE3adMatrix[X] . SE3adMatrix[Z1] + SE3adMatrix[Z2])];
];
Protect[H33];

H34[X_, Q_, Z_, eps_:SE3EPS] := Module[{x, Z1, Z2, Z3},
   x = Take[X, {1, 3}]; 
   Z1 = P[X, 1] . Z;
   Z2 = P[X, 2] . Z; 
   Z3 = P[X, 3] . Z;
   Return[(-a4bar[x, eps])*(Transpose[{Join[x, {0, 0, 0}]}] . {Q}) . (SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[Z] 
   + SE3adMatrix[X] . SE3adMatrix[Z2] + SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[Z1] + SE3adMatrix[Z3])];
];
Protect[H34];



(* SE3dexpInvHessian - Hessian of Q^T·dexp[X]^{-1}·Z
Input: X, Q, Z in R^6
Output: 6×6 Hessian
Reference: [1]  eq. ??? *)
SE3dexpInvHessian[X_, Q_, Z_] := H12Inv[X, Q, Z] + H14Inv[X, Q, Z] + 
   H22Inv[X, Q, Z] + Transpose[H22Inv[X, Q, Z]] +
   H24Inv[X, Q, Z] + Transpose[H24Inv[X, Q, Z]] +
   H32Inv[X, Q, Z] + Transpose[H32Inv[X, Q, Z]] +
   H34Inv[X, Q, Z] + Transpose[H34Inv[X, Q, Z]]
Protect[SE3dexpInvHessian];

(* Helper functions *)

H12Inv[X_, Q_, Z_, eps_:SE3EPS] := Module[{x, u, zero3}, 
   x = Take[X, {1, 3}];
   zero3 = 0*IdentityMatrix[3]; 
   Return[Q . P[X, 2] . Z*ArrayFlatten[{{b2bar[x, eps]*IdentityMatrix[3] + b2circ[x, eps]*Transpose[{x}] . {x}, zero3}, {zero3, zero3}}]];
];
Protect[H12Inv];

H14Inv[X_, Q_, Z_, eps_:SE3EPS] := Module[{x, u, zero3},
   x = Take[X, {1, 3}];
   zero3 = 0*IdentityMatrix[3]; 
   Return[Q . P[X, 4] . Z*ArrayFlatten[{{b4bar[x, eps]*IdentityMatrix[3] + b4circ[x, eps]*Transpose[{x}] . {x}, zero3}, {zero3, zero3}}]];
];
Protect[H14Inv];

H22Inv[X_, Q_, Z_, eps_:SE3EPS] := b2[Take[X, {1, 3}], eps]*SE3adBarMatrix[Q] . SE3adMatrix[Z]; 
Protect[H22Inv];

H24Inv[X_, Q_, Z_, eps_:SE3EPS] := Module[{x, u, Z1, Z2, Q1, Q2}, 
   x = Take[X, {1, 3}];
   Z1 = P[X, 1] . Z;
   Z2 = P[X, 2] . Z; 
   Q1 = Transpose[P[X, 1]] . Q;
   Q2 = Transpose[P[X, 2]] . Q; 
   Return[b4[x, eps]*(SE3adBarMatrix[Q2] . SE3adMatrix[Z] + SE3adBarMatrix[Q1] . (SE3adMatrix[Z1] + SE3adMatrix[X] . SE3adMatrix[Z]) 
   + SE3adBarMatrix[Q] . (SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[Z] + SE3adMatrix[X] . SE3adMatrix[Z1] + SE3adMatrix[Z2]))];
];
Protect[H24Inv];

H32Inv[X_, Q_, Z_, eps_:SE3EPS] := Module[{x, u, Z1, Z2, Q1, Q2},
   x = Take[X, {1, 3}];
   Return[b2bar[x, eps]*(Transpose[{Join[x, {0, 0, 0}]}] . {Q}) . (SE3adMatrix[Z] . SE3adMatrix[X] - 2*SE3adMatrix[X] . SE3adMatrix[Z])];
];
Protect[H32Inv];

H34Inv[X_, Q_, Z_, eps_:SE3EPS] := Module[{x, u, Z1, Z2, Z3},
   x = Take[X, {1, 3}];
   Z1 = P[X, 1] . Z;
   Z2 = P[X, 2] . Z;
   Z3 = P[X, 3] . Z;
   Return[(-b4bar[x, eps])*(Transpose[{Join[x, {0, 0, 0}]}] . {Q}) . (SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[Z] +
   SE3adMatrix[X] . SE3adMatrix[Z2] + SE3adMatrix[X] . SE3adMatrix[X] . SE3adMatrix[Z1] + SE3adMatrix[Z3])];
];
Protect[H34Inv];



(* Higher-Order Approximations of dexp
Reference: [1] eq. (18) *)

(* Order 1 *)
SE3dexpOrder1[X_] := IdentityMatrix[6] + (1/2)*SE3adMatrix[X]; 
Protect[SE3dexpOrder1];

(* Order 2 *)
SE3dexpOrder2[X_] := IdentityMatrix[6] + (1/2)*SE3adMatrix[X] + (1/6)*MatrixPower[SE3adMatrix[X], 2]; 
Protect[SE3dexpOrder2];

(* Order 3 *)
SE3dexpOrder3[X_] := IdentityMatrix[6] + (1/2)*SE3adMatrix[X] + (1/6)*MatrixPower[SE3adMatrix[X], 2] +
   (1/24)*MatrixPower[SE3adMatrix[X], 3]; 
Protect[SE3dexpOrder3];

(* Order 4 *)
SE3dexpOrder4[X_] := IdentityMatrix[6] + (1/2)*SE3adMatrix[X] + (1/6)*MatrixPower[SE3adMatrix[X], 2] + 
   (1/24)*MatrixPower[SE3adMatrix[X], 3] + (1/120)*MatrixPower[SE3adMatrix[X], 4]; 
Protect[SE3dexpOrder4];



(* Higher-Order Approximations of dexp^-1
Reference: [1] eq. (18) *)

(* Order 1 *)
SE3dexpInvOrder1[X_] := IdentityMatrix[6] - (1/2)*SE3adMatrix[X]; 
Protect[SE3dexpInvOrder1];

(* Order 2 *)
SE3dexpInvOrder2[X_] := IdentityMatrix[6] - (1/2)*SE3adMatrix[X] + (1/6)*(1/2!)*MatrixPower[SE3adMatrix[X], 2]; 
Protect[SE3dexpInvOrder2];

(* Order 3 *)
SE3dexpInvOrder3[X_] := IdentityMatrix[6] - (1/2)*SE3adMatrix[X] + (1/6)*(1/2!)*MatrixPower[SE3adMatrix[X], 2]; 
Protect[SE3dexpInvOrder3];

(* Order 4 *)
SE3dexpInvOrder4[X_] := IdentityMatrix[6] - (1/2)*SE3adMatrix[X] + (1/6)*(1/2!)*MatrixPower[SE3adMatrix[X], 2] - 
    (1/30)*(1/4!)*MatrixPower[SE3adMatrix[X], 4]; 
Protect[SE3dexpInvOrder4];



(* Higher-Order Approximations of Ddexp
Reference: [1] Table 2 *)

(* Order 0 *)
SE3DdexpOrder0[X_, U_] := (1/2)*SE3adMatrix[U]; 
Protect[SE3DdexpOrder0];

(* Order 1 *)
SE3DdexpOrder1[X_, U_] := Module[{adU = SE3adMatrix[U], adX = SE3adMatrix[X]}, 
   Return[(1/2)*adU + (1/6)*(adU . adX + adX . adU)];
];
Protect[SE3DdexpOrder1];

(* Order 2 *)
SE3DdexpOrder2[X_, U_] := Module[{adU = SE3adMatrix[U], adX = SE3adMatrix[X], P1, P2},
   P1 = SE3adMatrix[X]; 
   P2 = SE3adMatrix[X] . P1; 
   Return[(1/2)*adU + (1/6)*(adU . P1 + P1 . adU) + (1/24)*(P2 . adU + adU . P2 + P1 . adU . P1)];
];
Protect[SE3DdexpOrder2];

(* Order 3 *)
SE3DdexpOrder3[X_, U_] := Module[{adU = SE3adMatrix[U], P1, P2, P3},
   P1 = SE3adMatrix[X];
   P2 = SE3adMatrix[X] . P1;
   P3 = SE3adMatrix[X] . P2;
   Return[(1/2)*adU + (1/6)*(adU . P1 + P1 . adU) +
   (1/24)*(P2 . adU + adU . P2 + P1 . adU . P1) +
   (1/120)*(P3 . adU + adU . P3 + P1 . adU . P2 +
   P2 . adU . P1)];
];
Protect[SE3DdexpOrder3];



(* Higher-Order Approximations of Ddexp^-1
Reference: [1] Table 3 *)

(* Order 0 *)
SE3DdexpInvOrder0[X_, U_] := (-2^(-1))*SE3adMatrix[U]; 
Protect[SE3DdexpInvOrder0];

(* Order 1 *)
SE3DdexpInvOrder1[X_, U_] := Module[{adU = SE3adMatrix[U], adX = SE3adMatrix[X]}, 
   Return[(-1/2)*adU + (1/12)*(adU . adX + adX . adU)];
];
Protect[SE3DdexpInvOrder1];

(* Order 2 *)
SE3DdexpInvOrder2[X_, U_] := SE3DdexpInvOrder1[X, U]
Protect[SE3DdexpInvOrder2];

(* Order 3 *)
SE3DdexpInvOrder3[X_, U_] := Module[{adU = SE3adMatrix[U], adX = SE3adMatrix[X], P1, P2, P3}, 
   P1 = SE3adMatrix[X]; 
   P2 = SE3adMatrix[X] . P1;
   P3 = SE3adMatrix[X] . P2; 
   Return[(-1/2)*adU + (1/12)*(adU . P1 + P1 . adU) - (1/720)*(P3 . adU + adU . P3 + adX . adU . P2 + P2 . adU . adX)];
];
Protect[SE3DdexpInvOrder3];



(* Higher-Order Approximations of D2dexp
Reference: [1] Table 4 *)

(* Order 0 *)
SE3D2dexpOrder0[X_, U_, S_] := (1/6)*SE3adMatrix[S] . SE3adMatrix[U] + (1/6)*SE3adMatrix[U] . SE3adMatrix[S]; 
Protect[SE3D2dexpOrder0];


(* Order 1 *)
SE3D2dexpOrder1[X_, U_, S_] := Module[{adU = SE3adMatrix[U], adS = SE3adMatrix[S], P1}, 
   P1 = SE3adMatrix[X]; Return[(1/6)*SE3adMatrix[S] . SE3adMatrix[U] + (1/6)*SE3adMatrix[U] . SE3adMatrix[S] + 
      (1/24)*(P1 . adS + adS . P1) . adU + (1/24)*(adU . adS + adS . adU) . P1 + (1/24)*(P1 . adU + adU . P1) . adS]]
Protect[SE3D2dexpOrder1];

(* Order 2 *)
SE3D2dexpOrder2[X_, U_, S_] := Module[{adU = SE3adMatrix[U], adS = SE3adMatrix[S], P1, P2}, 
   P1 = SE3adMatrix[X]; P2 = SE3adMatrix[X] . P1; Return[SE3D2dexpOrder1[X, U, S] + 
      (1/120)*adS . (P2 . adU + adU . P2 + P1 . adU . P1) + (1/120)*adU . (P2 . adS + adS . P2 + P1 . adS . P1) + 
      (1/120)*P2 . (adS . adU + adU . adS) + (1/120)*P1 . (adS . (P1 . adU + adU . P1) + adU . (P1 . adS + adS . P1))]]
Protect[SE3D2dexpOrder2];



(* Higher-Order Approximations of D2dexp^-1
Reference: [1] Table 5 *)

(* Order 0 *)
SE3D2dexpInvOrder0[X_, U_, S_] := (1/12)*SE3adMatrix[S] . SE3adMatrix[U] + (1/12)*SE3adMatrix[U] . SE3adMatrix[S]; 
Protect[SE3D2dexpInvOrder0];

(* Order 2 *)
SE3D2dexpInvOrder2[X_, U_, S_] := Module[{adU = SE3adMatrix[U], adS = SE3adMatrix[S], P1, P2}, 
   P1 = SE3adMatrix[X]; P2 = SE3adMatrix[X] . P1; Return[SE3D2dexpInvOrder0[X, U, S] - 
      (1/720)*adS . (P2 . adU + adU . P2 + P1 . adU . P1) - (1/720)*adU . (P2 . adS + adS . P2 + P1 . adS . P1) - 
      (1/720)*P2 . (adS . adU + adU . adS) - (1/720)*P1 . (adS . (P1 . adU + adU . P1) + adU . (P1 . adS + adS . P1))]]
Protect[SE3D2dexpInvOrder2];



(* Higher-Order Approximations of Jacobian of dexp.Z
Reference: [1] Table 2 *)

(* Order 0 *)
SE3dexpJacOrder0[X_, Z_] := (-2^(-1))*SE3adMatrix[Z]; 
Protect[SE3dexpJacOrder0];

(* Order 1 *)
SE3dexpJacOrder1[X_, Z_] := Module[{Z1}, Z1 = SE3adMatrix[X] . Z; 
    Return[(-2^(-1))*SE3adMatrix[Z] - (1/6)*(SE3adMatrix[Z1] + SE3adMatrix[X] . SE3adMatrix[Z])];
];
Protect[SE3dexpJacOrder1];

(* Order 2 *)
SE3dexpJacOrder2[X_, Z_] := Module[{Z1, Z2, P1, P2}, P1 = SE3adMatrix[X]; P2 = SE3adMatrix[X] . P1; Z1 = P1 . Z; 
    Z2 = P2 . Z; Return[(-2^(-1))*SE3adMatrix[Z] - (1/6)*(SE3adMatrix[Z1] + P1 . SE3adMatrix[Z]) - 
      (1/24)*(SE3adMatrix[Z2] + P1 . SE3adMatrix[Z1] + P2 . SE3adMatrix[Z])]; ]
Protect[SE3dexpJacOrder2];



(* Higher-Order Approximations of Hessian of Q^T.dexp.Z
Reference: [1] Table 6 *)

(* Order 0 *)
SE3dexpHessianOrder0[X_, Q_, Z_] := Module[{H2bar}, H2bar = (1/6)*SE3adBarMatrix[Q] . SE3adMatrix[Z]; 
    Return[H2bar + Transpose[H2bar]]; ]
Protect[SE3dexpHessianOrder0];

(* Order 1 *)
SE3dexpHessianOrder1[X_, Q_, Z_] := Module[{H3bar, P1, adZ, adZ1, Q1, Z1}, 
   P1 = SE3adMatrix[X];
   adZ = SE3adMatrix[Z];
   Z1 = P1 . Z;
   adZ1 = SE3adMatrix[Z1];
   Q1 = Transpose[P1] . Q; 
   H3bar = SE3adBarMatrix[Q] . (P1 . adZ + adZ1) + SE3adBarMatrix[Q1] . adZ;
   Return[SE3dexpHessianOrder0[X, Q, Z] + (1/24)*(H3bar + Transpose[H3bar])];
];
Protect[SE3dexpHessianOrder1];

(* Order 2 *)
SE3dexpHessianOrder2[X_, Q_, Z_] := Module[{H4bar, P1, P2, adZ, adZ1, adZ2, Q1, Q2, Z1, Z2}, 
   P1 = SE3adMatrix[X]; 
   P2 = SE3adMatrix[X] . P1; 
   adZ = SE3adMatrix[Z];
   Z1 = P1 . Z;
   Z2 = P2 . Z;
   adZ1 = SE3adMatrix[Z1];
   adZ2 = SE3adMatrix[Z2];
   Q1 = Transpose[P1] . Q;
   Q2 = Transpose[P2] . Q;
   H4bar = SE3adBarMatrix[Q] . (P2 . adZ + P1 . adZ1 + adZ2) + SE3adBarMatrix[Q1] . (P1 . adZ + adZ1) + SE3adBarMatrix[Q2] . adZ;
   Return[SE3dexpHessianOrder1[X, Q, Z] + (1/120)*(H4bar + Transpose[H4bar])];
];
Protect[SE3dexpHessianOrder2];

SE3dexpHessianSwitchOrder2[X_, Q_, Z_, eps_:SE3EPS] := Module[{x, n2, n}, x = X[[1 ;; 3]]; n2 = x . x; n = Sqrt[n2]; 
    If[n < eps, Return[SE3dexpHessianOrder2[X, Q, Z]]]; 
    Return[H11[X, Q, Z, eps] + H12[X, Q, Z, eps] + H13[X, Q, Z, eps] + 
     H14[X, Q, Z, eps] + H22[X, Q, Z, eps] + Transpose[H22[X, Q, Z, eps]] + H23[X, Q, Z, eps] + 
     Transpose[H23[X, Q, Z, eps]] + H24[X, Q, Z, eps] + Transpose[H24[X, Q, Z, eps]] + H31[X, Q, Z, eps] + 
     Transpose[H31[X, Q, Z, eps]] + H32[X, Q, Z, eps] + Transpose[H32[X, Q, Z, eps]] + H33[X, Q, Z, eps] + 
     Transpose[H33[X, Q, Z, eps]] + H34[X, Q, Z, eps] + Transpose[H34[X, Q, Z, eps]]];
];
Protect[SE3dexpHessianSwitchOrder2];

(* Close up the open environments *)
End[];
EndPackage[];
