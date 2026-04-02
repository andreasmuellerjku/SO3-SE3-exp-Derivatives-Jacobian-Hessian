(* ::Package:: *)

(* 
Author: Andreas Mueller
Created: 2 January 2005
Last updated: 1 April 2026
Purpose: Computing on SE(3) in terms of canonical coordinates using the block partitioning according to the semi-direct product structure SE(3) = SO(3) |x R³
References: 
[1] A. Mueller: Closed Form Relations and Higher-Order Approximations of First and Second Derivatives of the Tangent Operator on SE(3), ZAMM - Journal of Applied Mathematics and Mechanics / Zeitschrift für Angewandte Mathematik und Mechanik, 2026, in press
[2] A. Mueller: Review of the exponential and Cayley map on SE(3) as relevant for Lie group integration of the generalized Poisson equation and flexible multibody systems,  Proc. Royal Soc. A, September 2021; 477 (2253): 20210303. 
https://doi.org/10.1098/rspa.2021.0303
https://arxiv.org/abs/2303.07928
*)

BeginPackage["SE3PExp`","SE3Core`","SO3Core`","SO3Exp`"];

Unprotect[SE3PExp,SE3PLog,SE3Pdexp,SE3PDdexp,SE3PdexpInv,SE3PDdexpInv]; 

SE3PExp::usage = "SE3PExp[X, eps] returns the SE(3) matrix exp(X) for twist X, with small-angle threshold eps.";
SE3PLog::usage = "SE3PLog[H, eps] returns the twist X such that H = exp[X] (matrix log mapped to twist). Small-angle threshold eps.";

(* ==============================================
   SE(3) using the SO(3)|xR³ block partitioning
   exp, log, dexp, its inverse, and derivatives
   ============================================== *)

SE3Pdexp::usage = "SE3Pdexp[X, eps] returns the 6×6 dexp at twist X.";
SE3PdexpInv::usage = "SE3PdexpInv[X, eps] returns the inverse of dexp at twist X.";
SE3PDdexp::usage = "SE3PDdexp[X, U, eps] returns the directional derivative D_X(dexp)(U) of dexp at X along U.";
SE3PDdexpInv::usage = "SE3PDdexpInv[X, U, eps] returns the directional derivative D_X(dexp^-1)(U) of dexp at X along U.";

(* -------- SyntaxInformation for public API -------- *)

SyntaxInformation[SE3PExp] = {"ArgumentsPattern" -> {_, _.}};
SyntaxInformation[SE3PLog] = {"ArgumentsPattern" -> {_, _.}};

SyntaxInformation[SE3Pdexp] = {"ArgumentsPattern" -> {_, _.}};
SyntaxInformation[SE3PdexpInv] = {"ArgumentsPattern" -> {_, _.}};
SyntaxInformation[SE3PDdexp] = {"ArgumentsPattern" -> {_, _, _.}};
SyntaxInformation[SE3PDdexpInv] = {"ArgumentsPattern" -> {_, _, _.}};


(* Begin private section of the package *)
Begin["`Private`"];

 
(* SE3Exp - exponential map on SE(3)
Input: X = {x, y} in R^6, eps > 0
Output: H = Exp(X^) in SE(3)
Method: Computes SO(3) and R³ entries separately 
Reference: [1] eq. (A.1), [2] eq. (2.27) *)
SE3PExp[X_, eps_:SE3EPS] := Module[{x, y, R, p, n, a, b, nx}, 
   x = Take[X, {1, 3}];
   y = Take[X, {4, 6}];
   n = Sqrt[x . x];
   If[n < eps, Return[RotPosToSE3[IdentityMatrix[3], y]]];
   a = Sinc[n];
   b = Sinc[n/2]^2; 
   nx = x/n;
   R = IdentityMatrix[3] + a*R3Toso3[x] + (b/2)*R3Toso3[x] . R3Toso3[x];
   p = ((b/2)*R3Toso3[x] - a*R3Toso3[nx] . R3Toso3[nx] + Transpose[{nx}] . {nx}) . y;
   Return[RotPosToSE3[R, p]];   
]; 
Protect[SE3PExp]; 

(* SE3Log - logarithm on SE(3)
Input: H in SE(3), eps > 0
Output: X = {x, y} in R^6 with H = exp(X^) *)
SE3PLog[H_,eps_:SE3EPS] := Module[{R, p, y, x}, 
   R = H[[1 ;; 3,1 ;; 3]]; 
   p = H[[1 ;; 3,4]]; 
   x = SO3Log[R,eps]; 
   y = SO3dexpInv[x,eps] . p; TwistTose3[Join[x, y]]]; 
Protect[SE3PLog];


(* SE3Pdexp - dexp on se(3)
Input: X = {x, y} in R^6, eps > 0
Output: 6×6 block matrix dexp_X 
Reference: [1] eq. (11), [2]  eq. (2.29) *)
SE3Pdexp[X_, eps_:SE3EPS] := Module[{n, x, y},
   x = Take[X, 3];
   y = Take[X, -3];
   n = Norm[x]; 
   If[n < eps, Return[IdentityMatrix[6]]];
   Return[Join[Join[SO3dexp[x], 0*IdentityMatrix[3], 2], Join[SO3Ddexp[x, y, eps], SO3dexp[x, eps], 2]]];
];
Protect[SE3Pdexp];


(* SE3PDdexp - derivative of dexp on se(3)
Input: X, U in R^6, eps > 0
Output: 6×6 matrix D(dexp)_X(U)
Reference: [1] eq. (A.5), [2]  eq. (2.47) *)
SE3PDdexp[X_, U_, eps_:SE3EPS] := Module[{n, x, y, u, v}, 
    x = Take[X, 3];
    y = Take[X, -3];
    u = Take[U, 3]; 
    v = Take[U, -3];
    n = Norm[x];
    If[n < eps, Return[(1/2)*SE3adMatrix[U]]]; 
    Return[ArrayFlatten[{{SO3Ddexp[x, u, eps], 0*IdentityMatrix[3]}, {SO3D2dexp[x, y, u, v, eps], SO3Ddexp[x, u, eps]}}]];
];
Protect[SE3PDdexp]; 

(* SE3PdexpInv - inverse dexp on se(3) (6×6)
Input: X = {w, v} in R^6, eps > 0
Output: 6×6 block matrix dexp_X^-1.
Reference: [1] eq. (11), [2] eq. (2.39) *)
SE3PdexpInv[X_, eps_:SE3EPS] := Module[{x, y}, 
   x = Take[X, {1, 3}];
   y = Take[X, {4, 6}]; 
   Return[Join[Join[SO3dexpInv[x, eps], 0*IdentityMatrix[3], 2], Join[SO3DdexpInv[x, y, eps], SO3dexpInv[x, eps], 2]]];
];
Protect[SE3PdexpInv]; 

(* SE3PDdexpInv - derivative of dexp^{-1} on se(3)
Input: X, U in R^6, eps > 0
Output: 6×6 matrix (D(dexp)_X)^1(U)
Reference: [1] eq. (A.6), [2] eq. (2.51) *)
SE3PDdexpInv[X_, U_, eps_:SE3EPS] := Module[{n, x, y, u, v},
   x = Take[X, {1, 3}];
   y = Take[X, {4, 6}]; 
   u = Take[U, {1, 3}];
   v = Take[U, {4, 6}];
   Return[ArrayFlatten[{{SO3DdexpInv[x, u, eps], 0*IdentityMatrix[3]}, {SO3D2dexpInv[x, y, u, v, eps], SO3DdexpInv[x, u, eps]}}]];
];
Protect[SE3PDdexpInv];

(* Close up the open environments *)
End[];
EndPackage[];
