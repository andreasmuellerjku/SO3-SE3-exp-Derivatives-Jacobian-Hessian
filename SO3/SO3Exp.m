(* ::Package:: *)

(* 
Author: Andreas Mueller
Created: 2 January 2005
Last updated: 1 April 2026
Purpose: Computing on SO(3) in terms of canonical coordinates
References: 
[1] A. Mueller: 
[2] A. Mueller: Review of the exponential and Cayley map on SE(3) as relevant for Lie group integration of the generalized Poisson equation and flexible multibody systems,  Proc. Royal Soc. A, September 2021; 477 (2253): 20210303. 
https://doi.org/10.1098/rspa.2021.0303
https://arxiv.org/abs/2303.07928
*)

BeginPackage["SO3Exp`","SO3Core`"];

Unprotect[SO3Exp, SO3Log, SO3dexp, SO3Ddexp, SO3dexpInv, SO3D2dexp, SO3DdexpInv, SO3D2dexpInv];

(* =========================
   SO(3) dexp and derivatives
   ========================= *)

SO3Exp::usage = "SO3Exp[x] returns the SO(3) rotation matrix exp(x) using the sinc-based closed form.";
SO3Log::usage = "SO3Log[R, eps] returns the logarithm vector x such that R = exp(x). Small-angle threshold eps.";
SO3dexp::usage = "SO3dexp[x, eps] returns the 3×3 dexp at x in so(3).";
SO3Ddexp::usage = "SO3Ddexp[x, y, eps] returns the linear map D(dexp)_x[y] on so(3).";
SO3dexpInv::usage = "SO3dexpInv[x, eps] returns the inverse of dexp at x in so(3).";
SO3D2dexp::usage = "SO3D2dexp[x, y, u, v, eps] returns the second Fréchet derivative of dexp at x acting on (y,u,v) (RSPA eq. 2.49).";
SO3DdexpInv::usage = "SO3DdexpInv[x, y, eps] returns the derivative of dexp^{-1} at x applied to y (RSPA eq. 2.40).";
SO3D2dexpInv::usage = "SO3D2dexpInv[x, y, u, v, eps] returns the second derivative of dexp^{-1} at x (RSPA eq. 2.53).";

(* -------- SyntaxInformation for public API -------- *)

SyntaxInformation[SO3Exp] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[SO3Log] = {"ArgumentsPattern" -> {_, _.}};
SyntaxInformation[SO3dexp] = {"ArgumentsPattern" -> {_, _.}};
SyntaxInformation[SO3Ddexp] = {"ArgumentsPattern" -> {_, _, _.}};
SyntaxInformation[SO3dexpInv] = {"ArgumentsPattern" -> {_, _.}};
SyntaxInformation[SO3D2dexp] = {"ArgumentsPattern" -> {_, _, _, _, _.}};
SyntaxInformation[SO3DdexpInv] = {"ArgumentsPattern" -> {_, _, _.}};
SyntaxInformation[SO3D2dexpInv] = {"ArgumentsPattern" -> {_, _, _, _, _.}};

(* Begin private section of the package *)
Begin["`Private`"];

(* SO3Exp — exponential map on SO(3)
Input: x in R³ (scaled rotation vector)
Output: R = exp(x~) in SO(3)
Well-conditioned for very small-angle through sinc; no explicit branching.
Reference: [2] eq. 2.8 *)
SO3Exp[x_?VectorQ] := IdentityMatrix[3] + Sinc[Sqrt[x . x]]*R3Toso3[x] + 
       (Sinc[Sqrt[x . x]/2]^2/2)*R3Toso3[x] . R3Toso3[x]
Protect[SO3Exp];
 
(* SO3Log — logarithm on SO(3)
Input: R in SO(3), eps > 0 (threshold)
Output: x in R³ so that R = exp(x~)
Small-angle: If ||x|| < eps, returns {0,0,0}. *)
SO3Log[R_, eps_:SO3EPS] := Module[{t}, 
	t = ArcCos[(R[[1,1]] + R[[2,2]] + R[[3,3]] - 1)/2]; 
	If[t < eps, Return[{0, 0, 0}]]; 
    Return[t*(so3ToR3[R - Transpose[R]]/(2*Sin[t]))]; ]
Protect[SO3Log]; 

(* SO3dexp - dexp on so(3)
Input: x in R³, eps > 0
Output: 3×3 matrix dexp_x
Reference: [2] eq. 2.13 *)
SO3dexp[x_,eps_:SO3EPS] := Module[{xtilde, nxtilde, n, n2, Alpha, Beta},
	n2 = x . x; 
	n = Sqrt[n2]; 
	If[n  <eps, Return[IdentityMatrix[3]]];
	Alpha = Sinc[n];
	Beta = Sinc[n/2]^2;
	xtilde = R3Toso3[x];
	nxtilde = R3Toso3[x/n];
	Return[IdentityMatrix[3] + 1/2*Beta*xtilde + (1 - Alpha)*nxtilde . nxtilde];
	];
Protect[SO3dexp]; 

(* SO3Ddexp - directional derivative of dexp (1st Fréchet derivative)
Input: x, y in R³, eps > 0
Output: 3×3 matrix D(dexp)_x(y) 
Reference: [2] eq. 2.33 *)
SO3Ddexp[x_,y_,eps_:SO3EPS]:=Module[{n,n2, nx, xtilde, nxtilde, ytilde, Alpha, Beta, Delta},
  n2 = x . x;
  n = Sqrt[n2];
  ytilde = R3Toso3[y];
  If[n<eps, Return[1/2*R3Toso3[y]]];
  Alpha =Sinc[n];
  Beta =(Sinc[n/2])^2;
  Delta =(1-Alpha)/n2;
  xtilde = R3Toso3[x];
  nx = x/n;
  nxtilde = R3Toso3[nx];
  Return[1/2*Beta*ytilde + Delta*(xtilde . ytilde + ytilde . xtilde) +
  (x . y)*(Alpha-Beta)/n*nxtilde+ 
  (x . y)*(Beta/2-3*Delta)*nxtilde . nxtilde];
];
Protect[SO3Ddexp];

(* SO3dexpInv - inverse of dexp on so(3)
Input: x in R³, eps > 0
Output: 3×3 matrix (dexp_x)^{-1}
Reference: [2] eq. 2.16 *)
SO3dexpInv[x_,eps_:SO3EPS] := Module[{xtilde, nxtilde, n, n2, Alpha, Beta, Gamma},
	n2 = x . x; 
	n = Sqrt[n2]; 
	If[n<eps,Return[IdentityMatrix[3]]];
	Alpha = Sinc[n];
	Beta = Sinc[n/2]^2; 
	Gamma = Alpha/Beta; 
	xtilde = R3Toso3[x];
	nxtilde = R3Toso3[x/n];
	Return[IdentityMatrix[3] - (1/2)*xtilde + (1 - Gamma)*nxtilde . nxtilde]
	];
Protect[SO3dexpInv];

(* SO3D2dexp - second Fréchet derivative of dexp.
Input: x, y, u, v in R³, eps > 0
Output: 3×3 matrix (D^2dexp_x)(y,u,v). 
Reference: [2] eq. 2.49 *)
SO3D2dexp[x_, y_, u_, v_, eps_:SO3EPS] := Module[{n, n2, Beta, Alpha, Gamma, Delta, xtilde, ytilde, utilde, vtilde},
   n2 = x . x;
   n = Sqrt[n2]; 
   ytilde = R3Toso3[y];
   utilde = R3Toso3[u];
   vtilde = R3Toso3[v];
   If[n < eps, Return[(1/2)*vtilde + (1/6)*(ytilde . utilde + utilde . ytilde)]]; 
   xtilde = R3Toso3[x];
   Alpha = Sinc[n];
   Beta = Sinc[n/2]^2;
   Gamma = Alpha/Beta; 
   Delta = (1 - Alpha)/n2; 
   Return[(1/2)*Beta*vtilde + x . u*(Alpha - Beta)*(ytilde/n2) + ((Beta/2 - 3*Delta)/n2)*(x . u*(ytilde . xtilde + xtilde . ytilde) 
   + x . y*(utilde . xtilde + xtilde . utilde)) 
   + Delta*(xtilde . vtilde + vtilde . xtilde + utilde . ytilde + ytilde . utilde) 
   + (1/n2)*((Alpha - Beta)*(x . v + u . y) - (1/2)*Beta*(x . u*x . y) 
   + ((1 - 5*Alpha + 4*Beta)/n2)*x . y*x . u)*xtilde 
   + (1/n2)*((Beta/2 - 3*Delta)*(x . v + y . u) + (1/n2)*(Alpha - (7/2)*Beta + 15*Delta)*x . y*x . u)*xtilde . xtilde + ((Alpha - Beta)/n2)*x . y*utilde];
];
Protect[SO3D2dexp];

(* SO3DdexpInv - derivative of dexp^{-1}.
Input: x, y in R³, eps > 0
Output: 3×3 matrix D(dexp^{-1})_x(y). 
Reference: [2] eq. 2.40 *)
SO3DdexpInv[x_, y_, eps_:SO3EPS] := Module[{n, n2, Alpha, Beta, Delta, Gamma, xtilde, ytilde}, 
	n2 = x . x; 
	n = Sqrt[n2];
	ytilde = R3Toso3[y];
    If[n < eps, Return[(-1/2)*ytilde]]; 
    xtilde = R3Toso3[x];
    Alpha = Sinc[n]; 
    Beta = Sinc[n/2]^2; 
    Gamma = Alpha/Beta; 
    Delta = (1 - Alpha)/n2; Return[(-1/2)*ytilde + ((1 - Gamma)/n2)*(xtilde . ytilde + 
        ytilde . xtilde) + (x . y/n2^2)*(1/Beta + Gamma - 2)*xtilde . xtilde]; 
    ];
Protect[SO3DdexpInv];

(* SO3D2dexpInv - second derivative of dexp^{-1}.
Input: x, y, u, v in R³, eps > 0
Output: 3×3 matrix (D^2dexp^{-1}_x)(y,u,v).
Reference: [2] eq. 2.53 *)
SO3D2dexpInv[x_, y_, u_, v_, eps_:SO3EPS] := Module[{n, n2, Beta, Alpha, Gamma, Delta, xtilde, ytilde, utilde, vtilde}, 
   n2 = x . x;
   n = Sqrt[n2];
   ytilde = R3Toso3[y];
   utilde = R3Toso3[u];
   vtilde = R3Toso3[v];
   If[n < eps, Return[(-1/2)*vtilde + (1/12)*(utilde . ytilde + ytilde . utilde)]];
   Alpha = Sinc[n];
   Beta = Sinc[n/2]^2;
   Gamma = Alpha/Beta;
   Delta = (1 - Alpha)/n2; 
   xtilde = R3Toso3[x];
   Return[-1/2*vtilde+1/n2*(1/4*(x . u)*(xtilde . ytilde + ytilde . xtilde) +
   (1-Gamma)*(utilde . ytilde+ytilde . utilde +
   vtilde . xtilde+xtilde . vtilde)) +
   1/n2^2*((x . u)*(Gamma-2+Gamma^2)*(xtilde . ytilde+ytilde . xtilde) +
   (1/Beta+Gamma-2)*(x . y)*(xtilde . utilde+utilde . xtilde) +
   ((1/Beta+Gamma-2)*(u . y+x . v)-1/4*(x . y)*(x . u))*xtilde . xtilde) +
   1/n2^3*(x . y)*(x . u)*(8-Gamma*(3+Gamma)-2/Beta*(1+Gamma))*xtilde . xtilde];
];
Protect[SO3D2dexpInv];

(* Close up the open environments *)
End[];
EndPackage[];
