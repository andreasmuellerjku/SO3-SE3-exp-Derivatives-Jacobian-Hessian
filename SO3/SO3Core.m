(* ::Package:: *)

(* 
Author: Andreas Mueller
Created: 2 January 2005
Last updated: 1 April 2026
Purpose: Core functions for manipulating on SO(3)
References: 
*)

BeginPackage["SO3Core`"];


(* =========================
   Basic Lie algebra helpers
   ========================= *)

R3Toso3::usage = "R3Toso3[x] returns the 3×3 skew-symmetric matrix for a 3-vector x.";
so3ToR3::usage = "so3ToR3[w] returns the axis vector corresponding to a 3×3 skew-symmetric matrix w (inverse of R3Toso3).";
SO3EPS::usage = "Default threshold for regarding a value to be zero.";

(* -------- SyntaxInformation for public API -------- *)

SyntaxInformation[R3Toso3] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[so3ToR3] = {"ArgumentsPattern" -> {_}};


(* Begin private section of the package *)
Begin["`Private`"];

(* Default threshold for regarding a value to be zero *)
SO3EPS = 1*^-12;

(* R3Toso3 — vector to skew-symmetric matrix
Input: x in R³
Output: x~ in so(3) such that x~ y = x × y
Notes: Implements the standard "tilde" operator. *)
Unprotect[R3Toso3]; 
R3Toso3[x_] := {{0, -x[[3]], x[[2]]}, {x[[3]], 0, -x[[1]]}, {-x[[2]], x[[1]], 0}}
Protect[R3Toso3]; 

(* so3ToR3 — skew-symmetric matrix to vector
Input: w in R^{3×3}, skew-symmetric
Output: x in R³ such that w = x~ (inverse of R3Toso3)
Notes: No validation that w is skew is performed. *)
Unprotect[so3ToR3]; 
so3ToR3[w_] :=  {w[[3, 2]], w[[1, 3]], w[[2, 1]]}
Protect[so3ToR3];

(* Close up the open environments *)
End[];
EndPackage[];
