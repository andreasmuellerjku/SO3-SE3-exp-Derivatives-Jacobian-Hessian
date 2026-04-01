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

AxisToso3::usage = "AxisToso3[x] returns the 3×3 skew-symmetric matrix for a 3-vector x.";
so3ToAxis::usage = "so3ToAxis[w] returns the axis vector corresponding to a 3×3 skew-symmetric matrix w (inverse of AxisToso3).";
SO3EPS::usage = "Default threshold for regarding a value to be zero.";

(* -------- SyntaxInformation for public API -------- *)

SyntaxInformation[AxisToso3] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[so3ToAxis] = {"ArgumentsPattern" -> {_}};


(* Begin private section of the package *)
Begin["`Private`"];

(* Default threshold for regarding a value to be zero *)
SO3EPS = 1*^-12;

(* AxisToso3 — vector to skew-symmetric matrix
Input: x in R³
Output: x~ in so(3) such that x~ y = x × y
Notes: Implements the standard "tilde" operator. *)
Unprotect[AxisToso3]; 
AxisToso3[x_] := {{0, -x[[3]], x[[2]]}, {x[[3]], 0, -x[[1]]}, {-x[[2]], x[[1]], 0}}
Protect[AxisToso3]; 

(* so3ToAxis — skew-symmetric matrix to vector
Input: w in R^{3×3}, skew-symmetric
Output: x in R³ such that w = x~ (inverse of AxisToso3)
Notes: No validation that w is skew is performed. *)
Unprotect[so3ToAxis]; 
so3ToAxis[w_] :=  {w[[3, 2]], w[[1, 3]], w[[2, 1]]}
Protect[so3ToAxis];

(* Close up the open environments *)
End[];
EndPackage[];
