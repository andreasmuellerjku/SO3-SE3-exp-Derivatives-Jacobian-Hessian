(* ::Package:: *)

(* 
Author: Andreas Mueller
Created: 2 January 2005
Last updated: 1 April 2026
Purpose: Core functions for manipulating on se(3) and SE(3)
References: 
*)

BeginPackage["SE3Core`","SO3Core`"];

Unprotect[se3ToTwist,TwistTose3,SE3adMatrix,RotPosToSE3,SE3Inverse]; 

(* ============================
   general operations on SE(3) 
   ============================ *)

TwistTose3::usage = "TwistTose3[X] maps a 6-vector twist X = {w,v} to its 4×4 matrix representation.";
se3ToTwist::usage = "se3ToTwist[Xhat] maps a 4×4 se(3) matrix to the 6-vector twist {w,v}.";
SE3adMatrix::usage = "SE3adMatrix[X] returns the 6×6 ad_X matrix for twist X (adjoint in se(3)).";
RotPosToSE3::usage = "RotPosToSE3[R, p] assembles rotation R and translation p into a 4×4 SE(3) matrix.";
SE3Inverse::usage = "SE3Inverse[H] returns the inverse of an SE(3) matrix H.";
SE3EPS::usage = "Default threshold for regarding a value to be zero.";

(* -------- SyntaxInformation for public API -------- *)

SyntaxInformation[TwistTose3] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[se3ToTwist] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[SE3adMatrix] = {"ArgumentsPattern" -> {_}};
SyntaxInformation[RotPosToSE3] = {"ArgumentsPattern" -> {_, _}};
SyntaxInformation[SE3Inverse] = {"ArgumentsPattern" -> {_}};

(* Begin private section of the package *)
Begin["`Private`"];

(* Default threshold for regarding a value to be zero *)
SE3EPS = 1*^-12;

(* se3ToTwist - matrix to twist vector
Input: X in se(3)
Output: {x, y} in R^6. *)
se3ToTwist[X_] := Module[{x, y}, x = so3ToAxis[Take[X, {1, 3}, {1, 3}]]; y = Flatten[Take[X, {1, 3}, {4, 4}]]; 
   Return[Join[x, y]];
]; 
Protect[se3ToTwist];

(* TwistTose3  twist vector to 6x6 se(3) matrix (hat)
Input: X = {w, v} in R^6
Output: X^ in se(3) *)
TwistTose3[X_] := Module[{w = Take[X, 3]}, 
       Return[ArrayFlatten[{{{{0, -w[[3]], w[[2]]}, {w[[3]], 0, -w[[1]]}, {-w[[2]], w[[1]], 0}}, Transpose[{Take[X, -3]}]},{{{0, 0, 0}}, {{0}}}}]]]; 
Protect[TwistTose3]; 

(* SE3adMatrix - matrix representation of adjoint on se(3), encodes Lie bracket in 6-vectors
Input: X = {x, y} in R^6
Output: ad_X such that ad_X Y = [X, Y]
*)
SE3adMatrix[X_] := Module[{adMat, x = Take[X, {1, 3}], y = Take[X, {4, 6}]}, 
   adMat = ArrayFlatten[{{AxisToso3[x], 0*IdentityMatrix[3]}, {AxisToso3[y], AxisToso3[x]}}];
   Return[adMat];
]; 
Protect[SE3adMatrix]; 

(* Assemble a rotation matrix R and translation vector p to a SE(3) matrix *)
RotPosToSE3[R_, p_] := ArrayFlatten[{{R, Transpose[{p}]}, {{{0, 0, 0}}, {{1}}}}]
Protect[RotPosToSE3]; 

(* SE3Inverse - inverse of H in SE(3) *)
SE3Inverse[H_] := Module[{R = H[[1 ;; 3, 1 ;; 3]], p = H[[1 ;; 3, 4]]}, 
       Return[RotPosToSE3[Transpose[R], -Transpose[R] . p]];
]; 
Protect[SE3Inverse];

(* Close up the open environments *)
End[];
EndPackage[];
