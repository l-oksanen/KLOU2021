(* Symmetric outer product *)

otimes[u_, v_] := Outer[Times, u, v] + Outer[Times, v, u];

(* Minkowski metric and Q and P symbols *)

G = {{-1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
Q[xi_] := 1/(xi.G.xi);
P[g_, xi_] := xi.G.g.G.xi;

(* Outgoing direction *)
xi = {1, 1, 0, 0};

Print["Use the following Pythagorean quadruples as incoming \
directions"];
pythaquads = 
 Flatten[Table[{p^2 + m^2 + n^2, 2 m p, 2 n p, p^2 - (m^2 + n^2)}, {m,
     1, 3}, {n, 1, 2}, {p, 1, 1}], 2]
bss = Subsets[pythaquads, {4}];

Print["Trying " <> ToString[Length[bss]] <> 
  " choices of 4-tuples of incoming directions"]
symbs = {};
For[jj = 1, jj <= Length[bss], jj++,
 Print["\nChoice " <> ToString[jj]];
 
 (* The incoming directions and the change of basis matrix *)
 
 bs = bss[[jj]];
 B = Transpose[bs];
 Print["The incoming directions b_j as column vectors"];
 Print[MatrixForm[B]];
 
 (* Check that the incoming vectors form a basis *)
 
 If[Det[B] == 0, 
  Print["SKIPPING: the vectors b_j do not form a basis"]; 
  Continue[]];
 
 (* The coefficients c_j satisfying xi_j = c_j b_j *)
 
 cs = Inverse[B].xi;
 xis = DiagonalMatrix[cs].bs;
 
 (* Check that no xi_j is zero  *)
 
 If[MemberQ[Map[Norm, xis], 0], 
  Print["SKIPPING: one of the vectors xi_j is zero"]; Continue[]];
 
 (* Function to substitute xj, j=1,2,.., with the xi's given by the \
indices *)
 
 indlabels[es_] := MapThread[List, {Table[i, {i, Length[es]}], es}];
 subs[is_] := Map[Function[pair,
    Symbol["x" <> ToString[pair[[1]]]] -> xis[[pair[[2]]]]], 
   indlabels[Flatten[is]]];
 dosubs[e_, is_] := Activate[Inactivate[e] /. subs[is]];
 SetAttributes[dosubs, HoldFirst];
 
 (* The 2nd order interactions g_jk *)
 gjk[jk_] := dosubs[
   -2*Q[x1 + x2]*otimes[x1, x2]
   , jk];
 
 (* The 3rd order interactions phi_jkl *)
 
 phijkl[jkl_] := With[{g23 = gjk[jkl[[2]]]}, dosubs[
    Q[x1 + x2 + x3]*P[g23, x1]
    , jkl]];
 
 (* The 4th order type 1+3 interactions *)
 
 t13s = {{2, {1, {3, 4}}}, {1, {2, {3, 4}}}, {4, {3, {1, 
      2}}}, {3, {4, {1, 2}}}};
 gt13[is_] := With[{phi123 = phijkl[is[[2]]]}, dosubs[
    -2*phi123*otimes[x1, x2 + x3 + x4]
    , is]];
 
 (* The 4th order type 2+2 interactions from g(D,D) *)
 
 t22s = {{{1, 2}, {3, 4}}, {{3, 4}, {1, 2}}};
 gt22box[is_] := With[{g12 = gjk[is[[1]]], g34 = gjk[is[[2]]]}, dosubs[
    P[g12, x3 + x4]*g34
    , is]];
 
 (* Lowered Christoffel symbols *)
 
 gammalow[g_, x_, p_, k_, 
   q_] := (1/2)*(x[[p]]*g[[k, q]] + x[[q]]*g[[p, k]] - 
     x[[k]]*g[[p, q]]);
 
 (* The 4th order type 2+2 interactions from the Christoffel symbol \
producs *)
 
 gammacontraction[g1_, g2_, x1_, x2_] := 
  Table[Sum[G[[p, q]]*G[[r, s]]*(
      gammalow[g1, x1, p, r, j]*gammalow[g2, x2, q, s, k]
       + gammalow[g1, x1, p, r, j]*gammalow[g2, x2, q, k, s]
       + gammalow[g1, x1, p, r, k]*gammalow[g2, x2, q, j, s]
      ), {p, 1, 4}, {q, 1, 4}, {r, 1, 4}, {s, 1, 4}], {j, 1, 4}, {k, 
    1, 4}];
 gt22gamma = Simplify[
   2*(gammacontraction[gjk[{1, 2}], gjk[{3, 4}], xis[[1]] + xis[[2]], 
       xis[[3]] + xis[[4]]]
      + gammacontraction[gjk[{3, 4}], gjk[{1, 2}], 
       xis[[3]] + xis[[4]], xis[[1]] + xis[[2]]]
     )];
 
 full = gt22gamma + 
   Simplify@Apply[Plus, Join[ Map[gt13, t13s], Map[gt22box, t22s]]];
 
 Print["The full 4th order interaction"];
 Print[MatrixForm[full]];
 
 symbs = Join[
   symbs, {Join[full[[1]], full[[2, 2 ;; 4]], full[[3, 3 ;; 4]], 
     full[[4, 4 ;; 4]]]}];
 ]

Print["\nNumber of successful choices"]
Print[Length[symbs]]

Print["The dimension of the subspace spanned by the principal symbols"]
Print[MatrixRank[symbs]]

