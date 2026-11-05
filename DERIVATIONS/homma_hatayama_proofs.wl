(* ::Package:: *)
(* Proofs for the Homma & Hatayama (2014) collision coefficients.

   Run by hand after changing the derivation (about 3 minutes):
     math -script DERIVATIONS/homma_hatayama_proofs.wl
   It
     1. computes the Rosenbluth potentials of Homma's Eq. (1) directly and
        checks them against a direct 3D integral,
     2. proves that their derivatives are the printed Eqs. (2), (5), (6),
        prefactors included,
     3. proves that the scalar kernel functions of homma_hatayama_common.wl
        reproduce them, that the small-speed series has converged at the
        switch, and the Maxwellian average used for the B03/B07 oracle.
   homma_hatayama.wl then generates the Fortran include and oracle values
   from the proved kernel functions. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "homma_hatayama_common.wl"}]];
(* ------------------------------------------------------------------ *)
(* 1. Rosenbluth potentials of Eq. (1), Legendre components l = 0, 1  *)
(* ------------------------------------------------------------------ *)

f0[s_] := Exp[-s^2]/Pi^(3/2);
(* l = 1 part per unit |g|: f~ = f0 (2/5)(s^2 - 5/2) s cos(angle to g) *)
f1[s_] := Exp[-s^2]/Pi^(3/2) (2/5) (s^2 - 5/2) s;
pos = {v > 0};

(* 1/|u-u'| and |u-u'| expanded in Legendre polynomials; the angular
   integral of P_l over the sphere gives 4 pi/(2l+1) P_l(cos). *)
H0[v_] = Simplify[4 Pi (Integrate[s^2 f0[s], {s, 0, v}]/v
                        + Integrate[s f0[s], {s, v, Infinity}]), pos];
H1[v_] = Simplify[4 Pi/3 (Integrate[s^3 f1[s], {s, 0, v}]/v^2
                          + v Integrate[f1[s], {s, v, Infinity}]), pos];
G0[v_] = Simplify[4 Pi (Integrate[s^2 f0[s] (v + s^2/(3 v)), {s, 0, v}]
                        + Integrate[s^2 f0[s] (s + v^2/(3 s)), {s, v, Infinity}]), pos];
G1[v_] = Simplify[4 Pi/3 (Integrate[s^2 f1[s] (s^3/(5 v^2) - s), {s, 0, v}]
                          + Integrate[s^2 f1[s] (v^3/(5 s^2) - v), {s, v, Infinity}]), pos];

Print["H0 = ", InputForm[H0[v]]];
Print["H1 = ", InputForm[H1[v]]];
Print["G0 = ", InputForm[G0[v]]];
Print["G1 = ", InputForm[G1[v]]];

lap[g_, l_] := D[g, {v, 2}] + 2/v D[g, v] - l (l + 1)/v^2 g;
check[Simplify[lap[H0[v], 0] + 4 Pi f0[v], pos] === 0, "Laplace H0 = -4 pi f0"];
check[Simplify[lap[H1[v], 1] + 4 Pi f1[v], pos] === 0, "Laplace H1 = -4 pi f1"];
check[Simplify[lap[G0[v], 0] - 2 H0[v], pos] === 0, "Laplace G0 = 2 H0"];
check[Simplify[lap[G1[v], 1] - 2 H1[v], pos] === 0, "Laplace G1 = 2 H1"];
(* regular at the origin and no growth beyond the Maxwellian G ~ u *)
check[Limit[H1[v], v -> 0] === 0 && Limit[G1[v], v -> 0] === 0, "l = 1 parts vanish at u = 0"];
check[Limit[H1[v], v -> Infinity] === 0 && Limit[G1[v] v^2, v -> Infinity] === 1/10,
      "l = 1 parts decay at large u"];

(* Spot check of the Legendre route against a direct 3D integral of
   f~(u') |u - u'| for one point and one heat flux (no expansion used). *)
(* Spherical coordinates (rho, th, ph) centred at u, so |u - u'| = rho and
   the integrand is smooth. *)
Module[{u = {3/10, -7/10, 11/10}, g = {2/5, 1/5, -1/2}, direct, viaLegendre, rr, fb, ang},
  fb[w_?VectorQ] := Exp[-w.w]/Pi^(3/2) (1 + (2/5) (w.w - 5/2) g.w);
  ang[rho_?NumericQ] := NIntegrate[Sin[th] fb[u + rho {Sin[th] Cos[ph], Sin[th] Sin[ph], Cos[th]}],
    {th, 0, Pi}, {ph, 0, 2 Pi}, Method -> {"GaussKronrodRule", "Points" -> 15},
    WorkingPrecision -> 20, PrecisionGoal -> 13];
  direct = NIntegrate[rho^3 ang[rho], {rho, 0, 3/2, 4, 10}, WorkingPrecision -> 20, PrecisionGoal -> 12];
  rr = Norm[u];
  viaLegendre = N[G0[rr] + G1[rr] g.u/rr, 20];
  Print["G direct 3D = ", direct, ", via Legendre = ", viaLegendre];
  check[Abs[direct - viaLegendre] < 1*^-10, "Legendre potentials match direct 3D integral"]];

(* ------------------------------------------------------------------ *)
(* 2. Cartesian A, D and Homma's printed Eqs. (2), (5), (6)           *)
(* ------------------------------------------------------------------ *)

uu = {u1, u2, u3}; gg = {g1, g2, g3};
rr = Sqrt[uu.uu];
cartAssume = Element[uu | gg, Reals] && rr > 0;
Hc = H0[rr] + H1[rr] gg.uu/rr;
Gc = G0[rr] + G1[rr] gg.uu/rr;
aC = D[Hc, {uu}];
dC = D[Gc, {uu, 2}];

uo = Outer[Times, uu, uu];
gu = Outer[Times, gg, uu] + Outer[Times, uu, gg];

(* Homma Eq. (2): F = F0 + F_gradT; with A = F/m_a and the normalisation
   above, F0/m_a = P_A aM and F_gradT/m_a = P_A aT where (shown in
   section 3 below) *)
aMH = -(Erf[rr] - rr phiP[rr])/rr^3 uu;
aTH = -2/(5 Sqrt[Pi]) Exp[-rr^2] (gg - 2 gg.uu uu);
(* Homma Eq. (5) times 2 (P_D = Gamma n/(2 v_tb) vs Homma's n/(16 pi) L) *)
dMH = 1/2 (((2/rr - 1/rr^3) Erf[rr] + phiP[rr]/rr^2) id3
           + ((-2/rr^3 + 3/rr^5) Erf[rr] - 3/rr^4 phiP[rr]) uo);
(* Homma Eq. (6), with q expressed through g (factor 1/5, section 3) *)
dTH = 1/5 ((-3/(2 rr^5) Erf[rr] + (1/rr^2 + 3/(2 rr^4)) phiP[rr]) (gu + gg.uu id3)
           + (15/(2 rr^7) Erf[rr] + (-2/rr^2 - 5/rr^4 - 15/(2 rr^6)) phiP[rr]) gg.uu uo);

zeroQ[expr_] := AllTrue[Flatten[{Simplify[expr, cartAssume]}], # === 0 &];
gZero = Thread[gg -> 0];
check[zeroQ[(aC /. gZero) - aMH], "grad H (q = 0) = Homma Eq. (2) friction F0"];
check[zeroQ[aC - (aC /. gZero) - aTH], "grad H (part linear in q) = Homma Eq. (2) thermal force"];
check[zeroQ[(dC /. gZero) - dMH], "Hessian G (q = 0) = Homma Eq. (5)"];
check[zeroQ[dC - (dC /. gZero) - dTH], "Hessian G (part linear in q) = Homma Eq. (6)"];

(* ------------------------------------------------------------------ *)
(* 3. Prefactors in physical units against Homma's printed forms      *)
(* ------------------------------------------------------------------ *)
(* Symbols: ma, mb masses, qa, qb charges, n density, T temperature,
   ep vacuum permittivity, lnL Coulomb log, qv heat-flux vector. *)
physAssume = {ma > 0, mb > 0, qa > 0, qb > 0, n > 0, T > 0, ep > 0, lnL > 0};
vtb = Sqrt[2 T/mb];
Gam = qa^2 qb^2 lnL/(4 Pi ep^2 ma^2);
PA = Gam (1 + ma/mb) n/vtb^2;
PD = Gam n/(2 vtb);
CH = lnL (qa qb/ep)^2 (1 + mb/ma);          (* Homma's C *)
LH = (qa qb/(ma ep))^2 lnL;                 (* Homma's L_SI^(a/b) *)
qv = {qx, qy, qz};
gOfq = mb vtb qv/(n T^2);
(* Eq. (2): F0 = -C n/(4 pi 2 T) (Phi - x Phi')/x^3 x;
   F_gradT = -C/(10 pi^(3/2)) sqrt(mb/(2T)) exp(-x^2)/T^2 [q - 2(q.x) x] *)
F0H = -CH n/(4 Pi 2 T) (Erf[rr] - rr phiP[rr])/rr^3 uu;
FTH = -CH/(10 Pi^(3/2)) Sqrt[mb/(2 T)] Exp[-rr^2]/T^2 (qv - 2 qv.uu uu);
D0H = n/(16 Pi) Sqrt[mb/(2 T)] LH (dMH 2);
DTH = mb/(40 Pi T^2) LH ((-3/(2 rr^5) Erf[rr] + (1/rr^2 + 3/(2 rr^4)) phiP[rr])
        (Outer[Times, qv, uu] + Outer[Times, uu, qv] + qv.uu id3)
      + (15/(2 rr^7) Erf[rr] + (-2/rr^2 - 5/rr^4 - 15/(2 rr^6)) phiP[rr]) qv.uu uo);
physZeroQ[expr_] := AllTrue[Flatten[{Simplify[expr, Join[physAssume, {cartAssume}]]}], # === 0 &];
check[physZeroQ[ma PA aMH - F0H], "m_a P_A aM = Homma F0 (prefactor C n/(8 pi T))"];
check[physZeroQ[ma PA (aTH /. Thread[gg -> gOfq]) - FTH],
      "m_a P_A aT = Homma F_gradT (prefactor incl. sqrt(m_b/(2 T_b)))"];
check[physZeroQ[PD dMH - D0H], "P_D dM = Homma D0 (Eq. 5)"];
check[physZeroQ[PD (dTH /. Thread[gg -> gOfq]) - DTH], "P_D dT = Homma D_gradT (Eq. 6)"];
(* ------------------------------------------------------------------ *)
(* 4. Scalar kernel functions and the small-r series                  *)
(* ------------------------------------------------------------------ *)
check[zeroQ[-alphaF[rr] uu - aMH], "alpha reproduces aM"];
check[zeroQ[d1F[rr] id3 + d2F[rr] uo - dMH], "d1, d2 reproduce dM"];
check[zeroQ[t1F[rr] (gu + gg.uu id3) + t2F[rr] gg.uu uo - dTH], "t1, t2 reproduce dT"];
(* Truncation error at the switch (exact arithmetic, next three terms) and
   round-off of the closed form at the switch in machine precision. *)
Module[{tail, xs = xSwitch, rel, clos},
  Do[
    tail = Abs[Sum[SeriesCoefficient[kernelFuns[[i]][x], {x, 0, 2 k}] xs^(2 k),
                   {k, nSeries + 1, nSeries + 3}]];
    rel = N[tail/Abs[kernelFuns[[i]][xs]], 5];
    clos = Max[Table[Abs[(kernelFuns[[i]][N[xx]] - kernelFuns[[i]][SetPrecision[xx, 40]])
                         /kernelFuns[[i]][SetPrecision[xx, 40]]],
                     {xx, {xs, 1.01 xs, 1.1 xs, 1.5 xs}}]];
    Print[kernelNames[[i]], ": series truncation at switch ", rel,
          ", closed-form relative round-off just above switch ", N[clos, 3]];
    check[rel < 1*^-16, kernelNames[[i]] <> " series converged at switch"],
    {i, 5}]];

(* ------------------------------------------------------------------ *)
(* 5. Maxwellian average of the thermal-force kernel (B03, B07)       *)
(* ------------------------------------------------------------------ *)
(* <exp(-r^2)(g - 2(g.u)u)> over a Maxwellian in u of variance 1/(2 beta)
   per component is g (beta/(1+beta))^(5/2). *)
Module[{bet, sig2, lhs},
  sig2 = 1/(2 bet);
  lhs = Integrate[Exp[-(w1^2 + w2^2 + w3^2)] (1 - 2 w1^2) Exp[-(w1^2 + w2^2 + w3^2)/(2 sig2)]/(2 Pi sig2)^(3/2),
          {w1, -Infinity, Infinity}, {w2, -Infinity, Infinity}, {w3, -Infinity, Infinity},
          Assumptions -> bet > 0];
  check[Simplify[lhs - (bet/(1 + bet))^(5/2), bet > 0] === 0, "Maxwellian average of thermal-force kernel"]];

Print["all proofs passed"];
Exit[0];
