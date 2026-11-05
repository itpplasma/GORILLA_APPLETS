(* ::Package:: *)
(* Definitions shared by homma_hatayama.wl (generator) and
   homma_hatayama_proofs.wl (proofs). Normalisation: see homma_hatayama.wl. *)

$HistoryLength = 0;
fail[msg_] := (Print["FAILED: ", msg]; Exit[1]);
check[cond_, msg_] := If[TrueQ[cond], Print["ok: ", msg], fail[msg]];

phiP[x_] := 2/Sqrt[Pi] Exp[-x^2];   (* Homma's Phi'(x) *)
id3 = IdentityMatrix[3];

(* Scalar kernel functions, proved equal to Homma Eqs. (2), (5), (6) in
   homma_hatayama_proofs.wl:
   a  = -alpha(r) u - kT exp(-r^2) (g - 2 (g.u) u),     kT = 2/(5 sqrt(pi))
   d  = d1(r) I + d2(r) u u + t1(r) (g u + u g + (g.u) I) + t2(r) (g.u) u u *)
alphaF[x_] := (Erf[x] - x phiP[x])/x^3;
d1F[x_] := 1/2 ((2/x - 1/x^3) Erf[x] + phiP[x]/x^2);
d2F[x_] := 1/2 ((-2/x^3 + 3/x^5) Erf[x] - 3/x^4 phiP[x]);
t1F[x_] := 1/5 (-3/(2 x^5) Erf[x] + (1/x^2 + 3/(2 x^4)) phiP[x]);
t2F[x_] := 1/5 (15/(2 x^7) Erf[x] + (-2/x^2 - 5/x^4 - 15/(2 x^6)) phiP[x]);
kernelNames = {"alpha", "d1", "d2", "t1", "t2"};
kernelFuns = {alphaF, d1F, d2F, t1F, t2F};

(* Each function is even and analytic at 0: series in w = r^2. *)
nSeries = 12;          (* highest power of w kept *)
xSwitch = 1/2;         (* series for r < xSwitch, closed form above *)
