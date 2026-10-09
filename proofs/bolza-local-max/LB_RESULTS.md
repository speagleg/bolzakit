# LB_RESULTS: SU13-free certified lower bound for λ₁(Bolza)

2026-10-09. Plan: `LB_PLAN.md` (route (c), Lehmann–Goerisch with an H(div) flux). Labels as in PROOF.md.

## 0. Result

**CERTIFIED (computer-assisted; no SU13, no FBP multiplicity theorem):**

  **λ₁(B) ≥ 3.8388757010527243** (Arb: [3.83887570105272435423 ± 6.1e-38]; `lb_bolza.py 48` → `lb_bolza_48.json`).

  *Rounding (review_lb_1 fix F2).* The headline 3.8388757010527243 is the Arb lower endpoint
  3.83887570105272435423… rounded **downward** (truncated), so it is a valid lower bound. The json field
  `certified_lambda1_lower` = 3.8388757010527246 is `float(lb.lower())`, rounded to nearest, and lies ≈2.5e-16
  *above* the Arb lower endpoint; it must not be quoted as a bound. Quote either the headline or the Arb string
  `lb_enclosure` = [3.8388757010527243542306205221870856407 +/- 6.06e-38]. No downstream number is affected: the
  downstream λ-ball starts at 3.838875701026…, below both.

The same certificate gives λ₂, λ₃ ≥ 3.8388781485 (float evaluation of the k = 1, 2 bounds; only k = 3 is
Arb-certified, and it alone implies λ₁,₂,₃ ≥ 3.8388757).

- Required: λ₁ ≥ 3.8388029653 (D1 threshold, Lemma 8 at n = 96). Margin 7.3e-5.
- Desired: λ₁ ≥ 3.83881. Margin 6.6e-5.
- With the P2 upper bound: **λ₁ ∈ [3.8388757010, 3.8388872830]** (width 1.16e-5), consistent with SU13's
  [3.83888625884, 3.83888825884] (not used).
- Loss against the true value (≈ 3.83888726): 1.16e-5 at n = 48; 1.84e-4 at n = 24 (h⁴ convergence, as on the
  torus).

**Inputs used (complete list):**
1. Lemma 1(a): Spec(B) ∩ (0, 3.3] = ∅ and #Spec(B) ∩ [3.3, 4.0] ≤ 3 (trace-formula LP, `lp_count.py`; cites
   Buser/FBP trace formula normalisation, systole and 48 geodesics of FBP §7.4). Hence **λ₄ > 4.0**. Nothing
   else about the spectrum is used. (Not used: SU13, FBP Lemmas 7.2–7.3, any upper bound.)
2. Lehmann–Goerisch theorem (QUOTED): Liu–Plum arXiv:2512.23182 Thm 4.1/4.3, assumptions A1–A5 (§4.1), from
   Goerisch–Albrecht / Goerisch–Haunhorst [21, 22 there]. Statement quoted in LB_PLAN §1(c). A5 includes the
   requirement **ρ > 0**, met here with ρ = 5 (review_lb_1 row 1a; fix F1).
   *Second source (fix F1):* Behnke–Mertins–Plum–Wieners (2000), Thm 2.1, as reproduced in Vejchodský
   arXiv:1801.02561 §2. It is stated in exactly the shifted form used here: with shift γ = 1 their matrices are
   A0 − ρA1 and B = A0 − 2ρA1 + ρ²A2 (A2 = ‖σ‖² + ‖v + div σ‖²), and their conclusion for the unshifted problem is
   λ₁ ≥ 4 − 5/(1 − ν̄) (review_lb_1 row 1b, CONFIRMED). Both sources are secondary with respect to Goerisch's
   original papers, which were not read (by the author or the reviewer).
3. Min-max / spectral theorem for the Laplacian on a closed surface (compact resolvent; gives A2 of Liu–Plum).
4. The exact Bolza octagon and the mesh of `octagon.Mesh` / `cert_rr_p2.positions_jet` (unchanged, imported),
   and the hyperboloid-barycentric charts (PROOF §1).
5. python-flint 0.9.0 Arb ball arithmetic (prec 128).

## 1. The argument
- D = H¹(B), M(u,v) = q(u,v) + m(u,v) (positive definite), N = m. Eigenvalues λ̃ = λ + 1:
  λ̃₁ = 1, λ̃₂,₃,₄ = λ₁,₂,₃ + 1, λ̃₅ = λ₄ + 1 > 5 (input 1). Take **ρ = 5, m = 4, n = 3**.
- X = L²(TB) × L²(B), T u = (∇u, u), b_G the L² product: b_G(Tu,Tv) = M(u,v) (A3).
- Trial v₁, v₂, v₃: conforming P2 functions on the n = 48 mesh (float eigenvectors of the discrete problem;
  their coefficients are exact binary data).
- w_i = (σ_i, v_i + div σ_i), σ_i ∈ RT₂ (j = √g σ, raw space P₂² + x·P̃₂ in each reference simplex; DOFs = outward
  flux density per unit edge parameter s ∈ {0, ½, 1}, s measured from the endpoint with the smaller global label,
  plus 6 interior moments). Equal-and-opposite edge DOFs ⇒ σ ∈ H(div, B): the outward flux through an edge
  segment is intrinsic (ι_σ dA = j^a db − j^b da), and the edge point W = sP_A + (1−s)P_B is the same from both
  sides, also across the projective side pairings (isometries are linear on W and preserve s), exactly as for
  P2 conformity. Hence ∫σ·∇φ = −∫(div σ)φ for every φ ∈ H¹(B), and
  **b_G(w_i, Tφ) = ∫σ_i·∇φ + ∫(v_i + div σ_i)φ = N(v_i, φ) exactly (A4)**, for every σ_i. The σ_i (float
  minimisers of b_G(w_i,w_i), sparse SPD solve) only affect sharpness.
- A0 = [M(v_i,v_j)], A1 = [m(v_i,v_j)], A2 = [⟨σ_i,σ_j⟩ + ⟨v_i + div σ_i, v_j + div σ_j⟩], all summed over the
  18432 elements in Arb:
  - |σ|² dA = jᵀ(LLᵀ − QG)j/(D√Q) da db, (div σ)·√g = ∂_a j^a + ∂_b j^b, √g = D Q^{−3/2}.
    The metric identity g = (LLᵀ − QG)/Q² (hence g/√g = (LLᵀ − QG)/(D√Q)) and det(LLᵀ − QG) = Q D² were checked
    symbolically (`scratch_LB/gram_identity.py` → `.out`).
  - Element integrals ∫_simplex poly·Q^{−s} (s = ½, 3/2), as in cert_rr_p2.py, with K = 14 series terms plus the
    Lagrange remainder; exact monomial integration; exact rational RT₂/P2 basis tables.
- A = A0 − ρA1, B = A0 − 2ρA1 + ρ²A2. Certified in Arb: **B ≻ 0** (Sylvester) and **ν̄B − A ≻ 0** with
  ν̄ = −30.031942622361008 < 0, so all three ν_k < ν̄ < 0 (q = 3). Liu–Plum Thm 4.3 with k = 3:
  λ₁ + 1 = λ̃₂ ≥ ρ − ρ/(1 − ν̄), i.e. λ₁ ≥ 4 − 5/(1 − ν̄) = 3.83887570105272…

## 2. Runs

| run | n | elements | RT dofs | float P2 λ₁ | certified λ₁ ≥ | time |
|---|---|---|---|---|---|---|
| `lb_bolza.py 24` | 24 | 4608 | 48384 | 3.8388920838 | 3.8387028952 (fails D1 threshold) | 60 s |
| `lb_bolza.py 48` | 48 | 18432 | 193536 | 3.8388875620 | **3.8388757011** | 211 s |

Outputs: `lb_bolza_24.out/.json`, `lb_bolza_48.out/.json` (A0, A1, A2 enclosures are in the json).
ν (float, n = 48): −30.03241, −30.03241, −30.03194.

## 3. Validation and checks
- **Known answer (flat torus ℝ²/ℤ², same code with Q ≡ 1, M = JᵀJ):** λ₁ = 4π², multiplicity 4, λ₅ = 8π²;
  n = 4 trial functions, ρ = 79.9 ≤ 1 + 8π². `lb_torus.out`:

  | N | certified λ₁ ≥ | 4π² − bound |
  |---|---|---|
  | 4 | 39.2041077545 | 2.7e-1 |
  | 8 | 39.4582495177 | 2.0e-2 |
  | 16 | 39.4770992124 | 1.3e-3 |
  | 32 | 39.4783342462 | 8.3e-5 |

  All below 4π² = 39.4784176044, converging at O(h⁴).
- **Global Green identity** Σ_T ∫(σ·∇v + (div σ)v) = 0, evaluated exactly in Arb for the certified σ_i, v_j:
  |residual| ≤ 7.4e-36 (n = 48). It tests the edge-DOF combinatorics (orientation, s-flip, gluing).
- **Negative control (conformity is load-bearing):** forcing all edge signs to +1 (non-conforming σ) makes B
  indefinite on the torus (`scratch_LB/neg_sign.out`: "B not positive definite"); by Liu–Plum Thm 4.2 this
  cannot happen when A4 holds and ρ > Λ_n.
- **Independent element check:** P2 stiffness/mass vs `cert_rr_p2.element_p2` agree to 2e-24 (in each run);
  the RT₂ matrix vs brute-force Gauss quadrature with a finite-difference metric g_ij = −B(∂X,∂X) agrees to
  ≤ 5.7e-10 relative on 4 elements (`scratch_LB/xcheck_S.out`; float check, FD-limited).
- **ρ is an input, not a check:** with ρ = 6.5 > λ₄ + 1 ≈ 6.354 (hypothesis violated) the run still returns
  3.8388694 (`lb_bolza_24_rho6.5.out`), below the true value because the trial functions have negligible λ₄
  components. So the computation does not detect a wrong ρ; the λ₄ > 4.0 input from the LP is essential.
- Consistency: 3.8388757 ≤ P2 upper bound 3.8388872830.

## 4. Downstream certificates with the SU13-free ball
λ-ball [3.8388757010263, 3.8388872830264] = arb("3.83888149202636215", "5.791e-6"), which contains
[certified lower bound, P2 upper bound]. Copies `lb_cert_sign.py`, `lb_cert_first_order.py` differ from the
originals only in this line and the output file name (`diff cert_sign.py lb_cert_sign.py`: line 34 `LAM = arb(...)` and the json name; `diff cert_first_order.py lb_cert_first_order.py`: line 290 `lam1 = arb(...)` and the json name). Originals untouched.

| certificate | original (SU13 ball) | SU13-free ball | status |
|---|---|---|---|
| Lemma 8, g″(0) upper (n = 96, P2, SYM=1) | −12.850476 | **−9.341154** | still < 0: CERTIFIED |
| Lemma 5, b(e,e)/m(u) upper, n = 96 (P1) | −1.02824 | **−0.99892** | still < 0: CERTIFIED |
| Lemma 5, b(e,e)/m(u) upper, n = 192 (P1) | −3.087569 | **−3.038347** | still < 0: CERTIFIED (A ≠ 0) |

*SYM note (corrected after review_lb_1, fix F3).* The Lemma 8 row is SYM=1: `lb_cert_sign.py` defaults to SYM=1
and `lb_cert_sign_96.json` records `"sym": true`. The Lemma 5 rows were previously labelled "SYM=1", but
`lb_cert_first_order.py` defaults to **SYM=0** and its output (`lb_cert_first_order.out`,
`lb_cert_first_order_96_192.json`) does not record which setting was used; the label is therefore withdrawn. The
n = 192 value does not depend on it: the reviewer's re-run with SYM=0 (full computation, no symmetry shortcut,
3119 s) gives b(e,e)/m ≤ −3.0383476316928864 (`scratch_RLB/lb_cert_first_order_192.json`), against
−3.0383476242423058 in the author's run; the ~7e-9 difference comes from the untrusted float eigenvectors. The
reviewer's SYM=1 re-run of the sign certificate gives −9.341154566965994 (author: −9.341154560446736).

Outputs: `lb_cert_sign.out`, `lb_cert_sign_96.json`, `lb_cert_first_order.out`, `lb_cert_first_order_96_192.json`.

## 5. Consequence for PROOF.md (integrated 2026-10-09 as PROOF Lemma 1(d), after review_lb_1 ACCEPT)
D1 (SU13) can be replaced by: "λ₁ ≥ 3.8388757010, CERTIFIED by Lehmann–Goerisch (LB_RESULTS.md), using only
Lemma 1(a) and the quoted Lehmann–Goerisch theorem". Lemma 1(c)'s identification of SU13's eigenvalue is no
longer needed; FBP Lemmas 7.2–7.3 remain needed only for multiplicity 3 (λ₁ = λ₂ = λ₃), not for the lower
bound. Lemma 8's certified constant becomes −9.3411 (from −12.8504); Lemma 5's becomes −3.0383·m(u) (from −3.0875). (Integration note 2026-10-09: these constants, and the §4 table entries, were previously rounded to nearest, e.g. −9.3412 and −3.03835, which slightly overstates an upper bound; they are now rounded toward zero from the json values −9.341154560446736 and −3.0383476242423058.) The new external
dependency is the Lehmann–Goerisch theorem (classical, quoted from Liu–Plum 2025 §4.1; Goerisch's papers were
not read by us).

Not done / limitations: the Lehmann–Goerisch statement is QUOTED from secondary sources (Liu–Plum 2025 §4.1; BMPW
2000 Thm 2.1 via Vejchodský 2018); Goerisch's original papers were not read. n = 96 (sharper, ≈ 7e-7 loss
expected) not run.

**Second read:** `review_lb_1.md` — **ACCEPT**, with three cosmetic fixes F1–F3 (none changes a number), applied
above (F1: §0 input 2 and LB_PLAN §1(c); F2: §0 rounding note; F3: §4 SYM note).
