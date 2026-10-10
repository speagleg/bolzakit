# LEAN.md — Lean 4 / Mathlib skeleton of the item-23 proof

**In this repository.** This is the record of the original item-23 formalisation. The files are in
[`lean/BolzaKit/LocalMax/`](../../lean/BolzaKit/LocalMax/) (modules `BolzaKit.LocalMax.*`, root
`lean/BolzaKit/LocalMax.lean`); item 23 called the directory `BolzaLean/` and the root `BolzaLean.lean`, and the
record below keeps those names. Two differences from the record:
- the field `BolzaInputs.g_sign` now assumes the SU13-free constant g″(0) ≤ −9.3411 (Lemma 8 with the
  Lehmann–Goerisch λ₁ ball) instead of −12.8504; `Cert.sign_su13_implies_lg` checks that the SU13 value also
  discharges it. The Lean files do not check the Lehmann–Goerisch certificate's own arithmetic
  (−9.3411481 ≤ −9.3411); that link is made by hand, like the other certificate links;
- the build runs through `lake build` (CI) or `lean/build_local.sh`, which replaces item 23's `build.sh` (§4).

Project: `items/23-bolza-lambda1-extremality/lean/` (Lean 4.34.1, Mathlib v4.34.1, reused read-only from
`bolzakit/lean/.lake/packages`). About 2250 lines in 8 files. The build has **no `sorry`**, **no `axiom`
declarations** and **no `native_decide`**. Every analytic or geometric input is a named field (an explicit
hypothesis) of the structure `Bolza.BolzaInputs` in `BolzaLean/Assembly.lean` (here `lean/BolzaKit/LocalMax/Assembly.lean`). Each field's docstring cites
the PROOF.md lemma or dependency it stands for.

## 0. One-paragraph summary

`Bolza.bolza_strict_local_max (I : BolzaInputs) : ∀ᶠ X in 𝓝 I.B, X ≠ I.B → I.lam1 X < I.lam1 I.B` is
proved in Lean from the inputs. The inputs are F1, Lemma N equivariance, the chart of Lemma 9, A ≠ 0,
Lemma 11 invariance, E irreducible, the group relations, the w_B/F_B data, and the Ky Fan comparison plus
the certified g″(0) ≤ −12.8504. Everything finite-dimensional between those inputs and the conclusion is
machine-checked:
- Lemma 10, in C² form, including a uniform second-order Taylor theorem;
- Lemma 4 (tr A ≡ 0);
- Lemma 6(iv)–(vi): K₀ ⊂ ker A, ker A ⊂ K₀ and D²Λ|_{K₀} = c·I;
- Lemma 7's exact line;
- the step c < 0;
- the Ky Fan step (Λ∘γ)″ ≤ g″;
- the second-order chain rule that turns it into D²Λ(w_B, w_B) ≤ −12.8504.

The exact representation theory (the 48-element group, characters and multiplicities) is checked by kernel
evaluation. The certificate *final arithmetic* is checked by `norm_num`; the Arb computations themselves
are not. A toy instance of `BolzaInputs` (`Bolza.Model.toyInputs`) shows that the hypotheses are jointly
satisfiable, so the main theorem is not vacuous. The instance says nothing about the Bolza surface.

## 1. Findings while formalizing (read these first)

1. **The sign bound cannot be "assembled from its components" (`Cert.sign_rows_not_additive`).** The
   central upper end is −14.2742 + 2.23e-5. Adding the six printed one-at-a-time error rows (sum ≈ 1.4100)
   gives about **−12.864**. That is *below* the certified joint upper end **−12.8505**, whose error is
   ≈ 1.4237. So the row sum is not an upper bound, and only the joint Arb ball certifies the sign. This is
   consistent with PROOF.md Lemma 8 step 5, which already labels the rows "one-at-a-time". It does mean
   that no re-derivation of −12.85 from the printed table is valid. Lean checks only the joint value
   (`Cert.sign_joint`: −12.850476214196528 ≤ −12.8504 < 0).
2. **Lemma 10 needs less than stated.**
   - The hypothesis "A ≠ 0" is not needed for the criterion (`local_max_criterion` does not assume it).
   - Symmetry of M is needed only to restate the conclusion as an eigenvalue inequality
     (`exists_eigenvalue_lt`).
   - The constant c_A > 0 comes from an explicit elementary bound (`traceless_symm_exists_neg`): a
     traceless symmetric 3×3 matrix S has a {0,±1} test vector with xᵀSx ≤ −‖S‖_sup/4. Combined with the
     antilipschitz bound of A on a complement of ker A, this gives c_A.
   - In case (ii) the proof uses tr A(v) = 0 directly. PROOF.md's tr A(w) is the same quantity, since
     A(k) = 0.
   - No error was found in Lemma 10.
3. **Lemma 4 needs only r.** In TT block coordinates, r acts by the same block on T₋ and T₊, and that block
   has no fixed vector. So tr A ≡ 0 follows from Lemma 3(b) equivariance for the single octagon rotation r
   (`trace_A_zero`, `no_invariant_functional`). That equivariance is already covered by the relabelling line
   and does not need OPEN-1.
4. **What each step needs from the group, as the formalization shows it.**
   - W₀ is irreducible using only t and κ (`Wplus_irreducible`).
   - The invariant forms on K₀ are scalar using r, t and κ (`invariant_form_Kminus`).
   - So both Lemma 6(v) and Lemma 6(vi) genuinely use the non-octagon element t, as PROOF.md's OPEN-1
     discussion says. For A, that equivariance comes from Lemma 3(b), which is intrinsic on TT. For the
     Hessian it comes from Lemma 11, the OPEN-1 closure.
5. **Lemma 7 depends on s_mid being the right reflection type (`Rep.smid_line_in_Kminus`).**
   - In the 48-element image group there are two s_mid-type reflections (fixing 4 Weierstrass points) and
     two s_vert-type reflections. For **both** s_mid-type ones, the line {r v = −v, s v = v} lies in
     K₀ = T₋.
   - For the s_vert-type ones, the line lies in **W₀** instead (Python check while formalizing).
   - So the identification of F_B's symmetry as an s_mid-type reflection is load-bearing. This is PROOF.md's
     "corrected wording (review_theory_1 row 1d)".
   - In cert_chars.py the chosen `s_mid` is literally κ: its 6×6 matrix equals `real6(kp)`.
6. **SU13 slack.** From the guaranteed SU13 lower end 3.83888625884, the slack over the sign-certificate
   threshold 3.8388029653 is 8.329e-5, i.e. **83.3×** SU13's 1e-6 (`Cert.su13_slack`). PROOF.md says
   "≈ 8.4e-5 / 84×" in Lemma 8(6), which is measured from the upper end 3.8388872830. Its "≈ 80×" in
   §Dependencies is right. This is cosmetic.
7. **Group relations used by Lemma 6(iv).** These are a² = J, gag⁻¹ ∈ {a, aJ} for all g, and J central. They
   are hypotheses in Lean, because there is no Lean model of the 96-element group G. While formalizing I
   re-checked them with exact arithmetic using cert_chars.py's own `gmul`, in a scratch script outside the
   repo: all three hold, for both lifts of x ↦ −1/x̄.
8. **The step from characters to irreducibility is not formalized.** The statement "⟨χ,χ⟩ = 1 ⇒ absolutely
   irreducible" is not formalized (Mathlib has no such converse). The facts the proof actually *uses* are
   proved directly over ℝ instead:
   - Schur-type invariant forms;
   - W₀ has no invariant subspace;
   - the involution Schur lemma.

   So the character table is a consistency check, not a logical input.

No statement of PROOF.md was found to be false.

## 2. Status table

Status labels:
- **AXIOM** = an explicit hypothesis, i.e. a field of `BolzaInputs`. There are no Lean `axiom`s.
- **ARITHMETIC-CHECKED** = final rational arithmetic, or exact integer data checked by kernel evaluation.
  The Arb/Python computations upstream of it are not verified.

File paths are relative to `lean/BolzaLean/` (here `lean/BolzaKit/LocalMax/`).

| PROOF.md item | Lean status | file : declarations |
|---|---|---|
| Lemma 1(a): trace formula and LP soundness (D2, D3) | AXIOM (implicit: enters only through the isolation/gap used by F1, Lemma 5 and Lemma 8, i.e. via `M_C2`, `lam1_le_rayleigh`, `A_ne_zero`, `g_sign`) | — |
| Lemma 1(a): final LP numbers 0.912612 < 1, 3.869661 < 4 ⇒ count ≤ 3 | ARITHMETIC-CHECKED | Certificates : `lpA_upper_lt_one`, `lpB_upper_lt_four`, `count_le_three` |
| Lemma 1(b): FEM conformity and min-max (D9) | AXIOM (implicit, as above) | — |
| Lemma 1(b): μ₃ ≤ 3.8391848048, gap ≥ 0.1608 > r₀ | ARITHMETIC-CHECKED | Certificates : `p1_mu3`, `gap_ge`, `gap_ge_r0` |
| Lemma 1(c): SU13 (D1), FBP 7.2–7.3 (D4), multiplicity 3, E irreducible | AXIOM | Assembly : `BolzaInputs.lamStar`, `lam1_B`, `E_irred` |
| Lemma 1(c): interval [3.83888625884, 3.8388872830]; script ball covers it (G1) | ARITHMETIC-CHECKED | Certificates : `p2_mu3`, `lambda1_interval`, `script_ball_covers`, `kmp_consistent` |
| Lemma 2 (F1): analytic symmetric cluster matrix, M(0) = λ*I, λ₁ = λ_min M | AXIOM | Assembly : `M`, `M_C2`, `M_sym`, `M_zero`, `lam1_le_rayleigh` |
| Lemma 3 (Lemma N): A intrinsic and G-equivariant | AXIOM | Assembly : `ρr…ρJ`, `ρr_orth`, `equiv_r`, `equiv_t`, `equiv_k`, `equiv_a` |
| Lemma 4: tr A ≡ 0 | PROVED-IN-LEAN (from `equiv_r` only) | Assembly : `trace_A_zero`; RepReal : `no_invariant_functional`; RepTheory : `chiT_one_ip` |
| Lemma 4 proof 1: d(angle sum)/dR ≠ 0 | ARITHMETIC-CHECKED | Certificates : `dangle_ne_zero` |
| Lemma 5: A ≠ 0 | AXIOM + ARITHMETIC-CHECKED (b(e,e) ≤ −3.0876 < 0, joint and re-assembled) | Assembly : `A_ne_zero`; Certificates : `lemma5_joint_192`, `lemma5_assembled_192`, `lemma5_assembled_96`, `lemma5_threshold` |
| Lemma 6(i): curve model y² = x⁵ − x, \|G\| = 96, T_B ≅ Re Q(B) (D5, D6) | AXIOM (implicit: the identification of the matrices with the geometric action is QUOTED) | — |
| Lemma 6(i): exact generator matrices on T_B, change of basis to blocks, a is a word | ARITHMETIC-CHECKED (kernel `decide`, integer data from cert_chars.py) | RepTheory : `R6`, `T6x2`, `K6`, `A6`, `P`, `P_invertible`, `conj_r`, `conj_t`, `conj_k`, `conj_a`, `a_is_word` |
| Lemma 6(i): the image group has order 48 and is exactly the generated group | ARITHMETIC-CHECKED + PROVED | RepTheory : `Hlist_length`, `H_closed`, `H_words`, `H_proj_nodup`, `InH_iff` |
| Lemma 6(ii): χ_T values; ⟨χ_T,1⟩ = 0; ⟨χ_T,χ_T⟩ = 2; χ_T = 0 on reversing elements | ARITHMETIC-CHECKED | RepTheory : `chiT_values`, `chiT_one_ip`, `chiT_norm`, `chiT_reversing` |
| Lemma 6(iii): T₊/T₋ splitting, norms 1, 1, 0 | ARITHMETIC-CHECKED | RepTheory : `conj_a`, `chiM_norm`, `chiP_norm`, `chiM_chiP` |
| Lemma 6(iv): K₀ ⊂ ker A (antipodal argument) | PROVED-IN-LEAN, given the group-relation hypotheses | RepReal : `involution_schur`, `rhoE_a_scalar`; Assembly : `A_vanishes_K0` (inputs `J_sq`, `J_central`, `a_sq`, `a_central_modJ`) |
| Lemma 6(v): ker A = K₀ | PROVED-IN-LEAN | RepReal : `Wplus_irreducible`; Assembly : `kerA_sub_K0` |
| Lemma 6(vi): invariant forms on K₀ are c·g_inv | PROVED-IN-LEAN | RepReal : `invariant_form_Kminus` |
| Lemma 6 consistency: Sym₀²E ∌ K₀ and ⊇ W₀ once, for E ∈ {X₁, X₂, O, Oε}; candidates orthonormal | ARITHMETIC-CHECKED | RepTheory : `sym0_*` (8 theorems), `chiO_norm`, `chiOe_norm`, `candidates_orthogonal` |
| Lemma 7: exact line {r = −1, s_mid = +1} ⊂ K₀ | PROVED-IN-LEAN + ARITHMETIC-CHECKED (every s_mid-type reflection) | RepReal : `line_r_smid`; RepTheory : `smid_line_in_Kminus`, `gK_is_smid` |
| Lemma 7: w_B ≠ 0, relabelling line r·w_B = −w_B, s_mid·w_B = w_B | AXIOM | Assembly : `wB`, `wB_ne`, `wB_r`, `wB_smid` |
| Lemma 7: surface condition C²xy − C(x+y) − 1 = 0 at p = 0; pairing trace 11 + 8√2 in the certified ball | PROVED-IN-LEAN | Certificates : `FB_surface_condition`, `pairing_trace` |
| Lemma 8 step 1: operator Ky Fan Λ∘γ ≤ g, equality at 0 (D9) | AXIOM | Assembly : `g`, `g_diff`, `g_diff2`, `kyfan_le`, `kyfan_eq` |
| Lemma 8 step 1: Ky Fan ⇒ (Λ∘γ)″(0) ≤ g″(0) | PROVED-IN-LEAN | KyFan : `second_deriv_nonneg_of_min`, `second_deriv_le_of_le` |
| Lemma 8: finite-dimensional Ky Fan (QUOTED in D9) | PROVED-IN-LEAN | KyFan : `bathtub`, `bathtub_exists`, `ky_fan` |
| Lemma 8: certified g″(0) ≤ −12.850476 (Arb, P2, n = 96) | AXIOM + ARITHMETIC-CHECKED (rounding to −12.8504; non-additivity finding) | Assembly : `g_sign`; Certificates : `sign_joint`, `sign_rows_not_additive`, `su13_slack` |
| Lemma 8 "Consequence": (Λ∘γ)″ = D²Λ(w_B, w_B) because dΛ = 0 | PROVED-IN-LEAN | Assembly : `second_deriv_trace_comp`, `wB_sign` |
| Lemma 9(a)(b): rank J = 5, rank dπ = 6 (determinant enclosures exclude 0) | ARITHMETIC-CHECKED | Certificates : `rankJ_minor`, `dpi_minor`, `dpi_cross` |
| Lemma 9(c)(d): chart, homeomorphism onto a neighbourhood, σ, retraction p, γ = p∘x | AXIOM | Assembly : `T2`, `topT2`, `B`, `lam1`, `chart`, `chart_zero`, `chart_nhds`, `γ`, `γ_zero`, `γ_diff`, `γ_diff2`, `γ_tangent` |
| Lemma 10: local-max criterion | **PROVED-IN-LEAN** | LocalMax : `local_max_criterion`, `traceless_symm_exists_neg`, `case1_arith`, `case2_arith`, `exists_eigenvalue_lt`; Taylor : `taylor2_littleO`, `fderiv_isSym3`, `local_max_criterion_C2` |
| Lemma 11 (OPEN-1 closure): G-invariance of D²Λ in the chart | AXIOM | Assembly : `hess_inv_r`, `hess_inv_t`, `hess_inv_k` |
| Theorem: assembly | **PROVED-IN-LEAN** (conditional on `BolzaInputs`) | Assembly : `hess_K0_neg`, `bolza_strict_local_max` |
| Non-vacuity of the inputs | PROVED-IN-LEAN (toy model, not geometric) | Model : `toyInputs` |
| D7, D8, D10, D11, D13, D14 (Fischer–Marsden, elliptic regularity, Kato, Poincaré, software, analyticity of π) | NOT ATTEMPTED (they sit inside the AXIOM fields above) | — |

## 3. The inputs (fields of `Bolza.BolzaInputs`)

There are no Lean `axiom` declarations. The main theorem takes `I : BolzaInputs`, and the inputs are these
fields; the full docstrings are in `Assembly.lean`.

| Field(s) | Stands for |
|---|---|
| `T2`, `topT2`, `B`, `lam1` | T₂ and λ₁ (D12) |
| `chart`, `chart_zero`, `chart_nhds` | Lemma 9(c): π\|_S is a homeomorphism onto a neighbourhood of B |
| `M`, `lamStar`, `lam1_B`, `M_C2`, `M_sym`, `M_zero`, `lam1_le_rayleigh` | Lemma 2 (F1) on S in TT block coordinates, plus Lemma 1 |
| `ρr`, `ρt`, `ρk`, `ρa`, `ρJ`, `ρr_orth` | the representation of G on E (an m₀-orthonormal basis) |
| `equiv_r`, `equiv_t`, `equiv_k`, `equiv_a` | Lemma 3(b), Lemma N |
| `E_irred` | Lemma 1(c) via FBP 7.2–7.3 (D4) |
| `J_sq`, `J_central`, `a_sq`, `a_central_modJ` | group relations in G (D5 / Lemma 6(iii); exact in cert_chars) |
| `A_ne_zero` | Lemma 5 |
| `hess_inv_r`, `hess_inv_t`, `hess_inv_k` | Lemma 11 |
| `wB`, `wB_ne`, `wB_r`, `wB_smid` | Lemma 7 and the Lemma 3 relabelling line |
| `γ`, `γ_zero`, `γ_diff`, `γ_diff2`, `γ_tangent` | Lemma 9(d)1, 9(d)3: F_B in the chart |
| `g`, `g_diff`, `g_diff2`, `kyfan_le`, `kyfan_eq` | Lemma 8 step 1, Ky Fan (D9) |
| `g_sign` | Lemma 8, the certified number |

The coordinate convention is an input too. TT ≅ `TB = ℝ³ × ℝ³` (K₀ × W₀) is the cert_chars.py basis
conjugated by `P`, and that conjugation is kernel-checked (`Rep.conj_*`). The statement that this basis
*is* the geometric TT action is QUOTED (D5, D6).

What the Lean theorem does **not** cover:
- whether the fields hold for the actual Bolza surface. That is the content of F1, Lemma N, Lemma 9,
  Lemma 11, the Arb certificates, SU13 and FBP;
- the Arb computations themselves;
- the operator Ky Fan / min-max inequality in H¹.

## 4. Build

In this repository, build with `lake build` (as CI does) or with `lean/build_local.sh`, which compiles the same
files against prebuilt Mathlib oleans. The record below is the original item-23 build.

Verified route (read-only with respect to bolzakit). The oleans were copied to a local scratch directory for
speed, because memory-mapping them from a network or WSL-mounted filesystem is very slow:

```
cd items/23-bolza-lambda1-extremality/lean
PKGS=<dir containing <pkg>/lean copies of bolzakit/lean/.lake/packages/<pkg>/.lake/build/lib/lean> ./build.sh
# or simply ./build.sh  (reads bolzakit's oleans in place; same commands, much slower first load;
#  this variant was not re-run end-to-end)
```

`build.sh` runs `nice lean -o build/BolzaLean/<F>.olean` for
F = RepTheory, LocalMax, Taylor, RepReal, Certificates, KyFan, Assembly, Model, then `lean BolzaLean.lean`.
A clean build took 2 min 33 s, with 0 errors and 0 `sorry` (the full log, `lean/build.log`, stayed with item 23). The
`lakefile.lean` is provided for convenience and was not used; see the caution in it.

`#print axioms` output (from `BolzaLean.lean`):

```
'Bolza.bolza_strict_local_max' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.local_max_criterion' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.local_max_criterion_C2' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.taylor2_littleO' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.trace_A_zero' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.A_vanishes_K0' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.kerA_sub_K0' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.hess_K0_neg' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.wB_sign' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.Rep.InH_iff' depends on axioms: [propext, Quot.sound]
'Bolza.Rep.chiT_norm' depends on axioms: [propext]
'Bolza.Rep.sym0_O_K0' depends on axioms: [propext]
'Bolza.Rep.smid_line_in_Kminus' depends on axioms: [propext]
'Bolza.Rep.conj_t' depends on axioms: [propext]
'Bolza.invariant_form_Kminus' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.Wplus_irreducible' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.second_deriv_le_of_le' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.ky_fan' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.Cert.sign_joint' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.Cert.sign_rows_not_additive' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.Cert.pairing_trace' depends on axioms: [propext, Classical.choice, Quot.sound]
'Bolza.Model.toyInputs' depends on axioms: [propext, Classical.choice, Quot.sound]
```

Only Lean's three standard axioms appear. The kernel `decide` proofs use `decide +kernel`, which adds no
axiom (unlike `native_decide`).
