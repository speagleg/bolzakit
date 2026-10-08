/-
Item 05 — Golden geodesic and the systole of Γ(p).

Formalised here (complete proofs, no placeholders):
* `cayley_hamilton_two`, `trace_mul_self`, `trace_pow_add_two`, `trace_dvd_trace_pow_odd`:
  2×2 integer matrix facts (Cayley–Hamilton, tr(M²) = (tr M)² − 2 det M, the trace recursion,
  tr M ∣ tr M^(2i+1)).
* `trace_Q_pow`: tr Q^n = L_n for Q = [[1,1],[1,0]] and the Lucas numbers L_n (defined by recursion).
* `lucas_two_mul`: L_{2k} = L_k² − 2(−1)^k.
* `lucas_dvd_lucas_mul_odd`: L_r ∣ L_{r(2i+1)}  (Lemma 3 of the note).
* `lucas_strictMono`, `lucas_inj`: L is strictly increasing on n ≥ 1.
* `pow_two_of_lucas_prime`: if k is even and L_k is prime then k = 2^j, j ≥ 1 (Lemma 4).
* `trace_of_pmGamma`: A ∈ SL₂(ℤ), A ≡ εI (mod N), ε = ±1  ⇒  tr A = ε(2 − N²d) for some d ∈ ℤ.
* `systole_trace_bound` (C1, lower bound): such A with |tr A| > 2 has |tr A| ≥ N² − 2 (any N ≠ 0).
* `A0_det`, `A0_mem`, `A0_trace`: attainment by A₀ = [[1 − N², N], [−N, 1]] ∈ Γ(N).
* `exists_goldenExp`: for p prime some even n ≥ 2 has Q^n ∈ ±Γ(p) (finite-group argument),
  so `goldenN p` (= n(p), the least such n) is well defined.
* `golden_systolic_iff` (C2, trace form): for an odd prime p,
  L_{n(p)} = p² − 2  ↔  ∃ j ≥ 1, p = L_{2^j}.
* `golden_minimal_iff` (C2, extremal form): the golden lift has the least |trace| among all
  hyperbolic elements of ±Γ(p)  ↔  ∃ j ≥ 1, p = L_{2^j}.
* `goldenN_eq` : in that case n(p) = 2^{j+1}.
* `stab_golden` (axis stabiliser, algebraic part): A ∈ SL₂(ℤ) commuting with Q ⇒ A = ±Q^{2m} or
  ±(Q − 1)^{2m} = ±Q^{−2m}; the norm-±1 units of ℤ[φ] are handled by an integer descent (`lin_isPmQPow`).

NOT formalised (stated and proved on paper in work/note.tex):
* the geometric translation |tr| ↦ length 2·arccosh(|tr|/2), closed geodesics of Γ(p)\H ↔ primitive
  hyperbolic conjugacy classes, and the identification of ±Γ(p) and Γ(p) in PSL₂(ℤ);
* "A fixes φ ⇔ A commutes with Q" for A ∈ SL₂(ℤ), the no-endpoint-swap argument (torsion-freeness),
  and Stab_{Γ̄(p)}(axis) = ⟨Q^{n(p)}⟩ (needs Lemma 0);
* the Fibonacci description Q^n = [[F_{n+1},F_n],[F_n,F_{n−1}]], the rank of apparition α(p),
  n(p) = lcm(2, α(p)) (Lemma 0) and α(p) = 2^{j+1} in the systolic case.
-/
import Mathlib

namespace BolzaKit.GoldenGeodesic

open Matrix

/-- 2×2 integer matrices. -/
abbrev M2 := Matrix (Fin 2) (Fin 2) ℤ

/-- The golden matrix `Q = [[1,1],[1,0]]` (note `det Q = -1`, so `Q ∉ SL₂(ℤ)`). -/
def Q : M2 := !![1, 1; 1, 0]

/-- Lucas numbers: `L 0 = 2`, `L 1 = 1`, `L (n+2) = L (n+1) + L n`. -/
def lucas : ℕ → ℤ
  | 0 => 2
  | 1 => 1
  | (n + 2) => lucas (n + 1) + lucas n

/-- `A ∈ ±Γ(N)`: `A ≡ ε·I (mod N)` entrywise for a sign `ε = ±1`.
(Membership in `SL₂(ℤ)` is a separate hypothesis `A.det = 1`.) -/
def PmGamma (N : ℤ) (A : M2) : Prop :=
  ∃ ε : ℤ, (ε = 1 ∨ ε = -1) ∧ ∀ i j, N ∣ A i j - ε * (1 : M2) i j

/-! ### General 2×2 facts -/

theorem cayley_hamilton_two (M : M2) : M * M = M.trace • M - M.det • (1 : M2) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp only [Matrix.mul_apply, Fin.sum_univ_two, Matrix.trace_fin_two, Matrix.det_fin_two,
      Matrix.sub_apply, Matrix.smul_apply, smul_eq_mul, Matrix.one_apply_eq,
      Matrix.one_apply_ne (show (0 : Fin 2) ≠ 1 by decide),
      Matrix.one_apply_ne (show (1 : Fin 2) ≠ 0 by decide), Fin.zero_eta, Fin.mk_one,
      Fin.isValue] <;> ring

theorem trace_mul_self (M : M2) : (M * M).trace = M.trace ^ 2 - 2 * M.det := by
  rw [cayley_hamilton_two, Matrix.trace_sub, Matrix.trace_smul, Matrix.trace_smul,
    Matrix.trace_one]
  simp
  ring

theorem trace_pow_add_two (M : M2) (n : ℕ) :
    (M ^ (n + 2)).trace = M.trace * (M ^ (n + 1)).trace - M.det * (M ^ n).trace := by
  have h : M ^ (n + 2) = M.trace • M ^ (n + 1) - M.det • M ^ n := by
    rw [pow_add, pow_two, cayley_hamilton_two, mul_sub, mul_smul_comm, mul_smul_comm, mul_one,
      pow_succ]
  rw [h, Matrix.trace_sub, Matrix.trace_smul, Matrix.trace_smul, smul_eq_mul, smul_eq_mul]

theorem trace_dvd_trace_pow_odd (M : M2) (i : ℕ) : M.trace ∣ (M ^ (2 * i + 1)).trace := by
  induction i with
  | zero => simp
  | succ i ih =>
    have e : 2 * (i + 1) + 1 = (2 * i + 1) + 2 := by ring
    rw [e, trace_pow_add_two]
    exact dvd_sub (dvd_mul_right _ _) (dvd_mul_of_dvd_right ih _)

/-! ### Lucas numbers as traces of powers of `Q` -/

theorem Q_det : Q.det = -1 := by simp [Q, Matrix.det_fin_two]

theorem Q_trace : Q.trace = 1 := by simp [Q, Matrix.trace_fin_two]

theorem Q_pow_det (n : ℕ) : (Q ^ n).det = (-1) ^ n := by rw [Matrix.det_pow, Q_det]

theorem trace_Q_pow : ∀ n : ℕ, (Q ^ n).trace = lucas n
  | 0 => by simp [lucas]
  | 1 => by simp [lucas, Q_trace]
  | (n + 2) => by
    rw [trace_pow_add_two, Q_trace, Q_det, trace_Q_pow (n + 1), trace_Q_pow n, lucas]
    ring

theorem lucas_two_mul (k : ℕ) : lucas (2 * k) = lucas k ^ 2 - 2 * (-1) ^ k := by
  rw [← trace_Q_pow, ← trace_Q_pow, two_mul, pow_add, trace_mul_self, Q_pow_det]

/-- Lemma 3: `L_r ∣ L_{r m}` for odd `m = 2i+1`. -/
theorem lucas_dvd_lucas_mul_odd (r i : ℕ) : lucas r ∣ lucas (r * (2 * i + 1)) := by
  rw [← trace_Q_pow, ← trace_Q_pow, pow_mul]
  exact trace_dvd_trace_pow_odd _ _

theorem lucas_pos : ∀ n : ℕ, 1 ≤ lucas n
  | 0 => by simp [lucas]
  | 1 => by simp [lucas]
  | (n + 2) => by
    have h1 := lucas_pos (n + 1)
    have h2 := lucas_pos n
    rw [lucas]
    linarith

theorem lucas_lt_succ (n : ℕ) (hn : 1 ≤ n) : lucas n < lucas (n + 1) := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  have h := lucas_pos m
  rw [show m + 1 + 1 = m + 2 from rfl, lucas]
  linarith

theorem lucas_strictMono {m n : ℕ} (hm : 1 ≤ m) (h : m < n) : lucas m < lucas n := by
  induction n with
  | zero => omega
  | succ n ih =>
    by_cases h' : m < n
    · exact (ih h').trans (lucas_lt_succ n (by omega))
    · have : m = n := by omega
      subst this
      exact lucas_lt_succ m hm

theorem lucas_le_iff {m n : ℕ} (hm : 1 ≤ m) (hn : 1 ≤ n) : lucas m ≤ lucas n ↔ m ≤ n := by
  constructor
  · intro h
    by_contra h'
    exact absurd h (not_le.mpr (lucas_strictMono hn (by omega)))
  · intro h
    rcases Nat.lt_or_eq_of_le h with h | rfl
    · exact (lucas_strictMono hm h).le
    · exact le_rfl

theorem lucas_inj {m n : ℕ} (hm : 1 ≤ m) (hn : 1 ≤ n) (h : lucas m = lucas n) : m = n :=
  le_antisymm ((lucas_le_iff hm hn).mp h.le) ((lucas_le_iff hn hm).mp h.ge)

theorem lucas_two : lucas 2 = 3 := by simp [lucas]

theorem three_le_lucas {n : ℕ} (hn : 2 ≤ n) : 3 ≤ lucas n := by
  rw [← lucas_two]
  exact (lucas_le_iff (by norm_num) (by omega)).mpr hn

/-- Lemma 4: a prime Lucas number of even index has index `2^j`, `j ≥ 1`. -/
theorem pow_two_of_lucas_prime {k : ℕ} (hk : Even k) (hk0 : k ≠ 0) {p : ℕ} (hp : p.Prime)
    (h : lucas k = p) : ∃ j, 1 ≤ j ∧ k = 2 ^ j := by
  obtain ⟨j, m, hm, rfl⟩ := Nat.exists_eq_two_pow_mul_odd hk0
  have hj : 1 ≤ j := by
    rcases Nat.eq_zero_or_pos j with rfl | hj
    · simp only [pow_zero, one_mul] at hk
      exact absurd hm (Nat.not_odd_iff_even.mpr hk)
    · exact hj
  obtain ⟨i, rfl⟩ := hm
  have hdvd := lucas_dvd_lucas_mul_odd (2 ^ j) i
  rw [h] at hdvd
  have h2j : 2 ≤ 2 ^ j := by
    calc 2 = 2 ^ 1 := by norm_num
      _ ≤ 2 ^ j := Nat.pow_le_pow_right (by norm_num) hj
  have h3 := three_le_lucas h2j
  have hnat : (lucas (2 ^ j)).natAbs ∣ p := by
    have := Int.natAbs_dvd_natAbs.mpr hdvd
    simpa using this
  rcases hp.eq_one_or_self_of_dvd _ hnat with h1 | h1
  · omega
  · have hL : lucas (2 ^ j) = lucas (2 ^ j * (2 * i + 1)) := by
      rw [h]
      have : (0 : ℤ) ≤ lucas (2 ^ j) := by linarith
      rw [← Int.natAbs_of_nonneg this, h1]
    have hpos : 1 ≤ 2 ^ j * (2 * i + 1) := Nat.one_le_iff_ne_zero.mpr hk0
    have := lucas_inj (by omega) hpos hL
    exact ⟨j, hj, this.symm⟩

/-! ### C1: the systole trace of ±Γ(N) -/

/-- If `A ∈ SL₂(ℤ)` and `A ≡ ε I (mod N)` then `tr A = ε (2 − N² d)` for an integer `d`;
in particular `tr A ≡ ±2 (mod N²)`. -/
theorem trace_of_pmGamma {N : ℤ} (hN : N ≠ 0) {A : M2} (hdet : A.det = 1)
    (hA : PmGamma N A) :
    ∃ ε d : ℤ, (ε = 1 ∨ ε = -1) ∧ A.trace = ε * (2 - N ^ 2 * d) := by
  obtain ⟨ε, hε, h⟩ := hA
  obtain ⟨a, ha⟩ := h 0 0
  obtain ⟨b, hb⟩ := h 0 1
  obtain ⟨c, hc⟩ := h 1 0
  obtain ⟨e, he⟩ := h 1 1
  have h00 : A 0 0 = ε + N * a := by simp at ha; linarith
  have h01 : A 0 1 = N * b := by simp at hb; linarith
  have h10 : A 1 0 = N * c := by simp at hc; linarith
  have h11 : A 1 1 = ε + N * e := by simp at he; linarith
  have hε2 : ε * ε = 1 := by rcases hε with rfl | rfl <;> norm_num
  rw [Matrix.det_fin_two, h00, h01, h10, h11] at hdet
  have h1 : N * (ε * (a + e) + N * (a * e - b * c)) = 0 := by
    linear_combination hdet - hε2
  have h2 : ε * (a + e) + N * (a * e - b * c) = 0 :=
    (mul_eq_zero.mp h1).resolve_left hN
  refine ⟨ε, a * e - b * c, hε, ?_⟩
  rw [Matrix.trace_fin_two, h00, h11]
  linear_combination (ε * N) * h2 - (N * (a + e)) * hε2

/-- **C1 (lower bound).** A hyperbolic element of `±Γ(N)` has `|tr| ≥ N² − 2`. -/
theorem systole_trace_bound {N : ℤ} (hN : N ≠ 0) {A : M2} (hdet : A.det = 1)
    (hA : PmGamma N A) (hyp : 2 < |A.trace|) : N ^ 2 - 2 ≤ |A.trace| := by
  obtain ⟨ε, d, hε, htr⟩ := trace_of_pmGamma hN hdet hA
  have habs : |A.trace| = |2 - N ^ 2 * d| := by
    rw [htr, abs_mul]
    rcases hε with rfl | rfl <;> simp
  rw [habs] at hyp ⊢
  have hN2 : 0 < N ^ 2 := by positivity
  rcases lt_trichotomy d 0 with hd | hd | hd
  · exact le_abs.mpr (Or.inl (by nlinarith))
  · subst hd; simp at hyp
  · exact le_abs.mpr (Or.inr (by nlinarith))

/-- The attaining element `A₀ = [[1 − N², N], [−N, 1]]`. -/
def A0 (N : ℤ) : M2 := !![1 - N ^ 2, N; -N, 1]

theorem A0_det (N : ℤ) : (A0 N).det = 1 := by
  simp [A0, Matrix.det_fin_two]
  ring

/-- `A₀ ∈ Γ(N)` (sign `ε = +1`). -/
theorem A0_mem (N : ℤ) : ∀ i j, N ∣ A0 N i j - (1 : M2) i j := by
  intro i j
  fin_cases i <;> fin_cases j <;> simp [A0]

theorem A0_pmGamma (N : ℤ) : PmGamma N (A0 N) :=
  ⟨1, Or.inl rfl, fun i j => by simpa using A0_mem N i j⟩

theorem A0_trace (N : ℤ) (hN : 2 ≤ N ^ 2) : |(A0 N).trace| = N ^ 2 - 2 := by
  rw [show (A0 N).trace = 2 - N ^ 2 by simp [A0, Matrix.trace_fin_two]; ring]
  rw [abs_of_nonpos (by linarith)]
  ring

/-! ### The golden lift -/

/-- `n` is an admissible exponent: `n ≥ 2` even with `Q^n ∈ ±Γ(p)`. -/
def IsGoldenExp (p : ℕ) (n : ℕ) : Prop := 2 ≤ n ∧ Even n ∧ PmGamma (p : ℤ) (Q ^ n)

theorem exists_goldenExp (p : ℕ) (hp : p.Prime) : ∃ n, IsGoldenExp p n := by
  have : Fact p.Prime := ⟨hp⟩
  classical
  let φ := (Int.castRingHom (ZMod p)).mapMatrix (m := Fin 2)
  let R : Matrix (Fin 2) (Fin 2) (ZMod p) := φ Q
  have hQ : Q * (Q - 1) = 1 := by
    rw [mul_sub, mul_one, cayley_hamilton_two, Q_trace, Q_det]; simp
  have hQ' : (Q - 1) * Q = 1 := by
    rw [sub_mul, one_mul, cayley_hamilton_two, Q_trace, Q_det]; simp
  let u : (Matrix (Fin 2) (Fin 2) (ZMod p))ˣ :=
    ⟨R, φ (Q - 1), by rw [← map_mul, hQ, map_one], by rw [← map_mul, hQ', map_one]⟩
  have hpos := orderOf_pos u
  set n := 2 * orderOf u with hn
  have hRn : φ (Q ^ n) = 1 := by
    rw [map_pow]
    have : u ^ n = 1 := by rw [hn, pow_mul', pow_orderOf_eq_one, one_pow]
    have := congrArg Units.val this
    simpa [u, R] using this
  refine ⟨n, by omega, ⟨orderOf u, by omega⟩, 1, Or.inl rfl, fun i j => ?_⟩
  have hmap : φ (Q ^ n - 1) = 0 := by rw [map_sub, hRn, map_one, sub_self]
  have hij := congrFun (congrFun hmap i) j
  simp only [φ, RingHom.mapMatrix_apply, Matrix.map_apply, Matrix.zero_apply] at hij
  rw [one_mul]
  exact (ZMod.intCast_zmod_eq_zero_iff_dvd _ p).mp (by simpa using hij)

open Classical in
/-- `n(p)`: the least even `n ≥ 2` with `Q^n ∈ ±Γ(p)`. -/
noncomputable def goldenN (p : ℕ) (hp : p.Prime) : ℕ := Nat.find (exists_goldenExp p hp)

theorem goldenN_spec (p : ℕ) (hp : p.Prime) : IsGoldenExp p (goldenN p hp) := by
  classical
  unfold goldenN
  convert Nat.find_spec (exists_goldenExp p hp)

theorem goldenN_min (p : ℕ) (hp : p.Prime) {m : ℕ} (hm : IsGoldenExp p m) :
    goldenN p hp ≤ m := by
  classical
  unfold goldenN
  convert Nat.find_min' (exists_goldenExp p hp) hm

/-- The golden lift `Q^{n(p)}` lies in `SL₂(ℤ)`, in `±Γ(p)`, and is hyperbolic. -/
theorem golden_lift_props (p : ℕ) (hp : p.Prime) :
    (Q ^ goldenN p hp).det = 1 ∧ PmGamma p (Q ^ goldenN p hp) ∧
      2 < |(Q ^ goldenN p hp).trace| ∧ (Q ^ goldenN p hp).trace = lucas (goldenN p hp) := by
  obtain ⟨h2, hev, hmem⟩ := goldenN_spec p hp
  refine ⟨by rw [Q_pow_det, hev.neg_one_pow], hmem, ?_, trace_Q_pow _⟩
  rw [trace_Q_pow]
  have := three_le_lucas h2
  rw [abs_of_pos (by linarith)]
  linarith

/-- If `p = L_k` with `k` even then `Q^{2k} ≡ −I (mod p)`. -/
theorem Q_pow_two_mul_mem {k : ℕ} (hk : Even k) {N : ℤ} (hN : N = lucas k) :
    PmGamma N (Q ^ (2 * k)) := by
  have hdet : (Q ^ k).det = 1 := by rw [Q_pow_det, hk.neg_one_pow]
  have hsq : Q ^ (2 * k) = N • Q ^ k - 1 := by
    rw [two_mul, pow_add, cayley_hamilton_two, hdet, trace_Q_pow, hN, one_smul]
  refine ⟨-1, Or.inr rfl, fun i j => ⟨(Q ^ k) i j, ?_⟩⟩
  rw [hsq, Matrix.sub_apply, Matrix.smul_apply, smul_eq_mul]
  ring

/-- **C2 (trace form).** For an odd prime `p`, `L_{n(p)} = p² − 2` iff `p = L_{2^j}`, `j ≥ 1`. -/
theorem golden_systolic_iff (p : ℕ) (hp : p.Prime) (hodd : p ≠ 2) :
    lucas (goldenN p hp) = (p : ℤ) ^ 2 - 2 ↔ ∃ j, 1 ≤ j ∧ (p : ℤ) = lucas (2 ^ j) := by
  have hp3 : (3 : ℤ) ≤ p := by
    have := hp.two_le
    have : 3 ≤ p := by omega
    exact_mod_cast this
  obtain ⟨h2, hev, -⟩ := goldenN_spec p hp
  constructor
  · intro h
    obtain ⟨k, hk⟩ := hev
    have hk2 : goldenN p hp = 2 * k := by omega
    rw [hk2, lucas_two_mul] at h
    have hL := lucas_pos k
    rcases Nat.even_or_odd k with hke | hko
    · rw [hke.neg_one_pow] at h
      have hsq : (lucas k - p) * (lucas k + p) = 0 := by linear_combination h
      have hLp : lucas k = p := by
        rcases mul_eq_zero.mp hsq with h' | h'
        · linarith
        · linarith
      obtain ⟨j, hj, hkj⟩ :=
        pow_two_of_lucas_prime hke (by omega) hp hLp
      exact ⟨j, hj, by rw [← hkj, hLp]⟩
    · rw [hko.neg_one_pow] at h
      exfalso
      rcases le_or_gt (p : ℤ) (lucas k) with h' | h'
      · nlinarith
      · nlinarith
  · rintro ⟨j, hj, hpj⟩
    set k := 2 ^ j with hkdef
    have hke : Even k := (Nat.even_pow' (by omega)).mpr even_two
    have hk1 : 1 ≤ k := Nat.one_le_two_pow
    have hmem : IsGoldenExp p (2 * k) :=
      ⟨by omega, even_two_mul k, Q_pow_two_mul_mem hke hpj⟩
    have hle := goldenN_min p hp hmem
    obtain ⟨hdet, hpm, hhyp, htr⟩ := golden_lift_props p hp
    have hN : (p : ℤ) ≠ 0 := by exact_mod_cast hp.ne_zero
    have hb := systole_trace_bound hN hdet hpm hhyp
    have hL2k : lucas (2 * k) = (p : ℤ) ^ 2 - 2 := by
      rw [lucas_two_mul, hke.neg_one_pow, ← hpj]; ring
    rw [htr, abs_of_pos (by have := lucas_pos (goldenN p hp); linarith), ← hL2k] at hb
    have hge := (lucas_le_iff (by omega) (by omega)).mp hb
    have heq : goldenN p hp = 2 * k := le_antisymm hle hge
    rw [heq, hL2k]

/-- In the systolic case `n(p) = 2^{j+1}`. -/
theorem goldenN_eq (p : ℕ) (hp : p.Prime) {j : ℕ} (hj : 1 ≤ j) (hpj : (p : ℤ) = lucas (2 ^ j)) :
    goldenN p hp = 2 ^ (j + 1) := by
  have hodd : p ≠ 2 := by
    intro h2
    have h3 := three_le_lucas (show 2 ≤ 2 ^ j from
      calc 2 = 2 ^ 1 := by norm_num
        _ ≤ 2 ^ j := Nat.pow_le_pow_right (by norm_num) hj)
    rw [← hpj, h2] at h3; norm_num at h3
  have h := (golden_systolic_iff p hp hodd).mpr ⟨j, hj, hpj⟩
  have hL : lucas (2 * 2 ^ j) = (p : ℤ) ^ 2 - 2 := by
    rw [lucas_two_mul, (Nat.even_pow' (by omega)).mpr even_two |>.neg_one_pow, ← hpj]; ring
  have h2 := (goldenN_spec p hp).1
  rw [pow_succ, mul_comm]
  exact lucas_inj (by omega) (by have := Nat.one_le_two_pow (n := j); omega) (h.trans hL.symm)

/-- **C2 (extremal form).** For an odd prime `p`, the golden lift `Q^{n(p)}` has the least
`|trace|` among all hyperbolic elements of `±Γ(p)` (i.e. it represents a systole of `Γ(p)\H`,
via ℓ = 2 arccosh(|tr|/2), not formalised) iff `p = L_{2^j}` for some `j ≥ 1`. -/
theorem golden_minimal_iff (p : ℕ) (hp : p.Prime) (hodd : p ≠ 2) :
    (∀ A : M2, A.det = 1 → PmGamma p A → 2 < |A.trace| →
        |(Q ^ goldenN p hp).trace| ≤ |A.trace|) ↔
      ∃ j, 1 ≤ j ∧ (p : ℤ) = lucas (2 ^ j) := by
  rw [← golden_systolic_iff p hp hodd]
  obtain ⟨hdet, hpm, hhyp, htr⟩ := golden_lift_props p hp
  have hN : (p : ℤ) ≠ 0 := by exact_mod_cast hp.ne_zero
  have hp3 : (3 : ℤ) ≤ p := by
    have := hp.two_le
    have : 3 ≤ p := by omega
    exact_mod_cast this
  have hb := systole_trace_bound hN hdet hpm hhyp
  have hpos : 0 < lucas (goldenN p hp) := by have := lucas_pos (goldenN p hp); linarith
  rw [htr, abs_of_pos hpos] at hb
  constructor
  · intro H
    have hA := H (A0 p) (A0_det p) (A0_pmGamma p)
      (by rw [A0_trace p (by nlinarith)]; nlinarith)
    rw [htr, abs_of_pos hpos, A0_trace p (by nlinarith)] at hA
    linarith
  · intro h A hA hmem hhypA
    rw [htr, abs_of_pos hpos, h]
    exact systole_trace_bound hN hA hmem hhypA

/-! ### The stabiliser of the golden axis

A matrix of `SL₂(ℤ)` fixes `φ` (equivalently, for `A ∈ SL₂(ℤ)`, commutes with `Q`; this equivalence is
proved on paper only) iff it is `±Q^{2m}`, `m ∈ ℤ`. We prove: `A.det = 1 ∧ A * Q = Q * A →
A = ±Q^{2m}` or `A = ±(Q - 1)^{2m}` (`Q - 1 = Q⁻¹`). The unit-group step (norm-±1 elements of `ℤ[φ]`
are `±φ^n`) is done by an integer descent, with no appeal to Dirichlet's unit theorem. -/

/-- `d·I + c·Q` (the image of `d + cφ ∈ ℤ[φ]`). -/
def lin (d c : ℤ) : M2 := d • (1 : M2) + c • Q

theorem Q_mul_Q_sub_one : Q * (Q - 1) = 1 := by
  rw [mul_sub, mul_one, cayley_hamilton_two, Q_trace, Q_det]; simp

theorem Q_sub_one_mul_Q : (Q - 1) * Q = 1 := by
  rw [sub_mul, one_mul, cayley_hamilton_two, Q_trace, Q_det]; simp

theorem lin_eq (d c : ℤ) : lin d c = !![d + c, c; c, d] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp only [lin, Matrix.add_apply, Matrix.smul_apply, smul_eq_mul, Q] <;> simp

theorem lin_det (d c : ℤ) : (lin d c).det = d ^ 2 + d * c - c ^ 2 := by
  rw [lin_eq]; simp [Matrix.det_fin_two]; ring

theorem lin_eq_Q_mul (d c : ℤ) : lin d c = Q * lin (c - d) d := by
  unfold lin
  rw [mul_add, mul_smul_comm, mul_one, mul_smul_comm, cayley_hamilton_two, Q_trace, Q_det]
  module

theorem lin_eq_Q_sub_one_mul (d c : ℤ) : lin d c = (Q - 1) * lin c (c + d) := by
  unfold lin
  rw [mul_add, mul_smul_comm, mul_one, mul_smul_comm, sub_mul, one_mul, cayley_hamilton_two,
    Q_trace, Q_det]
  module

/-- `A = ±Q^n` or `A = ±(Q - 1)^n` for some `n : ℕ`, i.e. `A = ±Q^z`, `z ∈ ℤ`. -/
def IsPmQPow (A : M2) : Prop :=
  ∃ s : ℤ, (s = 1 ∨ s = -1) ∧ ∃ n : ℕ, A = s • Q ^ n ∨ A = s • (Q - 1) ^ n

theorem IsPmQPow.Q_mul {A : M2} (h : IsPmQPow A) : IsPmQPow (Q * A) := by
  obtain ⟨s, hs, n, h | h⟩ := h <;> subst h
  · exact ⟨s, hs, n + 1, Or.inl (by rw [mul_smul_comm, pow_succ'])⟩
  · cases n with
    | zero => exact ⟨s, hs, 1, Or.inl (by rw [mul_smul_comm, pow_zero, mul_one, pow_one])⟩
    | succ n =>
      refine ⟨s, hs, n, Or.inr ?_⟩
      rw [mul_smul_comm, pow_succ', ← mul_assoc, Q_mul_Q_sub_one, one_mul]

theorem IsPmQPow.Q_sub_one_mul {A : M2} (h : IsPmQPow A) : IsPmQPow ((Q - 1) * A) := by
  obtain ⟨s, hs, n, h | h⟩ := h <;> subst h
  · cases n with
    | zero =>
      exact ⟨s, hs, 1, Or.inr (by rw [mul_smul_comm, pow_zero, mul_one, pow_one])⟩
    | succ n =>
      refine ⟨s, hs, n, Or.inl ?_⟩
      rw [mul_smul_comm, pow_succ', ← mul_assoc, Q_sub_one_mul_Q, one_mul]
  · exact ⟨s, hs, n + 1, Or.inr (by rw [mul_smul_comm, pow_succ'])⟩

/-- Integer descent: if `d² + dc − c² = ±1` then `d·I + c·Q = ±Q^z`. -/
theorem lin_isPmQPow_aux : ∀ n : ℕ, ∀ d c : ℤ, d ^ 2 + c ^ 2 ≤ n →
    (d ^ 2 + d * c - c ^ 2 = 1 ∨ d ^ 2 + d * c - c ^ 2 = -1) → IsPmQPow (lin d c) := by
  intro n
  induction n with
  | zero =>
    intro d c hle hN
    push_cast at hle
    have hd : d = 0 := by nlinarith [sq_nonneg d, sq_nonneg c]
    have hc : c = 0 := by nlinarith [sq_nonneg d, sq_nonneg c]
    subst hd hc
    norm_num at hN
  | succ n ih =>
    intro d c hle hN
    by_cases hc : c = 0
    · subst hc
      have hd : d = 1 ∨ d = -1 := by
        have : (d - 1) * (d + 1) = 0 := by rcases hN with h | h <;> nlinarith [sq_nonneg d]
        rcases mul_eq_zero.mp this with h | h
        · left; linarith
        · right; linarith
      exact ⟨d, hd, 0, Or.inl (by simp [lin])⟩
    by_cases hd : d = 0
    · subst hd
      have hc' : c = 1 ∨ c = -1 := by
        have : (c - 1) * (c + 1) = 0 := by rcases hN with h | h <;> nlinarith [sq_nonneg c]
        rcases mul_eq_zero.mp this with h | h
        · left; linarith
        · right; linarith
      exact ⟨c, hc', 1, Or.inl (by simp [lin])⟩
    rcases lt_or_gt_of_ne (mul_ne_zero hd hc) with hdc | hdc
    · -- opposite signs: multiply by φ, (d, c) ↦ (c, c + d)
      have hlt : c ^ 2 + (c + d) ^ 2 < d ^ 2 + c ^ 2 := by
        by_contra hcon
        rcases lt_or_gt_of_ne hd with hd' | hd' <;> rcases lt_or_gt_of_ne hc with hc' | hc' <;>
          rcases hN with h | h <;> nlinarith
      have hN' : c ^ 2 + c * (c + d) - (c + d) ^ 2 = 1 ∨ c ^ 2 + c * (c + d) - (c + d) ^ 2 = -1 := by
        rcases hN with h | h
        · right; linear_combination -h
        · left; linear_combination -h
      have := ih c (c + d) (by push_cast at hle ⊢; linarith) hN'
      rw [lin_eq_Q_sub_one_mul]
      exact this.Q_sub_one_mul
    · -- same sign: multiply by φ⁻¹, (d, c) ↦ (c − d, d)
      have hlt : (c - d) ^ 2 + d ^ 2 < d ^ 2 + c ^ 2 := by
        by_contra hcon
        rcases lt_or_gt_of_ne hd with hd' | hd' <;> rcases lt_or_gt_of_ne hc with hc' | hc' <;>
          rcases hN with h | h <;> nlinarith
      have hN' : (c - d) ^ 2 + (c - d) * d - d ^ 2 = 1 ∨ (c - d) ^ 2 + (c - d) * d - d ^ 2 = -1 := by
        rcases hN with h | h
        · right; linear_combination -h
        · left; linear_combination -h
      have := ih (c - d) d (by push_cast at hle ⊢; linarith) hN'
      rw [lin_eq_Q_mul]
      exact this.Q_mul

theorem lin_isPmQPow {d c : ℤ} (hN : d ^ 2 + d * c - c ^ 2 = 1 ∨ d ^ 2 + d * c - c ^ 2 = -1) :
    IsPmQPow (lin d c) :=
  lin_isPmQPow_aux (d ^ 2 + c ^ 2).toNat d c (by rw [Int.toNat_of_nonneg (by positivity)]) hN

theorem eq_lin_of_commute {A : M2} (h : A * Q = Q * A) : A = lin (A 1 1) (A 0 1) := by
  have e := fun i j => congrFun (congrFun h i) j
  have e00 := e 0 0
  have e01 := e 0 1
  simp [Q, Matrix.mul_apply, Fin.sum_univ_two] at e00 e01
  ext i j
  rw [lin_eq]
  fin_cases i <;> fin_cases j <;> simp <;> linarith

/-- **Stabiliser of the golden axis.** If `A ∈ SL₂(ℤ)` commutes with `Q` then `A = ±Q^{2m}` or
`A = ±Q^{-2m}` (`Q⁻¹ = Q - 1`) for some `m : ℕ`. -/
theorem stab_golden {A : M2} (hdet : A.det = 1) (h : A * Q = Q * A) :
    ∃ s : ℤ, (s = 1 ∨ s = -1) ∧ ∃ m : ℕ, A = s • Q ^ (2 * m) ∨ A = s • (Q - 1) ^ (2 * m) := by
  have hA := eq_lin_of_commute h
  have hN : (A 1 1) ^ 2 + A 1 1 * A 0 1 - (A 0 1) ^ 2 = 1 := by
    rw [← lin_det, ← hA, hdet]
  obtain ⟨s, hs, n, hn⟩ := lin_isPmQPow (Or.inl hN)
  rw [← hA] at hn
  have hs2 : s ^ 2 = 1 := by rcases hs with rfl | rfl <;> norm_num
  have hQ1 : (Q - 1).det = -1 := by simp [Q, Matrix.det_fin_two]
  have hev : Even n := by
    have key : ((-1 : ℤ)) ^ n = 1 := by
      rcases hn with hn | hn
      · have := congrArg Matrix.det hn
        rw [Matrix.det_smul, Matrix.det_pow, Q_det, hdet] at this
        simp only [Fintype.card_fin] at this
        rw [hs2, one_mul] at this
        exact this.symm
      · have := congrArg Matrix.det hn
        rw [Matrix.det_smul, Matrix.det_pow, hQ1, hdet] at this
        simp only [Fintype.card_fin] at this
        rw [hs2, one_mul] at this
        exact this.symm
    exact (neg_one_pow_eq_one_iff_even (by norm_num)).mp key
  obtain ⟨m, hm⟩ := hev
  exact ⟨s, hs, m, by rw [show 2 * m = n by omega]; exact hn⟩

end BolzaKit.GoldenGeodesic
