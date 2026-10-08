import Mathlib

/-! Salem-type norm identity: the Bolza quartic is the norm from Q(sqrt 2) of the systolic characteristic polynomial. -/

theorem item15_bolza_norm (x : ℝ) :
    (x ^ 2 - (2 + 2 * Real.sqrt 2) * x + 1) * (x ^ 2 - (2 - 2 * Real.sqrt 2) * x + 1)
      = x ^ 4 - 4 * x ^ 3 - 2 * x ^ 2 - 4 * x + 1 := by
  have h : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  linear_combination (-4 * x ^ 2) * h
