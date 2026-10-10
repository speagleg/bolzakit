/-
Bolza surface: the defining relation of the Bolza Fuchsian group, checked exactly.

Model (ported from the original research code).  Let ζ = e^{iπ/4}, μ = √(2+2√2).
A group element is the SU(1,1) matrix [[A, μB], [μ·conj B, conj A]] with A, B ∈ ℤ[ζ]; it is stored as
the pair (A, B).  Since μ² = 2 + 2√2 ∈ ℤ[ζ], the product of such matrices is again of this form:
  (A₁,B₁)(A₂,B₂) = (A₁A₂ + μ² B₁ conj B₂,  A₁B₂ + B₁ conj A₂).
ℤ[ζ] = ℤ[x]/(x⁴+1) is represented by 4 integer coefficients of 1, ζ, ζ², ζ³.

The generators g_k (k = 0..7, g_{k+4} = g_k⁻¹) are the hyperbolic translations of length
2·arccosh(1+√2) along the axes at angles kπ/4: g_k = (1+√2, ζ^k).  We prove
  * g_k · g_{k+4} = 1 for all k, and
  * g₀ g₁⁻¹ g₂ g₃⁻¹ g₀⁻¹ g₁ g₂⁻¹ g₃ = 1
as identities in ℤ[ζ]², by kernel evaluation (`decide`), no placeholder proofs, no compiled-code evaluation.
Also: det g_k = |A|² − μ²|B|² = 1 and the half-trace Re A = 1 + √2.
-/
import Mathlib

namespace BolzaKit.Group

/-- ℤ[ζ₈] with ζ⁴ = -1: a₀ + a₁ζ + a₂ζ² + a₃ζ³. -/
structure Z8 where
  a0 : Int
  a1 : Int
  a2 : Int
  a3 : Int
  deriving DecidableEq, Repr

namespace Z8

def add (x y : Z8) : Z8 := ⟨x.a0 + y.a0, x.a1 + y.a1, x.a2 + y.a2, x.a3 + y.a3⟩

def neg (x : Z8) : Z8 := ⟨-x.a0, -x.a1, -x.a2, -x.a3⟩

/-- negacyclic convolution (ζ⁴ = -1) -/
def mul (x y : Z8) : Z8 :=
  ⟨x.a0 * y.a0 - x.a1 * y.a3 - x.a2 * y.a2 - x.a3 * y.a1,
   x.a0 * y.a1 + x.a1 * y.a0 - x.a2 * y.a3 - x.a3 * y.a2,
   x.a0 * y.a2 + x.a1 * y.a1 + x.a2 * y.a0 - x.a3 * y.a3,
   x.a0 * y.a3 + x.a1 * y.a2 + x.a2 * y.a1 + x.a3 * y.a0⟩

/-- complex conjugation ζᵏ ↦ ζ⁻ᵏ = -ζ^{4-k} -/
def conj (x : Z8) : Z8 := ⟨x.a0, -x.a3, -x.a2, -x.a1⟩

def zero : Z8 := ⟨0, 0, 0, 0⟩
def one : Z8 := ⟨1, 0, 0, 0⟩
def zeta : Z8 := ⟨0, 1, 0, 0⟩
/-- √2 = ζ - ζ³ -/
def sqrt2 : Z8 := ⟨0, 1, 0, -1⟩
/-- μ² = 2 + 2√2 -/
def mu2 : Z8 := ⟨2, 2, 0, -2⟩

end Z8

/-- An element [[A, μB], [μ conj B, conj A]] of SU(1,1), stored as (A, B). -/
structure El where
  A : Z8
  B : Z8
  deriving DecidableEq, Repr

namespace El
open Z8

def mul (g h : El) : El :=
  ⟨Z8.add (Z8.mul g.A h.A) (Z8.mul mu2 (Z8.mul g.B (Z8.conj h.B))),
   Z8.add (Z8.mul g.A h.B) (Z8.mul g.B (Z8.conj h.A))⟩

def one : El := ⟨Z8.one, Z8.zero⟩

/-- determinant |A|² − μ²|B|² (an element of ℤ[√2] ⊂ ℤ[ζ]) -/
def det (g : El) : Z8 :=
  Z8.add (Z8.mul g.A (Z8.conj g.A)) (Z8.neg (Z8.mul mu2 (Z8.mul g.B (Z8.conj g.B))))

/-- 1 + √2 -/
def aGen : Z8 := ⟨1, 1, 0, -1⟩

/-- ζᵏ for k = 0..7 -/
def zpow : Fin 8 → Z8
  | 0 => ⟨1, 0, 0, 0⟩ | 1 => ⟨0, 1, 0, 0⟩ | 2 => ⟨0, 0, 1, 0⟩ | 3 => ⟨0, 0, 0, 1⟩
  | 4 => ⟨-1, 0, 0, 0⟩ | 5 => ⟨0, -1, 0, 0⟩ | 6 => ⟨0, 0, -1, 0⟩ | 7 => ⟨0, 0, 0, -1⟩

/-- the side-pairing generator g_k = (1+√2, ζᵏ); g_{k+4} = g_k⁻¹ -/
def gen (k : Fin 8) : El := ⟨aGen, zpow k⟩

def word : List (Fin 8) → El
  | [] => one
  | k :: ks => mul (gen k) (word ks)

end El

open El

/-- every generator has determinant 1 -/
theorem gen_det (k : Fin 8) : (gen k).det = Z8.one := by
  revert k; decide

/-- g_k · g_{k+4} = 1 -/
theorem gen_inv (k : Fin 4) :
    El.mul (gen ⟨k.val, by omega⟩) (gen ⟨k.val + 4, by omega⟩) = El.one := by
  revert k; decide

/-- half-trace: 2·Re A = A + conj A = 2(1 + √2) -/
theorem gen_trace (k : Fin 8) :
    Z8.add (gen k).A (Z8.conj (gen k).A) = ⟨2, 2, 0, -2⟩ := by
  revert k; decide

/-- The Bolza surface relation g₀ g₁⁻¹ g₂ g₃⁻¹ g₀⁻¹ g₁ g₂⁻¹ g₃ = 1, exactly. -/
theorem bolza_relation : word [0, 5, 2, 7, 4, 1, 6, 3] = El.one := by
  decide

/-- The adjacent-side pairing relation abABcdCD (used by the broken source code) does NOT hold for
these generators. -/
theorem adjacent_relation_fails : word [0, 1, 4, 5, 2, 3, 6, 7] ≠ El.one := by
  decide

end BolzaKit.Group
