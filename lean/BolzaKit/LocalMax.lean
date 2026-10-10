/-
Root of the LocalMax skeleton (item 23) (Bolza λ₁ local maximality).  See proofs/bolza-local-max/LEAN.md.
Building this file prints the axioms used by the main theorems.
-/
import BolzaKit.LocalMax.RepTheory
import BolzaKit.LocalMax.LocalMax
import BolzaKit.LocalMax.Taylor
import BolzaKit.LocalMax.RepReal
import BolzaKit.LocalMax.Certificates
import BolzaKit.LocalMax.KyFan
import BolzaKit.LocalMax.Assembly
import BolzaKit.LocalMax.Model

-- main theorem, conditional on the explicit hypothesis structure `BolzaInputs`
#print axioms Bolza.bolza_strict_local_max
-- Lemma 10 (expansion form and C² form)
#print axioms Bolza.local_max_criterion
#print axioms Bolza.local_max_criterion_C2
#print axioms Bolza.taylor2_littleO
-- derived steps of the assembly
#print axioms Bolza.trace_A_zero
#print axioms Bolza.A_vanishes_K0
#print axioms Bolza.kerA_sub_K0
#print axioms Bolza.hess_K0_neg
#print axioms Bolza.wB_sign
-- representation theory (kernel-evaluated)
#print axioms Bolza.Rep.InH_iff
#print axioms Bolza.Rep.chiT_norm
#print axioms Bolza.Rep.sym0_O_K0
#print axioms Bolza.Rep.smid_line_in_Kminus
#print axioms Bolza.Rep.conj_t
#print axioms Bolza.invariant_form_Kminus
#print axioms Bolza.Wplus_irreducible
-- Ky Fan
#print axioms Bolza.second_deriv_le_of_le
#print axioms Bolza.ky_fan
-- certificate arithmetic
#print axioms Bolza.Cert.sign_joint
#print axioms Bolza.Cert.sign_su13_implies_lg
#print axioms Bolza.Cert.sign_rows_not_additive
#print axioms Bolza.Cert.pairing_trace
-- non-vacuity model
#print axioms Bolza.Model.toyInputs
