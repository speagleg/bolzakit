# NOVELTY: is "Bolza is a strict local maximum of λ₁ on T₂" already known?

Search date 2026-10-08. This file extends `work/LITERATURE.md` and does not repeat it: quotes already recorded there
(SU13 §5, Str16 pp. 22–23, KMP Thm 1.1 and §4.5, Bonifacio, FBP-K, FBP-LP, Cook Conj. 5.1, Petri HDR, Kalvin 2026)
are only cross-referenced here. Every new quote below was read in the source text (arXiv PDF → pdftotext) unless it
is marked otherwise. Line breaks have been removed from the quotes.

## 1. Verdict

**NEW as a theorem; the statement itself was already asserted numerically.**

- **No proof found.** No paper, preprint or thesis we could reach proves that B is a local maximum (strict or not)
  of λ₁ on T₂ or M₂. None proves anything that implies it, at either the hyperbolic or the all-metrics level.
- **Global maximality is still open.** As of 2026-10 it is still described as conjectural:
  - Hide–Petri, Compositio 2026 / arXiv:2412.15111 (Dec 2024);
  - FBP-LP v3 (Feb 2026);
  - Bonifacio v2 (May 2025).
  So local maximality cannot be read off from a global theorem.
- **Prior numerical assertions.** SU13 and Str16 assert local maximality on numerical grounds. Str16 says it most
  explicitly: "The reason for it being locally maximal is however its degeneracy". No analysis was published; SU13's
  "the corresponding analysis will be discussed elsewhere" never appeared that we could find. Cook 2018 proposes a
  route but proves nothing.
- **What a paper should claim.** It should claim the *first proof*, computer-assisted, of a numerically observed
  fact. It should not claim to discover the fact.
- **Residual risk (UNCLEAR parts), listed in §6.**
  - Maintrot's EPFL thesis (2012, advisor Buser) on FEM for the Laplacian on Riemann surfaces could not be opened
    (WAF block).
  - The GEMSTONE minicourse notes by Fortier Bourque (Laval) could not be opened (DNS failure).
  - Neither is likely to contain a proof, but either could contain a numerical Teichmüller study or a remark.

## 2. Provenance of the conjecture

Nobody stated *local* maximality as a named conjecture. The literature conjectures *global* maximality. Local
maximality appears only as a numerical assertion (Str16), an unpublished "analysis" (SU13), and a proposed proof step
(Cook).

| Source | What it says (verbatim) | Kind |
|---|---|---|
| Jenni 1984 (Comment. Math. Helv. 59) | not read; SU13 says Jenni estimated λ₁ ∈ [3.83, 3.85] | — |
| SU13, arXiv:1110.2150 §5 | "we conjecture that this is the maximal value of the first eigenvalue … in genus 2. This conjecture is supported by our calculations in Teichmüller space (the corresponding analysis will be discussed elsewhere)" (LITERATURE.md) | **origin of the global conjecture** |
| Str16, arXiv:1604.02722 pp. 22–23 | "Numerical evidence suggests that this is the global maximum … The reason for it being locally maximal is however its degeneracy. … Numerically one can see that no matter what perturbation one chooses, none of the eigenvalues λ1, λ2 and λ3 has an extremal value at the Bolza surface." (re-read this session, lines 1392–1398 of the text) | **numerical assertion of local max** |
| Cook thesis 2018, arXiv:2108.11825 Ch. 5 | "Conjecture 5.1 Among all compact Riemann surfaces of genus 2 with constant negative curvature, the Bolza surface globally maximizes the first positive eigenvalue … It is thought that by analyzing the nodal lines of eigenfunctions relating to χ8 under perturbations in Teichmüller space, one would be able to show that the Bolza surface is a local maximum for λ1. One would then be able to use this as a building block for proving that it is a global maximum." (re-read, p. 98) | global conjecture + local max as a proposed step |
| Bonifacio 2107.09674 | "For a hyperbolic surface with g = 2, the largest value of λ1 is conjectured to be that of the Bolza surface" | restates SU13 |
| Bonifacio 2111.13215v2 (30 May 2025) | "this provides strong evidence for the conjecture that the Bolza surface is the closed hyperbolic surface with the largest spectral gap [4]" (§4) | restates SU13 |
| KMP 2111.12716v3 §4.5 | "giving very strong evidence for the conjecture that the Bolza surface maximizes λ1 among all g = 2 orbifolds" | restates SU13 (extended to orbifolds) |
| FBP-K 2111.14699 §1 | "That B maximizes λ1 among hyperbolic surfaces was conjectured in [SU13]" | restates |
| FBP-LP 2302.02540v3 §9 | "these surfaces are the conjectured maximizers in these genera"; "it is not known if the supremum of λ1 is attained"; "Our intuition is that m1 should be maximized among local maximizers of λ1." | restates; mentions local maximisers in general but makes no claim about B |
| Lee 2311.02632 §1 (Fricke–Macbeath) | "In a work by Strohmaier and Uski [SU13], they conjectured that the Bolza surface maximizes λ1 among all hyperbolic surfaces of genus 2 based on numerical experiments. A recent study [KMP21] provides additional evidence supporting the conjecture" | restates |
| Hide–Petri 2412.15111 §1 | "Λ2 = sup { λ1 (X) : X a closed hyperbolic surface of genus 2 }. The latter conjecturally equals the spectral gap of the Bolza surface"; also "There are however isolated points outside of this interval, for instance corresponding to the Bolza surface and the Klein quartic" (this is about isolated points of the *arithmetic* bass-note set, not local maximality on T₂) | restates |
| Petri HDR | "it's natural to conjecture that these two surfaces maximize the first eigenvalue in their moduli spaces" (LITERATURE.md) | restates |
| Mondal 1406.1080 | "It is known that the Bolza surface has λ1 approximately 3.8 (see [S-U] for more details)"; proves {λ₁ > 1/4} is unbounded and disconnects M₂ | context only |
| Strohmaier, private communication (via FBP-K) | conjectures B uniquely maximises m₁ in genus 2 | m₁, not λ₁ |

## 3. Has anything since 2023 settled it (globally or locally)?

**No.** The relevant post-2023 work:

- **Radcliffe, *Non-saturation of bootstrap bounds by hyperbolic orbifolds*, arXiv:2404.14479 (v2, Nov 2025).**
  This is the most relevant new item. It is directly about whether the KMP/Bonifacio genus-2 bound can be closed.
  - §3.3: "The bound on λ1 for genus two surfaces (2.12) is almost saturated by λ1 on the Bolza surface … the
    Bolza surface does not saturate λ̃1^∞", with "λ̃1^∞ − λ1 ≥ Δ∞ ≈ 5.59683 × 10−6".
  - §4: "the bounds on λ1, (1.1) and (2.12) are not saturated by the [0; 2, 3, 7] orbifold or the Bolza surface
    respectively, under reasonable assumptions about the convergence of the numerics."
  - It adds that richer correlator systems "might give tighter bounds".
  - **Consequence.** The idea "KMP bound + only Bolza nearly attains it ⇒ Bolza is the global max" *cannot* work
    with the existing single-correlator bootstrap. Even at infinite derivative order, the bound stays about 5.6e-6
    above λ₁(B) (numerically, not rigorously). So the global-max route is not closed, and local maximality does
    not follow from any bound.
- **KMP v3 (Jan 2024; Comm. AMS).** Unchanged: rigorous λ₁ < 3.8388976481 for all genus-2 orbifolds, 1.04e-5 above
  λ₁(B). Only "very strong evidence".
- **Bonifacio 2111.13215v2 (May 2025).** Unchanged: numerical λ₁ ≤ 3.8388977; "evidence".
- **Gesteau–Pal–Simmons-Duffin–Xu, arXiv:2311.13330.** Spin bootstrap. "the Bolza surface nearly saturate[s] the
  bounds at genus … 2"; §5: "the most symmetric surfaces … sit almost exactly at the boundary, or even the kink, of
  the allowed region." This is evidence only.
- **Adve, arXiv:2509.17935.** Converse theorem: solutions of the bootstrap equations come from hyperbolic surfaces.
  It is structural and has no λ₁ extremality result.
- **FBP-LP v3 (Feb 2026); FB–Petri et al., *Three counterexamples …*, arXiv:2312.03504v2 (Jun 2026).**
  Multiplicity and LP bounds only. They contain no local-max statement about λ₁(B).
- **Hatcher–Petri, arXiv:2610.07193 (5 Oct 2026); Hide–Petri 2412.15111; Hide–Petri–Roig Sanchis–Thomas
  2512.13139; Monk–Naud survey 2601.13988; Monk 2601.15157.** These are about covers and large genus; Bolza appears
  only as a base surface or example.
- **Kalvin, arXiv:2608.01611 (Aug 2026).** Critical point of det Δ only (in LITERATURE.md).

**Has the global maximum been proved? No.**

## 4. Why Nayatani–Shoda (and the conformal / all-metrics theory) does not bear on this question

- **Correction to the task statement.** The sharp genus-2 all-metrics constant is **16π**, not 8π. 8π is Hersch's
  sphere constant.
  - Nayatani–Shoda, arXiv:1704.06384 (CRAS 2019), abstract: "we settle in the affirmative the
    Jakobson-Levitin-Nadirashvili-Nigam-Polterovich conjecture, stating that a certain singular metric on the Bolza
    surface, with area normalized, should maximize the first eigenvalue of the Laplacian."
  - Their introduction: "When γ = 2, the inequality (1) becomes Λ(ds²) ≤ 16π", with Λ = λ₁·Area.
  - JLNNP (math/0509398) abstract: the maximiser is "a singular surface which is realized as a double branched
    covering over a sphere. The six ramification points are chosen in such a way that this surface has a complex
    structure of the Bolza surface."
- **Reasons this is irrelevant to the hyperbolic local problem:**
  1. **Different functional, different maximiser.**
     - NS maximise λ̄₁ = λ₁·Area over *all* metrics. The maximiser is a branched spherical (cone) metric,
       conformal to Bolza but not hyperbolic.
     - The hyperbolic Bolza metric has λ̄₁ = 4π·3.8388873 ≈ 48.241 < 16π ≈ 50.265. It is not the maximiser in
       the full space, and not even in its own conformal class: Λ₁([g_B]) = 16π is attained by the NS/JLNNP
       metric.
  2. **Global max on a big space says nothing about local maxima on a subspace.**
     - Hyperbolic metrics, normalised to area 4π, form a 6-dimensional slice T₂ inside the space of all
       metrics.
     - The NS maximiser lies *off* this slice. So the restriction of λ̄₁ to T₂ can have its own maxima, and
       "B maximises on the big space" does not apply, because B does not.
     - The only consequence for hyperbolic surfaces is λ₁ ≤ 4 (FBP-K §6, quoted in LITERATURE.md). That is far
       above 3.8389 and useless locally.
  3. **The critical-point theory is of a different kind.**
     - Nadirashvili / El Soufi–Ilias criticality (minimal immersions or harmonic maps to spheres by first
       eigenfunctions) is criticality with respect to *all* or *conformal* metric variations.
     - Our problem has the constraint K ≡ −1. Admissible variations are trace-free divergence-free (TT) tensors,
       i.e. holomorphic quadratic differentials. Pointwise conformal factors are excluded.
     - In the Teichmüller setting, criticality of a symmetric cluster is forced by symmetry: B/Aut(B) is a rigid
       triangle orbifold, so there is no Aut-invariant quadratic differential. Kalvin uses the same argument for
       det Δ, and Bavard uses it for eutaxy of systoles. That is our tr A ≡ 0.
     - No harmonic-map condition enters. We found no claim that the hyperbolic Bolza metric is El Soufi–Ilias
       critical in its conformal class. It would not imply anything here anyway, since conformal and
       Teichmüller directions are complementary.
  4. **The second-order theory for the all-metrics problem only appeared in 2025, and it is not about B.**
     - Karpukhin, *Second variation formula for eigenvalue functionals on surfaces*, arXiv:2508.10465 (Aug 2025),
       abstract: the Nadirashvili / El Soufi–Ilias results "only contain first order information and cannot be
       used to determine whether a given critical metric a local maximiser or not. In the present paper we write
       down the second variation formula for critical metrics and show that the flat metric on the non-rhombic
       torus can never be a conformal maximiser".
     - §1 ends: "We hope that the present paper will inspire other researchers … to proving local maximality
       properties of the critical metrics."
     - Even in that (unconstrained / conformal) setting, second-order local-max proofs are therefore recent and
       rare. It treats flat tori and Steklov problems, not hyperbolic genus 2.
  5. **KPS25 (arXiv:2505.05293)** proves existence of λ̄₁-maximisers over all metrics. It is the same
     non-hyperbolic problem (LITERATURE.md).

## 5. Closest prior results and how ours differs

| Result | Object / functional | Method | Relation to our claim |
|---|---|---|---|
| SU13 §5, Str16 pp. 22–23 | λ₁ on T₂ at B | numerical perturbations (MPS) in Teichmüller space; unpublished | Same statement, numerical only. Str16's mechanism ("no matter what perturbation … none of λ1, λ2, λ3 has an extremal value") is incomplete: on the 3-dimensional kernel K (e.g. F_B), all three branches are first-order stationary and λ₁ has a strict quadratic maximum (README §3). |
| SU13 §7.1, Str16 (lines 2034–2061) | det Δ at B | numerical: "the Bolza surface is indeed a local maximum for the spectral determinant. Note that the Bolza surface is known to be a critical point by symmetry considerations." | Analogous numerical local max of a *different*, smooth spectral functional |
| Kalvin 2608.01611 | det Δ | symmetry ⇒ critical point for every quasiplatonic surface | First-order only, smooth functional. Our tr A ≡ 0 is the λ₁-cluster analogue; the second-order part is new. |
| Schmutz 1993 (GAFA; not read directly, cited via Str16: "the unique global maximum of the length of the systole. This was shown by Schmutz-Schaller") | hyperbolic systole | Voronoi-type perfect + eutactic criterion for a max of finitely many length functions | Systole is a min of *smooth* functions, so first-order (perfection + eutaxy) suffices. λ₁ is a min of eigenvalue branches with a nonzero kernel K of the first-order map, so first order *cannot* suffice. A genuine second-order certificate on K is needed (our Lemmas 8–11). |
| Fortier Bourque–Martínez-Granado–Vargas Pallete, arXiv:2105.03871 (AHL 2024) | extremal-length systole | "We prove that the extremal length systole of genus two surfaces attains a strict local maximum at the Bolza surface"; Bavard's generalised-systole framework: "any isolated fixed point of G is eutactic", plus perfection by an injectivity computation (§7.2) | Closest *structural* analogue: symmetry gives eutaxy, a computer-aided rank computation gives perfection, together a strict local max. Ours is the same pattern (symmetry ⇒ tr A ≡ 0; certified rank A = 3) but must add a second-order step because rank A = 3 < 6 = dim T₂ (no "perfection"). |
| Fortier Bourque–Rafi, *Local maxima of the systole function* (JEMS) | hyperbolic systole | eutaxy/perfection | Not about λ₁ |
| Berger 1973 / Nadirashvili 1996; Kao–Lai–Osting 1405.4944 | λ₁ (λ_k) on flat tori | explicit lattice computation; "among flat tori, Λk has a local maximum" (KLO §3.3) | Genus-1 analogue, where the eigenfunctions are explicit. No hyperbolic analogue exists. |
| KMP / Bonifacio / GPSX / Radcliffe | sup λ₁ over genus-2 orbifolds | bootstrap SDP | Upper bound only, 1.04e-5 above λ₁(B). Radcliffe shows this route alone will not reach λ₁(B). It gives no local information. |
| FBP-K Thm 7.1 | m₁(B) = 3 | trace-formula LP + representation theory | An input to our Lemma 1, not a competitor |

**What is genuinely new in our work:**
1. A proof (computer-assisted, certified) that B is a strict local maximum of λ₁ on T₂.
2. The structure behind it:
   - the first-order cluster map A has rank 3;
   - K = ker A is a 3-dimensional absolutely irreducible Isom(B)-module;
   - λ₁ decreases linearly off K and quadratically on K;
   - D²Λ|_K = c·g_inv with c < 0.
3. A correction to the only published description of the mechanism (Str16).

We found no prior second-order (Hessian on a degenerate cluster) analysis of any Laplace eigenvalue on Teichmüller
space at a symmetric surface.

## 6. What a human should still check

1. **Maintrot, *Finite Element Method on Riemann Surfaces and Applications to the Laplacian Spectrum*, EPFL thesis
   TH5337 (2012, advisor Buser).** The infoscience PDF returned "WAF: bot suspect activity" and the item page
   returned 429. It may contain genus-2 Teichmüller scans of λ₁. Check whether it states local maximality or
   computes derivatives at B.
2. **GEMSTONE minicourse 5 (Laval; Fortier Bourque, "extremal problems on hyperbolic surfaces"), Lectures 1–5 at
   archimede.mat.ulaval.ca/agirouard/GEMSTONE/MiniCourse5/Lectures/.** The host would not resolve (ENOTFOUND).
   Lecture notes on extremal problems are a plausible place for a remark such as "B is known/believed to be a local
   max of λ₁".
3. **Jenni 1984 (Comment. Math. Helv. 59, 193–203).** Not read (paywalled / no open copy). It is unlikely to
   contain a deformation result, but this should be confirmed.
4. **Schmutz 1993 GAFA "Local maxima of the systole function"** and **Bavard 2005**. These were cited for the
   method comparison only and not read directly.
5. **Strohmaier and Uski's promised "analysis … discussed elsewhere" (SU13).** We found no such publication; Uski
   has no later arXiv paper on this. It is worth emailing A. Strohmaier (Leeds) to ask whether an unpublished
   Teichmüller analysis exists, and to tell him about the K-direction behaviour that contradicts Str16's wording.
6. **Cadavid–Osorno–Ruiz 2012 (Sci. J. Physics, doi 10.7237/sjp/128)**, on critical points of λ₁-eigenfunctions in
   genus 2. Known only from the repository abstract and assumed irrelevant.
7. **MathSciNet / zbMATH full-text search** for "Bolza" + "first eigenvalue" + "local" (not available to us). Also
   Google Scholar "cited by" for SU13 to catch theses. Semantic Scholar's list was used instead (see §7).
8. **Ask the experts directly** before claiming priority: Fortier Bourque, Petri, Strohmaier, Karpukhin.

## 7. Search log (this session; LITERATURE.md has the earlier ones)

**WebSearch:**
1. `Bolza surface local maximum first eigenvalue Laplacian Teichmüller space` (extended)
2. `"Bolza surface" maximizes first eigenvalue genus 2 proof 2025`
3. `"extremal" λ1 hyperbolic surfaces Teichmüller "critical point" first eigenvalue Fortier Bourque Petri`
4. `critical points Laplace eigenvalue hyperbolic metrics Teichmüller space holomorphic quadratic differential eigenfunctions extremal` (extended)
5. `arXiv 2024 2025 "local maximum" "first eigenvalue" hyperbolic surface genus two Bolza Laplacian`
6. `"hyperbolic metrics" eigenvalue functional critical "Teichmüller space" first variation eigenvalue Hopf differential eigenfunctions "extremal" surfaces constant curvature` (extended)
7. `Schmutz Schaller Bolza surface strict local maximum systole perfect eutactic Teichmüller`
8. `Bolza surface kissing number maximum genus 2 Fortier Bourque Petri; Bolza local maximum Arakelov invariant / Bergman / Faltings delta genus 2`
9. `"Bolza surface" "local maximum" λ1 OR "spectral gap" OR "first eigenvalue" proof computer-assisted 2026` (extended)
10. `first eigenvalue Laplacian hyperbolic surfaces moduli space "not differentiable" multiplicity local maximum symmetric surfaces critical in the sense of Clarke`
11. `Uski thesis "Eigenvalues and Spectral Determinants on Compact Hyperbolic Surfaces" Loughborough Teichmüller Bolza`
12. `"On the Critical Point Structure of Eigenfunctions Belonging to the First Nonzero Eigenvalue of a Genus Two Closed Hyperbolic Surface"`
13. `"Finite Element Method on Riemann Surfaces and Applications to the Laplacian Spectrum" thesis`
14. `Maintrot Buser finite element Riemann surfaces genus 2 first eigenvalue maximum Teichmüller Bolza conjecture`
15. `Fortier Bourque Petri 2025 arXiv hyperbolic surfaces first eigenvalue new paper kissing number "first eigenvalue" bound`
16. `quasiplatonic OR "triangle surface" "critical point" first eigenvalue λ1 Teichmüller symmetric hyperbolic surface rigid "no invariant quadratic differential"`
17. `equilateral torus local maximum first eigenvalue flat tori Berger 1973 moduli space analogue hyperbolic surfaces genus 2`
18. `numerical optimization λ1 over moduli space genus 2 hyperbolic surfaces Hypermodes Fenchel-Nielsen maximum spectral gap search` (extended)

**arXiv API (export.arxiv.org, sorted by date, read 2021–2026):**
- `abs:Bolza` (all 91 hits listed; only the hyperbolic-geometry ones are relevant, and none is about λ₁ local max);
- `abs:"spectral gap" AND abs:hyperbolic AND abs:bootstrap`;
- `abs:"first eigenvalue" AND abs:hyperbolic AND abs:"genus 2"` and the same with `"genus two"`;
- `abs:eigenvalue AND abs:hyperbolic AND abs:"local maximum"`;
- `abs:eigenvalue AND abs:"hyperbolic surfaces" AND abs:"local maxima"`;
- `abs:Teichmuller AND abs:"first eigenvalue"`;
- `abs:"Laplace eigenvalue" AND abs:"hyperbolic surfaces"`;
- `abs:"first positive eigenvalue" AND abs:hyperbolic`;
- `abs:"Klein quartic"`;
- `abs:"maximal spectral gap"`;
- `abs:"spectral gap" AND abs:"genus 2"`;
- author listings: `au:Fortier_Bourque`, `au:Strohmaier`, `au:Karpukhin`, `au:Mazac`, `au:Kravchuk_P` (2022–2026);
- Petri's homepage publication list (fetched).

**Citation sweeps (Semantic Scholar API):** all papers citing:
- KMP 2111.12716;
- FBP-LP 2302.02540;
- SU13 1110.2150;
- Bonifacio 2111.13215.

The FBP-K citation list was rate-limited (429).

**Full texts downloaded and grepped** for `bolza|local max|locally max|critical point`:
- 2404.14479, 2507.10988, 2311.02632, 2307.06646, 2302.02540v3, 2111.12716v3, 2111.13215v2, 2107.09674,
  2111.14699v4;
- 2412.15111, 2403.00928, 2108.11825, 1604.02722, 2311.13330, 2608.01611, 2610.07193, 2601.13988, 2508.10465,
  2512.13139, 2312.03504v2, 2610.03023;
- 1405.4944, 1406.1080, 2312.11428, 1704.06384, math/0509398, 2105.03871.

The only hits relevant to λ₁ local maximality at B are the Str16 and Cook passages quoted in §2.

**Failed fetches:**
- archimede.mat.ulaval.ca (DNS);
- infoscience.epfl.ch (WAF / 429);
- arXiv 1309.5057 (HTML returned instead of the PDF; it is the earlier Nayatani–Shoda preprint and is superseded by
  1704.06384);
- the two CMS/GAPDE abstract PDFs were read but contain only generic talk abstracts.
