# 02 - Degenerate-mode resonators: prior art and landscape

Scope: read-only web check (WebSearch/WebFetch, 2026-10). Only things I actually saw are quoted. Items marked [UNVERIFIED] come from my background knowledge, not from a page I read. Search snippets are weaker evidence than fetched pages; I say which is which.

## 1. How splits are modeled and trimmed in the MEMS/HRG literature

Core rule (n-mode pair splits only through the 2n-th Fourier harmonic of the imperfection) is well established in the literature found:

- Ma, Chen, Jin, Jin, Zheng, Jin, "Effects of Structural Dimension Variation on the Vibration of MEMS Ring-Based Gyroscopes", Micromachines 2021, doi:10.3390/mi12121483 (fetched, https://pmc.ncbi.nlm.nih.gov/articles/PMC8703920/). Quote: "the 4θ variations in the width of both the ring and the supporting beams cause modal coupling to occur and the degenerate 2θ mode pair to split in frequency." FEM tool used: COMSOL 5.6. Also reports analytic <10 Hz residual split after geometric compensation.
- Hosseini-Pishrobat and Tatar (Bilkent), "Analysis of quadrature and frequency split in a MEMS vibrating ring gyroscope with structural imperfections", Transducers 2021, doi:10.1109/Transducers50396.2021.9495619 (fetched abstract, https://repository.bilkent.edu.tr/items/bed90b20-d9ff-4bfb-b6ae-7bb8697018e7). Quote: "a 4° mode-shape rotation and 15Hz frequency split in a 60kHz VRG generate an approximate quadrature error of 1000°/s." Width non-uniformity rotates the mode shapes; suspension flaws split frequency. A longer journal version is at https://repository.bilkent.edu.tr/items/da9726e0-d3b3-4e1d-a405-9f9e53836f5e (search result only, not fetched).
- Hemispherical resonators: Wang, Pan, Qu et al., Sensors 2018, 18(11), 3772, doi:10.3390/s18113772 (fetched, https://pmc.ncbi.nlm.nih.gov/articles/PMC6263873/). Quote: "This paper mainly focuses on reducing M4, which is the cause of the frequency split of the resonator and the main sources of bias and random errors in HRGs." Search snippets (not verified on the page) say the first three mass-defect harmonics matter an order of magnitude less than the 4th, and that trimming algorithms use Fourier expansion and harmonic superposition. Trimming methods named in results: mechanical, laser ablation, chemical etching, ion-beam etching.
- McWilliam and co-workers (Nottingham) gave an analytic frequency split of ring resonators as a Fourier representation of anisotropic elastic properties (search snippet only; I did not open the paper).
- Electrostatic tuning: results state that DC voltage on electrodes lowers the frequency of the corresponding mode via negative electrostatic stiffness, which achieves mode matching (snippet from https://www.mdpi.com/2076-3417/11/3/1129/htm and PMC9319629; not fetched).
- Mass-stiffness decoupled trimming for ring gyros: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9319629/ and PMC10609602 (search results only).

Not retrieved: primary Shkel-lab, Ayazi-lab and Lynch papers. I did not find them in search results, so I cannot cite them. Lynch's "Vibratory gyro analysis by the method of averaging" and Rourke/McWilliam/Burdess ring analyses [UNVERIFIED titles from memory]; verify before citing. Anyone writing this up should pull Ayazi's disk-resonator-gyro and Shkel's ring/wineglass trimming papers directly.

## 2. Open tools for imperfection-to-split maps

Found:
- Bindel (Cornell/Berkeley), HiQLab, a FE tool for resonant MEMS with parameter-dependent eigensolvers, thesis "Structured and Parameter-Dependent Eigensolvers for Simulation-Based Design of Resonant MEMS", UCB EECS-2006-108 (https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-108.html). Focus is damping/anchor loss and continuation, not a degenerate-cluster split matrix as far as the abstract shows.
- arXiv:1008.0944, "Analysis of the manufacturing variation in a coupled micro resonators array based on its designed value and measured eigenfrequencies": perturbation theory for fabricated-system eigenvalues (snippet only; closest academic neighbor, read it before claiming novelty).
- Kratos Multiphysics reportedly has eigenvalue sensitivities (snippet from a search; nonlinear structural context, not degenerate clusters).
- Commercial: Abaqus imperfection keywords (buckling-mode superposition) are for buckling, not splitting.

No open tool found that returns a first-order splitting matrix plus degeneracy-preserving kernel for arbitrary geometry. Search coverage was shallow (3 queries on this point); absence of evidence only. The practice visible in the papers above is ad hoc COMSOL (Ma et al.) plus analytic Fourier models. Generic "eigenvalue derivative with repeated eigenvalues" is classical (Kato, Lancaster, Nelson-type results) [UNVERIFIED from memory]; the tool's contribution would be packaging it behind FE adapters plus the kernel and second-order pieces.

## 3. Other domains

- Photonic crystals: search results state degenerate cavity modes split under perturbation and that symmetry analysis predicts it (https://arxiv.org/pdf/2202.07620 and others; snippets). Quadratic degeneracies are "robust against perturbations that preserve crystal symmetry" but split into linear degeneracies when 90-degree rotational symmetry is broken (https://ar5iv.arxiv.org/html/0803.1854, snippet). This is the same A-matrix/kernel structure, but Maxwell eigenproblems are generalized Hermitian with different FE (edge elements, curl-curl), so adapters differ; Bloch k-dependence adds parameters (k.p is the first-order split matrix in k).
- Dirac points / graphene: k.p method yields Dirac equation from a two-band degeneracy (arXiv:1105.1351, snippet). The degeneracy-preserving kernel = symmetry-allowed perturbations.
- Acoustic metamaterials, civil structures (symmetric towers, tanks, domes with repeated frequencies), quantum dots/molecular Jahn-Teller: I did not find or verify specific references. Plausible but unsupported by this search.

## 4. Proposed one-day example: ring 2D plane stress, n=2 pair

Setup in scikit-fem: thin ring, mean radius R, width w << R, plane stress, free-free (or lightly constrained to kill rigid modes), mesh ~ 2-3 quadratic elements across width. Imperfection parameterized as radial width/edge perturbation w(theta)=w0(1+sum_k a_k cos(k theta)+b_k sin(k theta)), k=1..6, or the simpler mass-density field rho(theta)=rho0(1+sum ...) (no remeshing, K and M derivatives are cheap elementwise). Start with density/stiffness fields to avoid mesh-morphing derivatives.

Known result (literature above): only the 4theta harmonic splits the 2theta (n=2) pair at first order; 1,2,3,5,6 do not. For a mass perturbation delta_rho = rho0 * eps * cos(4(theta - phi)), my own first-order derivation (not from a source; verify numerically) gives A = -(omega^2 eps/4) [[cos4phi, sin4phi],[sin4phi, -cos4phi]] in the cos/sin 2theta basis, i.e. eig(A) = +-omega^2 eps/4, so Delta f = f * eps/4 (and mean shift comes from the 0th harmonic). Stiffness harmonics follow the same 2n rule with different coefficients. Analytic checks:
  1. Singular values of A vs harmonic index: nonzero only at k=4 (to ~1e-10 relative in the discretized ring, up to mesh symmetry effects: use a mesh with exact C4 or better Cn symmetry).
  2. Delta f = f*eps/4 for pure 4theta mass imperfection; direction phi sets the mode-rotation angle (phi mod 90 degrees).
  3. Kernel of the splitting map on the (a_k,b_k) parameter space contains all k != 0 mod 4 harmonics (and k=0 trivially); dimension check against the count.
  4. Compare to finite-difference eigenvalues of the perturbed ring at eps~1e-3; second derivative of the cluster sum vs central differences.
  5. Optional: trim demo: choose a_4,b_4 corrections in the kernel-complement to null the split (the "trimming" story), citing Wang 2018 for the HRG analog.
Pitfall: k=8 splits the n=4 pair, not n=2; and the n=2 pair also couples at second order via other harmonics (e.g. 2+2), which the exact second derivative should expose: a nice added panel beyond the standard first-order rule.

## 5. Recommendation

Value: moderate. The first-order rule is textbook-level in this field, so the example validates the tool but is not itself new. Real utility is for non-ring geometries (disk resonator gyros with spokes, wineglass shells, anisotropic silicon) where closed forms fail and people use COMSOL parameter sweeps; a kernel (trim directions that preserve degeneracy, or the complement that cancels the split) and exact second-order cluster-sum are the differentiators. Novelty: mostly packaging plus the kernel/second-order/Ky Fan analysis; I did not find a competing open FE-plus-sensitivity tool, but search was shallow, and Bindel's HiQLab and arXiv:1008.0944 are closest prior art. Suggest: build the ring example, publish as validation, and only then decide on pursuing the disk-gyro design-for-trim application. Priority: worth a day; not a standalone research direction without a design-optimization result.
