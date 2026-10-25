# TCAE 26.03 port – implicit mesh-skewness correction schemes

Port target: **TCAE 26.03** (OpenFOAM-dev from the OpenFOAM Foundation,
`openfoam.org`).
Reference paper: Hill, Deising, Acher, Klein, Bothe, Marschall,
*Boundedness-preserving implicit correction of mesh-induced errors for
VoF based heat and mass transfer*, arXiv:1705.11063.

The original code was written against **ESI OpenFOAM-v2412**
(`openfoam.com`).  These two forks have diverged enough that the code
does **not** compile against TCAE 26.03 without changes.  This document
records the port, reviews the numerics, and adds a matching **implicit
advection** scheme for turbomachinery (§3.3).

Ported source lives beside the original:
```
sngrads/
├── src/            (untouched ESI v2412 code)
└── src-tcae/       (TCAE 26.03 port + new advection scheme)
    └── finiteVolume/
        ├── finiteVolume/snGradSchemes/                 (diffusion)
        │   ├── implicitCorrectedSnGrad/
        │   ├── implicitStabCorrectedSnGrad/
        │   └── newCorrectedSnGrad/
        └── interpolation/surfaceInterpolation/schemes/ (advection)
            └── implicitSkewCorrected/
```

Build with `./Allwmake-tcae` after sourcing the TCAE bashrc.  The
resulting library is `$FOAM_USER_LIBBIN/libSnGradTCAE.so`; add it via
`libs ("libSnGradTCAE.so");` in `system/controlDict`.

---

## 1. Port: ESI OpenFOAM-v2412 → OpenFOAM-org dev (TCAE 26.03)

### 1.1 API changes applied

| ESI (openfoam.com, v2412) | Foundation (openfoam.org, dev) |
|---|---|
| `GeometricField<Type, fvPatchField, volMesh>` | `VolField<Type>` |
| `GeometricField<Type, fvsPatchField, surfaceMesh>` | `SurfaceField<Type>` |
| `mesh.gradScheme("grad(x)")` | `mesh.schemes().grad("grad(x)")` |
| `virtual bool corrected() const noexcept` | `virtual bool corrected() const` |
| `ssf.setOriented()` | *(no equivalent; removed)* |
| `IOobject::AUTO_WRITE, false` for delta field | Use `surfaceScalarField::New(...)` |
| `snGrad<T>::correction` returns `nullptr` default | Same, but base class returns tmp<>(nullptr) |

### 1.2 Files ported

The following files were rewritten for TCAE. Semantics fixed where
noted (details in §2).

| File | Ported | Fixes vs. ESI original |
|---|---|---|
| `implicitCorrectedSnGrad.[HC]` | yes | 3 bug fixes (see §2.1) |
| `implicitStabCorrectedSnGrad.[HC]` | yes | 1 fix (coupled patches) |
| `newCorrectedSnGrad.[HC]` | yes | Cleanup / API port |
| `implicitSkewCorrected.[HC]` (**advection**) | **new** in the TCAE port | Implicit non-conjunctionality correction for `div(phi,U)` (paper Eq. 14 with vector reduction; see §3.3) |
| `skewImpCorrectedSnGrad.*` | **not ported** | Body was `NotImplemented`. Skipping until you decide what it should do. |
| `oldSchemes/*` | **not ported** | Older WIP; not needed. |

The library was compiled cleanly against TCAE 26.03:
```
/p/home/jusers/ammad1/jureca/OpenFOAM/ammad1-dev-26.03-cfdsupport/
     platforms/linux64Gcc132DPInt32Opt/lib/libSnGradTCAE.so   ✔
```

All five OpenFOAM primitive types are registered (scalar, vector,
sphericalTensor, symmTensor, tensor), confirmed via `nm | c++filt`.

---

## 2. Correctness review vs. arXiv:1705.11063

### 2.1 Derivation (used by all three schemes)

OpenFOAM decomposes the face area vector as
```
Sf = Δ + k        (Eq. 7 in the paper)
```
with the **over-relaxed** choice
```
Δ = |Sf|²/(Sf·d) · d,   |Δ|/|d| = |Sf|/(|d| cos θ) = |Sf|·nonOrthDeltaCoeffs.
```
The face-normal gradient with both non-orthogonal (NO) and
non-conjunctional (NC / skew) corrections is (paper Eq. 17, divided by
|Sf|):
```
snGrad_f = |Δ|/(|Sf||d|) · (φ_N − φ_P)                    ← orthogonal
         + (k/|Sf|) · (∇φ)_{f'}                            ← NO
         + (Sf/|Sf|) · [ ∇(∇φ)_{f'} · m ]                  ← NC (skew)
        = nonOrthDeltaCoeffs · (φ_N − φ_P)
         + nonOrthCorrectionVectors · (∇φ)_f
         + n̂ · ( H_f · m )
```
This matches the definitions in `surfaceInterpolation.C` (`unitArea −
delta·nonOrthDeltaCoeffs` = `k/|Sf|`) and `skewCorrectionVectors.C`
(m = C_f − f′ where f′ is the intersection of the P–N line with the
face). So the arithmetic in your code correctly maps to Eq. 17. ✔

The **implicit** treatment (paper Eq. 18/19) then folds the two
correction terms into a modified delta:
```
deltaEff = nonOrthDeltaCoeffs
         + [ nonOrthCorrectionVectors · (∇φ)_f + n̂ · (H_f · m) ] / (φ_N − φ_P)
```
`snGrad = deltaEff · (φ_N − φ_P)` then plugs directly into any
implicit Laplacian – no source term is needed.  That is exactly the
idea used in `implicitCorrectedSnGrad::deltaCoeffs()`.

### 2.2 Bugs found in the ESI original (fixed in the port)

**Bug 1 – divide-by-SMALL blow-up** (`implicitCorrectedSnGrad.C`):
```cpp
if (mag(delta) > SMALL) tNonOrth /= delta;
else                    tNonOrth /= SMALL;   // ← blows up
deltaCoeffs[facei] += tNonOrth()[facei];
```
Whenever two neighbouring cells have nearly identical field values
(happens *everywhere* at the start of a run and in almost-uniform
regions), the correction term is divided by `SMALL ≈ 1e-15`, producing
values on the order of `1/SMALL` that get added into the diagonal.
The Laplacian matrix is then no longer physical and either diverges
immediately or converges to garbage.

Fixed in `src-tcae/.../implicitCorrectedSnGrad.C`:
```cpp
if (mag(dPhi) > SMALL)
{
    const scalar dPhiSafe = sign(dPhi)*max(mag(dPhi), SMALL);
    deltaCoeffs[facei] += corr[facei]/dPhiSafe;
}
// else: leave deltaCoeffs alone.  If (φ_N − φ_P) is really zero
// the correction is meaningless anyway; the orthogonal term is
// already correct.
```

**Bug 2 – coupled patches not corrected**
(`implicitCorrectedSnGrad.C`, `implicitStabCorrectedSnGrad.C`):
Both files only walk the internal-face `forAll(owner, facei)` loop and
never touch `deltaCoeffs.boundaryField()`.  On processor and cyclic
patches the effective delta is therefore the uncorrected orthogonal
value while adjacent internal faces are corrected – i.e. **parallel
runs and periodic cases will give a different answer to a serial /
non-periodic run of the same case.**

Fixed in both files by looping over `deltaBf`, filtering with
`pDelta.coupled()`, and using `patchInternalField() /
patchNeighbourField()` to get φ_P / φ_N across the coupled patch. On
non-coupled patches `nonOrthCorrectionVectors` is zero by construction
in `surfaceInterpolation.C::makeNonOrthCorrectionVectors()`, so
skipping them is correct.

**Bug 3 – laplacian diagonal / GAMG stability**

Unbounded signed corrections (paper Eq. 19 without positivity) can drive
`deltaCoeffs` to (near) zero on skewed industrial meshes.  The pressure
Poisson laplacian then gets a zero cell diagonal and GAMG crashes in
`scale()` with `FOAM_SIGFPE` (divide by `D[i]`).

Fixed in `src-tcae/.../implicitCorrectedSnGrad.C` without reverting to
the ESI global `max(deltaCoeffs, 0)` limiter (which the paper §5.2
reports hurts VoF boundedness):

1. Each component increment uses `max(dC, -delta)` so `deltaCoeffs`
   never crosses zero while accumulating.
2. After all components, faces with `deltaCoeffs < 0.01*deltaOrth` fall
   back to the orthogonal `nonOrthDeltaCoeffs` on that face only.

This keeps the implicit correction wherever it is well-behaved and
reverts locally on pathological faces.  For a fully positive
construction use `implicitStabCorrected` (§2.5).

### 2.3 Non-orthogonality: correctness

`mesh.nonOrthCorrectionVectors() & tGradF()` is dimensionally
`(dimless) & (dimQuantity/dimLength) = dimQuantity/dimLength`, which
matches `nonOrthDeltaCoeffs · (φ_N − φ_P)`. ✔

For the *linear* interpolation of the cell gradient onto the face,
the paper uses the same convention (Eq. 16 / Eq. 17: `(∇φ)_{f'} =
f(∇φ_C, ∇φ_D)` interpolated). The linear interpolation used in the
code is standard second-order interpolation. ✔

The paper (§5.1) reports that **least-squares gradient outperforms
Gauss** for the correction. This is *not* hard-coded – the user picks
the gradient by setting `gradSchemes { grad(<name>) leastSquares; }`
in `fvSchemes`, and the code correctly picks that scheme up via
`mesh.schemes().grad("grad(" + vf.name() + ')')`. Add the following to
your case's `fvSchemes` to activate the paper's recommendation:
```
gradSchemes
{
    default              leastSquares;
    grad(alpha.water)    leastSquares;
}
```

### 2.4 Non-conjunctionality (skew): correctness

The construction
```cpp
tmp<VolField<HessianType>> tHvol( fvc::grad(tGradVol()) );
tmp<SurfaceField<HessianType>> tHf( linear<HessianType>(mesh).interpolate(tHvol()) );
corr += tNhat() & (tHf() & tSkewVec());
```
implements `n̂ · (H_f · m)` in Eq. 17.

Two subtleties:

1. **`scv()` has dimless dimensions but stores metres.** OpenFOAM's
   `skewCorrectionVectors` field carries `dimless` dimensions (see
   `skewCorrectionVectors.C` L60) even though its numerical values are
   in metres.  The multiplication by `unitLength = 1 [m]` restores
   proper length dimensions so that `H · m` has units of gradient
   ([φ/m]), matching the other terms.  ✔ (This is what the ESI code
   did too.)

2. **Nested Gauss gradients amplify noise.**  The Hessian is computed
   as `fvc::grad(tGradVol())`.  On strongly distorted meshes this
   can be noisy, particularly at boundaries.  The paper (§5.2) reports
   that **NC correction gives no measurable improvement** even for
   scalar fields with strong gradient curvature (interface jumps),
   which is why the standard OpenFOAM `corrected` scheme omits it.  So
   your NC correction is **correct** but the paper suggests you can
   often turn it off.  In the port you can turn it off by using a
   perfectly orthogonal mesh (then `scv.skew()` returns false and the
   Hessian is not computed) or by adding a runtime switch (see §3
   below).

### 2.5 The stabilised variant

`implicitStabCorrectedSnGrad.C` adds `|corr|` – the absolute value –
to `nonOrthDeltaCoeffs`.  This is **not** the paper's method: the
paper's Eq. 19 uses the signed correction divided by (φ_N − φ_P).
What the stabilised scheme does instead is guarantee
`deltaEff ≥ nonOrthDeltaCoeffs`, so the Laplacian is always
diagonally dominant and cannot break the solver on ugly meshes.

Trade-off:
- Pro: never divides by (φ_N − φ_P); always positive; robust.
- Con: it is a genuinely different discretisation, not a
  convergent-in-space correction. The order of accuracy on
  systematically distorted meshes will be lower than the
  `implicitCorrected` scheme.

Suggestion: rename the scheme to something that makes this clear (e.g.
`nonOrthMagStabilised`) so nobody confuses it with the paper's SISC.

### 2.6 Comparison with the paper's schemes

| Paper name             | Term        | This code name                | Kind                          | What it does                                                                       |
|---                     |---          |---                            |---                            |---                                                                                 |
| explicit EC            | diffusion   | `newCorrected`                | `snGradScheme`                | Deferred correction added to snGrad (paper Eq. 3 / 17)                             |
| NO / NO/NC (implicit)  | diffusion   | `implicitCorrected`           | `snGradScheme`                | Signed correction folded into deltaCoeff (paper Eq. 18–19)                         |
| n/a (heuristic)        | diffusion   | `implicitStabCorrected`       | `snGradScheme`                | `|correction|` folded into deltaCoeff – always positive; robust but not convergent |
| SISC (implicit, unlim.)| advection   | `implicitSkewCorrected`       | `surfaceInterpolationScheme`  | `m · ∇φ / (φ_N − φ_P)` folded into interpolation weight (paper Eq. 14, no TVD)     |

The paper's "NO only" (no NC) variant of the diffusion correction is
not exposed here – it is easy to add (see §3.1).  The paper's TVD
limiter on the advection weight (Eq. 15) is intentionally omitted
from `implicitSkewCorrected` – see §3.3 for why.

---

## 3. Recommendations / follow-ups

### 3.1 Independent switch for NC correction

Right now `useSkew` is decided by
`skewCorrectionVectors::New(mesh).skew()`.  That covers "does the mesh
need NC correction at all?" but there is no way to say *"I have a
skewed mesh but I only want NO correction"*.  The paper found (§5.2,
Fig. 14–17) that NC has essentially no measurable effect on
convergence, so being able to turn it off cheaply is useful for
benchmarking and to avoid the cost of computing the Hessian.

Suggested API (both `.H` files):
```cpp
implicitCorrectedSnGrad(const fvMesh& mesh, Istream& is)
:
    snGradScheme<Type>(mesh),
    // optional integer flag: 0 = NO only, 1 = NO + NC (default = 1)
    nonConjCorr_(is.eof() ? 1 : readLabel(is))
{}
```
Then in `deltaCoeffs()`:
```cpp
if (nonConjCorr_ && scv.skew()) { ... Hessian block ... }
```

### 3.2 Add the paper's "NO only" scheme

For a clean comparison with the paper you probably want three
schemes: `implicitCorrected_NO`, `implicitCorrected_NC`,
`implicitCorrected_NOandNC`.  Trivial once the switch above is in.

### 3.3 Advection: `implicitSkewCorrected`

An **implicit** skew (non-conjunctionality) correction for the
advection term was added in
`src-tcae/finiteVolume/interpolation/surfaceInterpolation/schemes/implicitSkewCorrected/`.
It plugs into `divSchemes`, not `snGradSchemes` – the correction acts
on the face **interpolation weight**, not on the surface-normal
gradient.

Usage in `system/fvSchemes`:
```
divSchemes
{
    default              none;
    div(phi,U)           Gauss implicitSkewCorrected linearUpwindV grad(U);
    div(phi,k)           Gauss implicitSkewCorrected linearUpwind  grad(k);
    div(phi,omega)       Gauss implicitSkewCorrected linearUpwind  grad(omega);
}
gradSchemes
{
    default              leastSquares;
    grad(U)              leastSquares;
}
```

**Derivation.** The paper (Eq. 14, unlimited variant) rewrites the
non-conjunctional correction `φ_f += m · ∇φ_{f'}` into a modified
scalar interpolation weight

```
w_eff = w_base + m · (∇φ)_{f'} / (φ_N − φ_P)   ← paper Eq. 14 (scalar)
```

so that `φ_f = φ_P + w_eff · (φ_N − φ_P)` fits into OpenFOAM's
implicit convection matrix without an explicit deferred-correction
term.  This is *the* reason the paper calls it "implicit".

**Vector generalisation (needed for U).** OpenFOAM's
`gaussConvectionScheme::fvmDiv` only consumes **one scalar weight per
face** for all components of a `VolField<Type>`:

```cpp
tmp<surfaceScalarField> tweights = tinterpScheme_().weights(vf);
fvm.lower() = -weights*faceFlux;
fvm.upper() = fvm.lower() + faceFlux;
```

Applying Eq. 14 component-wise gives three different weights per face
for a vector, which cannot go into `fvm::div` implicitly.  The
implemented reduction is the **least-squares scalar** that best
represents the component-wise correction, obtained by dotting Eq. 14
with `(U_N − U_P)`:

```
w_eff = w_base
      + [ (m · ∇U)_{f'} · (U_N − U_P) ] / |U_N − U_P|²
```

For a scalar field this reduces exactly to the paper's Eq. 14
(unlimited).  For a vector, the correction component orthogonal to
`(U_N − U_P)` is dropped – this is the price of keeping the scheme
fully implicit.

**Design notes.**
- Wraps any base surface interpolation scheme (linear, upwind,
  linearUpwind, linearUpwindV, ...).
- Passes through the base scheme's `corrected()` / `correction()`, so
  e.g. `linearUpwind`'s explicit second-order upwind reconstruction is
  preserved on top of the implicit skew correction.
- `useSkew` guard: if `skewCorrectionVectors::New(mesh).skew()` is
  false, `weights()` returns the base weights unchanged (zero cost).
- Sign-preserving floor: any face where `|vf_N − vf_P|` (in the
  Type-aware sense) is below `SMALL` is skipped, so nearly uniform
  fields do not destabilise the diagonal.
- Coupled patches (processor / cyclic) get the same correction using
  `patchInternalField()` / `patchNeighbourField()`.

**Not implemented on purpose:**
- TVD limiter on `w_eff` (paper Eq. 15). Boundedness is meaningful
  for α ∈ [0, 1] but *not* for velocity or turbulence quantities in a
  turbomachinery context – applying TVD component-wise to a vector
  breaks Galilean invariance.
- Virtual upwind reconstruction (paper Eq. 12–13). Needed only for
  bounded VoF advection.
- Non-orthogonality correction of the interpolation. Advection depends
  only on the face-normal *flux* `phi = Sf · u`, so non-orthogonality
  does not enter face interpolation – it only enters the diffusion
  surface-normal gradient (already handled by `implicitCorrected`).

### 3.4 Test suggestions

The paper is very specific about how their test cases are set up:

1. **Systematically-distorted 2D mesh** (their Fig. 6). Random
   distortion causes error cancellation and hides skewness errors –
   don't use it as your primary test.
2. **Species diffusion** (Fig. 13): 0.04 m tall 2D strip, D_g = 1e-1,
   D_l = 1e-5, H ∈ {0.033, 0.2, 1, 5, 30}, run to t = 0.5 s, evaluate
   over the vertical line L3.  Compare against a 1D FDM reference on
   the same domain (their "Exact" line).
3. **Convergence** (Fig. 16–17): plot the average molar concentration
   in the gas phase at t = 0.1 s and t = 0.5 s against Δx.  Uncorrected
   should give roughly zero-order; NO and NO/NC should recover
   first-order (H ≪ 1 or H ≫ 1) or second-order (H = 1).

The paper's §5.1 advection test (VoF plug flow, circle translation)
targets scalar α with strict boundedness [0, 1] and TVD limiting.
Because `implicitSkewCorrected` intentionally drops the TVD limiter
(inappropriate for U / turbomachinery), it should **not** be used for
VoF-style α advection – use one of OpenFOAM's dedicated
compressive VoF schemes there.  Reasonable tests for the
turbomachinery use case are:

* Poiseuille / Couette on a systematically skewed 2D channel:
  verify that the analytic linear-in-y velocity profile is preserved
  when `divSchemes{ div(phi,U) Gauss implicitSkewCorrected linearUpwindV grad(U); }`
  is used, and compare against `Gauss linearUpwindV grad(U)` and
  `Gauss skewCorrected linearUpwindV grad(U)`.
* Uniform-flow benchmark on a distorted mesh: initialise
  U = (1, 0, 0), fix at inlet, zero-gradient at outlet, wall on top
  and bottom.  The exact solution is U = (1, 0, 0) everywhere;
  measure `max|U − (1, 0, 0)|` for the three schemes.  Only the
  skew-corrected variants should reach machine zero.
* Turbomachinery cascade with a periodic mesh: turn on periodic
  patches, run with each scheme, check that residuals converge and
  that the pressure/velocity fields near blade leading edges are
  smooth (no zig-zag in the wall-normal direction – a classic
  skewness artefact).

### 3.4 Minor cleanups (already applied in the port)

* Removed obsolete commented-out code (`fullGradCorrection`,
  duplicated `computeVolumeGradient`, ESI-specific `setOriented`,
  duplicated comment banners).
* Kept templated `<Type>` uniform across the three files.
* Header/source pairing consistent with base `snGradScheme.H`
  (`~implicitCorrectedSnGrad() {}` – base has `= 0` for `type()` only).

---

## 4. Answering your specific concerns

> "Normally in OpenFOAM the non-orthogonal part is treated explicitly
> and here I want to treat it implicitly"

Yes – that is exactly what `implicitCorrectedSnGrad::deltaCoeffs()`
does. Because the scheme's `corrected()` returns **false**, the base
`snGradScheme::snGrad(vf)` does *not* add an explicit correction; the
NO / NC contributions are entirely inside the delta coefficient, and
thus inside the implicit Laplacian matrix.  This is the paper's Eq.
18–19, verbatim.

> "and there is no non-conjunctional correction (skew correction)
> which I am adding in this scheme"

The NC / skew term `n̂ · (H_f · m)` is added conditionally on
`skewCorrectionVectors::New(mesh).skew()`. The Hessian is computed
via nested `fvc::grad`. The math is correct (paper Eq. 17). The paper
itself concludes NC adds little value for the diffusion case, but
**your code does compute it correctly** and you can turn it off
easily by running on a mesh that OpenFOAM classifies as non-skew.

> "port it to TCAE 26.03"

Done. `libSnGradTCAE.so` builds cleanly, no warnings. Both implicit
schemes and the explicit `newCorrected` are registered for all five
primitive types.
