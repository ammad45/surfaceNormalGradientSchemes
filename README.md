# Surface Normal Gradient Schemes

**Implicit mesh-quality corrections for OpenFOAM** — surface-normal gradients (`snGrad`) and face interpolation (`divSchemes`) that fold non-orthogonality and skewness into the matrix coefficients instead of deferred explicit sources.

Developed for **transonic and supersonic flows**, as well as **injection moulding** and **turbomachinery**, where tet/poly meshes and strong gradients make classical deferred corrections noisy or unstable. In high-Mach regimes the acoustic (pressure–velocity) coupling is comparatively weak, so truncation errors from non-orthogonality and skewness in the discrete gradient and convective operators are less strongly damped by the continuum equations; absorbing those geometric corrections *implicitly* into the face coefficients therefore improves robustness and consistency of the discrete system.

| Build | Target | Library |
| --- | --- | --- |
| `./Allwmake` | ESI OpenFOAM (**v2412** / openfoam.com) | `libSnGrad.so` |
| `./Allwmake-tcae` | TCAE / Foundation (**OpenFOAM-dev**) | `libSnGradTCAE.so` |

Theory follows Hill, Deising, Acher, Klein, Bothe & Marschall,  
[*Boundedness-preserving implicit correction of mesh-induced errors…*](https://arxiv.org/abs/1705.11063), arXiv:1705.11063.

---

## Why implicit correction?

On a non-orthogonal, skewed mesh the face flux / Laplacian needs more than the simple neighbour difference. OpenFOAM usually adds that as an **explicit** deferred correction. Here the same terms are absorbed into an **effective coefficient**, so:

- the Laplacian / convection matrix stays consistent with the correction,
- no large deferred source appears every iteration,
- boundedness and diagonal dominance are easier to protect (limiters / fallbacks).

---

## Mathematical background

### Geometry

For a face with area vector $\mathbf{S}_f$, decompose (paper Eq. 7, over-relaxed form)

$$
\mathbf{S}_f = \boldsymbol{\Delta} + \mathbf{k},
\qquad
\boldsymbol{\Delta}
=
\frac{|\mathbf{S}_f|^2}{\mathbf{S}_f\cdot\mathbf{d}}\,\mathbf{d}.
$$

OpenFOAM’s `nonOrthDeltaCoeffs` and `nonOrthCorrectionVectors` correspond to

$$
\delta_{\perp}
=
\frac{|\boldsymbol{\Delta}|}{|\mathbf{S}_f|\,|\mathbf{d}|},
\qquad
\mathbf{k}_f
=
\frac{\mathbf{k}}{|\mathbf{S}_f|}.
$$

The **skew (non-conjunctionality) vector** is the offset of the face centre $C_f$ from the intersection $f'$ of the owner–neighbour line with the face plane:

$$
\mathbf{m} = C_f - f'.
$$

### Corrected surface-normal gradient

Paper Eq. 17 (normalised by $|\mathbf{S}_f|$):

$$
\begin{aligned}
(\nabla\phi)_f\cdot\hat{\mathbf{n}}
&=
\delta_{\perp}\,(\phi_N-\phi_P)
+
\mathbf{k}_f\cdot(\nabla\phi)_{f'}
+
\hat{\mathbf{n}}\cdot\bigl[\mathbf{H}_{f'}\cdot\mathbf{m}\bigr],
\end{aligned}
$$

where $\mathbf{H}=\nabla(\nabla\phi)$ is the Hessian (second derivatives) used for the skew term.

### Implicit form (diffusion)

Fold corrections into an effective delta (paper Eqs. 18–19):

$$
\delta_{\mathrm{eff}}
=
\delta_{\perp}
+
\frac{
\mathbf{k}_f\cdot(\nabla\phi)_{f'}
+
\hat{\mathbf{n}}\cdot(\mathbf{H}_{f'}\cdot\mathbf{m})
}{\phi_N-\phi_P},
$$

so that

$$
\mathrm{snGrad}(\phi)_f
=
\delta_{\mathrm{eff}}\,(\phi_N-\phi_P)
$$

enters `fvm::laplacian` with **no explicit correction source** (`correction()` returns `nullptr`).

### Implicit form (advection / interpolation)

For convection, Hill et al. Eq. 14 folds skew into the interpolation weight of a base scheme $w_{\mathrm{base}}$:

**Scalar**

$$
w_{\mathrm{eff}}
=
w_{\mathrm{base}}
+
\frac{\mathbf{m}\cdot(\nabla\phi)_{f'}}{\phi_N-\phi_P}.
$$

**Vector / tensor** (single scalar weight required by `gaussConvectionScheme`)

$$
w_{\mathrm{eff}}
=
w_{\mathrm{base}}
+
\frac{
\bigl(\mathbf{m}\cdot(\nabla\mathbf{U})_{f'}\bigr)\cdot(\mathbf{U}_N-\mathbf{U}_P)
}{|\mathbf{U}_N-\mathbf{U}_P|^2}.
$$

If $|\phi_N-\phi_P|$ (or $|\mathbf{U}_N-\mathbf{U}_P|$) is negligible, the correction is skipped on that face.

---

## Schemes

### Diffusion — `snGradSchemes`

| Keyword | Role | Explicit source? |
| --- | --- | --- |
| `implicitCorrected` | Full NO + NC folded into $\delta_{\mathrm{eff}}$; positivity / GAMG-safe limiter | No |
| `implicitStabCorrected` | Stabilised variant: adds $|\mathrm{corr}|$ so $\delta_{\mathrm{eff}}\ge\delta_{\perp}$ | No |
| `newCorrected` | Explicit NO + NC correction (classical deferred form) | Yes |
| `skewImpCorrected` | ESI-only experimental / WIP | — |

**Recommended**

- Vectors / rough meshes → **`implicitStabCorrected`**
- Scalars when you want the signed paper form → **`implicitCorrected`** (with the built-in face fallbacks)

### Advection — `divSchemes` / interpolation

| Keyword | Role |
| --- | --- |
| `implicitSkewCorrected <baseScheme> …` | Wraps a base interpolant (`linear`, `linearUpwind`, `linearUpwindV`, …) and applies Eq. 14 |

Passes through the base scheme’s `corrected()` / `correction()` (e.g. upwind gradient reconstruction still works).

---

## Build

### ESI OpenFOAM (openfoam.com)

```bash
source $WM_PROJECT_DIR/etc/bashrc   # e.g. OpenFOAM-v2412
cd surfaceNormalGradientSchemes
./Allwmake
# → $FOAM_USER_LIBBIN/libSnGrad.so
```

Clean:

```bash
./Allwclean
```

### TCAE / Foundation (openfoam.org)

```bash
source $WM_PROJECT_DIR/etc/bashrc   # OpenFOAM-dev / TCAE
cd surfaceNormalGradientSchemes
./Allwmake-tcae
# → $FOAM_USER_LIBBIN/libSnGradTCAE.so
```

Clean:

```bash
./Allwclean-tcae
```

> Source the matching OpenFOAM environment **before** building. The TCAE script exits if `$WM_PROJECT_DIR` is unset.

---

## Case setup

### 1. Load the library

`system/controlDict`:

```c++
libs
(
    "libSnGrad.so"        // ESI
    // "libSnGradTCAE.so" // TCAE
);
```

### 2. Surface-normal gradients

`system/fvSchemes`:

```c++
snGradSchemes
{
    default         orthogonal;
    snGrad(T)       implicitCorrected;
    snGrad(U)       implicitStabCorrected;
}

gradSchemes
{
    // leastSquares is preferred with these corrections (paper §5.1)
    default         leastSquares;
}
```

### 3. Skew-aware advection

```c++
divSchemes
{
    default         none;
    div(phi,U)      Gauss implicitSkewCorrected linearUpwindV grad(U);
    div(phi,k)      Gauss implicitSkewCorrected linearUpwind  grad(k);
    div(phi,omega)  Gauss implicitSkewCorrected linearUpwind  grad(omega);
}
```

---

## Repository layout

```text
surfaceNormalGradientSchemes/
├── Allwmake / Allwclean              # ESI → libSnGrad
├── Allwmake-tcae / Allwclean-tcae    # TCAE → libSnGradTCAE
├── REVIEW-TCAE.md                    # Port notes & numerics review
├── src/                              # ESI OpenFOAM-v2412
│   └── finiteVolume/
│       ├── finiteVolume/snGradSchemes/
│       └── interpolation/.../implicitSkewCorrected/
└── src-tcae/                         # Foundation / TCAE 26.03
    └── finiteVolume/
        ├── finiteVolume/snGradSchemes/
        └── interpolation/.../implicitSkewCorrected/
```

ESI `src/` also keeps `oldSchemes/` and `skewImpCorrectedSnGrad` for reference; they are not part of the TCAE library.

---

## Practical notes

- Prefer **`leastSquares`** (or a limited least-squares) for the gradients that feed the corrections.
- On highly skewed wall tets, **mesh quality** still dominates; these schemes reduce numerical noise, they do not replace prism layers.
- Coupled patches (processor / cyclic) are corrected in the current `implicitCorrected` / `implicitStabCorrected` / `implicitSkewCorrected` implementations.
- Near-uniform regions skip the $\Delta\phi$ denominator so the scheme does not blow up ($\div\mathrm{SMALL}$).

See [`REVIEW-TCAE.md`](REVIEW-TCAE.md) for the ESI↔TCAE API map and the bug-fix history.

---

## Reference

> R. Hill, D. Deising, H. Acher, H. Klein, D. Bothe, H. Marschall,  
> *Boundedness-preserving implicit correction of mesh-induced errors for VoF based heat and mass transfer*,  
> [arXiv:1705.11063](https://arxiv.org/abs/1705.11063).

Especially Eqs. **7**, **14**, **17–19**.

---

## License

Distributed under the **GNU General Public License v3** (or later), consistent with OpenFOAM.
