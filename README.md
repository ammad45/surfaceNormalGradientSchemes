# Surface Normal Gradient Schemes

**Implicit mesh-quality corrections for OpenFOAM** — surface-normal gradients (`snGrad`) and face interpolation (`divSchemes`) that fold non-orthogonality and skewness into the matrix coefficients instead of deferred explicit sources.

Developed for **transonic and supersonic flows**, as well as **injection moulding** and **turbomachinery**, where tet/poly meshes and strong gradients make classical deferred corrections noisy or unstable. In high-Mach regimes the acoustic (pressure–velocity) coupling is comparatively weak, so truncation errors from non-orthogonality and skewness in the discrete gradient and convective operators are less strongly damped by the continuum equations; absorbing those geometric corrections *implicitly* into the face coefficients therefore improves robustness and consistency of the discrete system.

| Build | Target | Library |
| --- | --- | --- |
| `./Allwmake` | ESI OpenFOAM (**v2412** / openfoam.com) | `libSnGrad.so` |
| `./Allwmake-tcae` | TCAE / Foundation (**OpenFOAM-dev**) | `libSnGradTCAE.so` |

Based on Hill et al., [*Boundedness-preserving implicit correction of mesh-induced errors…*](https://arxiv.org/abs/1705.11063), arXiv:1705.11063.

---

## Idea

On a non-orthogonal or skewed mesh, OpenFOAM usually applies geometric corrections as **explicit** deferred terms. These schemes fold the same corrections into an **effective face coefficient**, so:

- the Laplacian / convection matrix stays consistent with the correction,
- large deferred sources are avoided each iteration,
- boundedness and diagonal dominance are easier to protect.

---

## Schemes

### Diffusion — `snGradSchemes`

| Keyword | Role | Explicit source? |
| --- | --- | --- |
| `implicitCorrected` | Non-orthogonal + skew correction folded into the face delta; positivity / GAMG-safe limiter | No |
| `implicitStabCorrected` | Stabilised variant that keeps the effective delta at least as large as the orthogonal one | No |
| `newCorrected` | Explicit non-orthogonal + skew correction (classical deferred form) | Yes |
| `skewImpCorrected` | ESI-only experimental / WIP | — |

**Recommended:** vectors / rough meshes → `implicitStabCorrected`; scalars → `implicitCorrected`.

### Advection — `divSchemes`

| Keyword | Role |
| --- | --- |
| `implicitSkewCorrected <baseScheme> …` | Wraps a base interpolant (`linear`, `linearUpwind`, `linearUpwindV`, …) and adds an implicit skew correction to its weights |

The base scheme’s own explicit correction (e.g. upwind gradient reconstruction) is passed through unchanged.

---

## Build

### ESI OpenFOAM (openfoam.com)

```bash
source $WM_PROJECT_DIR/etc/bashrc   # e.g. OpenFOAM-v2412
cd surfaceNormalGradientSchemes
./Allwmake
# → $FOAM_USER_LIBBIN/libSnGrad.so
```

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

```bash
./Allwclean-tcae
```

Source the matching OpenFOAM environment **before** building. The TCAE script exits if `$WM_PROJECT_DIR` is unset.

---

## Case setup

### Load the library

`system/controlDict`:

```c++
libs
(
    "libSnGrad.so"        // ESI
    // "libSnGradTCAE.so" // TCAE
);
```

### Surface-normal gradients

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
    default         leastSquares;
}
```

### Skew-aware advection

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
└── src-tcae/                         # Foundation / TCAE 26.03
```

---

## Notes

- Prefer `leastSquares` (or limited least-squares) for gradients used by these corrections.
- Mesh quality still matters; these schemes reduce numerical noise but do not replace good near-wall layers.
- Coupled patches (processor / cyclic) are corrected.
- Near-uniform faces skip the correction so the scheme does not blow up.

See [`REVIEW-TCAE.md`](REVIEW-TCAE.md) for the ESI↔TCAE port details.

---

## Reference

> R. Hill, D. Deising, H. Acher, H. Klein, D. Bothe, H. Marschall,  
> *Boundedness-preserving implicit correction of mesh-induced errors for VoF based heat and mass transfer*,  
> [arXiv:1705.11063](https://arxiv.org/abs/1705.11063).

---

## License

Distributed under the **GNU General Public License v3** (or later), consistent with OpenFOAM.
