# p3c

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16836495.svg)](https://doi.org/10.5281/zenodo.16836495)

**p3c** is a Python library that computes a cubic polynomial transformation, based on the complex plane, for regularization purposes.  

Given a collocation point and the three interpolation points defining a quadratic boundary element, it returns:

- **Main pole**: `zeta0`, `eta0`  
- **Secondary pole** (only for quadratic elements): `zeta1`, `eta1`  
- **Transformed pole**: `zetatilde_0`, `etatilde_0`  
- **Iteration info**: `the_n_iteraciones` (number of Newton–Raphson iterations) and `the_xerr` (final error)  
- **Cubic coefficients**: `coef_a`, `coef_b`, `coef_c`, `coef_d`  
- **Callable functions**:
  - `p3c(xitilde)` – cubic polynomial: `a*xitilde**3 + b*xitilde**2 + c*xitilde + d`
  - `jac_p3c(xitilde)` – derivative of `p3c(xitilde)` with respect to `xitilde`

---

## Authors

- [Juan J. Granados](https://orcid.org/0000-0002-0707-2897)  
- [Rafael Gallego](https://orcid.org/0000-0002-7260-0940)  

---

## Questions, Bug Reports, and Feature Requests

If you encounter a problem, have a question, or would like to request a new feature, please use the  
[GitHub Issues page](https://github.com/juan-j-granados/p3c/issues).

We provide dedicated issue templates for:
- **Bug reports** – to help us reproduce and fix problems quickly.
- **Feature requests** – to suggest improvements or new capabilities.

Before opening an issue:
1. Check the [README](https://github.com/juan-j-granados/p3c#readme) and examples to see if your question is already answered.
2. Search existing issues to avoid duplicates.

---

## Background

In many boundary element method (BEM) formulations, certain integrals become nearly singular when the collocation point lies close to the element being integrated.  
The **p3c transformation** maps the integration domain onto a cubic polynomial in a transformed coordinate (`xitilde`) such that the singularity is regularized, improving numerical stability and accuracy.

This approach is particularly useful for:
- Quadratic boundary elements in 2D potential and elasticity problems
- Cases where the collocation point is close to the element
- Applications requiring high-precision evaluation of nearly singular integrals

The method computes the normalized pole in the complex plane and applies a Newton–Raphson–based refinement to determine the transformed pole, from which the cubic mapping and its Jacobian are derived.

---

## Installation

From the project root directory, you can either install normally or in editable mode (recommended during development), and then run the example script:

```bash
cd /path/to/project/root

# Standard installation
pip install .

# Editable installation (for development)
pip install -e .

# Run example
python -m examples.basic_usage

