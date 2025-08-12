# p3c

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

