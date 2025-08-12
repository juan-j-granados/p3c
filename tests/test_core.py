# To run tests/test_core.py and check that everything is working:
# 1. Install pytest (if not already installed):
#     pip install pytest
# 2. Navigate to the project root (the directory containing pyproject.toml).
# 3. Run pytest:
#     pytest -q
# The -q flag runs tests in “quiet” mode, showing only minimal output.

import math
from p3c.core import build_p3c

def test_build_and_callable():
    X1, Y1 =   1.,   1.
    X2, Y2 =   3.,   2.5
    X3, Y3 =   5.,   2. 

    Xc, Yc = 4.41, 2.13

    zeta0, eta0, zeta1, eta1, \
        zetatilde_0, etatilde_0, the_n_iteraciones, the_xerr,  \
        coef_a, coef_b, coef_c, coef_d, \
        p3c_of_xitilde, jac_p3c_of_xitilde  = build_p3c(Xc,Yc, X1,Y1, X2,Y2, X3,Y3)


    p3c_of_zeta0 = p3c_of_xitilde(zeta0)
    jac_p3c_of_zeta0 = jac_p3c_of_xitilde(zeta0)

    p3c_of_zetatilde_0 = p3c_of_xitilde(zetatilde_0)
    jac_p3c_of_zetatilde_0 = jac_p3c_of_xitilde(zetatilde_0)


    zeta0_d, eta0_d =  0.7524705856764565, -0.09447435736471754
    zeta1_d, eta1_d = -0.2524705856764565, -1.9055256426352825

    zetatilde_0_d, etatilde_0_d =  0.4375586278465583, -0.47346508373606316
    the_n_iteraciones_d, the_xerr_d = 9, 0.0

    p3c_of_zeta0_d, jac_p3c_of_zeta0_d = 0.860625147740979, 0.43171690846965405

    p3c_of_zetatilde_0_d, jac_p3c_of_zetatilde_0_d = zeta0_d, 0.29930725815902925


    assert math.isclose(zeta0, zeta0_d, rel_tol=1e-9, abs_tol=0.0)
    assert math.isclose(eta0, eta0_d, rel_tol=1e-9, abs_tol=0.0)
    assert math.isclose(zeta1, zeta1_d, rel_tol=1e-9, abs_tol=0.0)
    assert math.isclose(eta1, eta1_d, rel_tol=1e-9, abs_tol=0.0)
    assert math.isclose(zetatilde_0, zetatilde_0_d, rel_tol=1e-9, abs_tol=0.0)
    assert math.isclose(etatilde_0, etatilde_0_d, rel_tol=1e-9, abs_tol=0.0)
    assert the_n_iteraciones == the_n_iteraciones_d
    assert math.isclose(the_xerr, the_xerr_d, rel_tol=1e-9, abs_tol=0.0)
    assert math.isclose(p3c_of_zeta0, p3c_of_zeta0_d, rel_tol=1e-9, abs_tol=0.0)
    assert math.isclose(jac_p3c_of_zeta0, jac_p3c_of_zeta0_d, rel_tol=1e-9, abs_tol=0.0)
    assert math.isclose(p3c_of_zetatilde_0, p3c_of_zetatilde_0_d, rel_tol=1e-9, abs_tol=0.0)
    assert math.isclose(jac_p3c_of_zetatilde_0, jac_p3c_of_zetatilde_0_d, rel_tol=1e-9, abs_tol=0.0)

    # a, b, c, p2c, jac = build(7, 3)
    # assert a == 10
    # assert b == 4
    # assert c == 11
    # assert p2c(2) == a * 4 + b * 2 + c
    # assert jac(2) == 2 * a * 2 + b
    
    