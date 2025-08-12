"""
Example script demonstrating how to use the `p3c` module to build a cubic polynomial transformation.
It also shows how to evaluate the polynomial and its Jacobian at specific points.
"""

from p3c import build_p3c

# Interpolation points of the boundary element.
# Three points are always required, even for centered straight elements.
X1, Y1 =   1.,   1.
X2, Y2 =   3.,   2.5
X3, Y3 =   5.,   2. 

# Collocation point
Xc, Yc = 4.41, 2.13

# Compute the cubic polynomial transformation
(   zeta0, eta0, zeta1, eta1,
    zetatilde_0, etatilde_0, the_n_iteraciones, the_xerr,
    coef_a, coef_b, coef_c, coef_d,
    p3c_of_xitilde, jac_p3c_of_xitilde
) = build_p3c(Xc, Yc, X1, Y1, X2, Y2, X3, Y3)


print("\nComputing the cubic polynomial p3c...")
print(  "-------------------------------------")
    
print("\nNormalized collocation point in the complex plane: ")
print("    zeta0, eta0 = ", zeta0, eta0)
print("    zeta1, eta1 = ", zeta1, eta1, " (secondary solution in the quadratic element case—not necessary)")

print("\nTransformed pole xi_tilde, number of iterations needed, and error achieved:")
print("    zetatilde_0, etatilde_0 = ", zetatilde_0, etatilde_0)
print("    Number of iterations = ", the_n_iteraciones)
print("    Error = ", the_xerr)

print("\nCoefficients of the polynomial p3c: a·x^3 + b·x^2 + c·x + d")
print("    a coef. =", coef_a)
print("    b coef. =:", coef_b)
print("    c coef. =:", coef_c)
print("    d coef. =:", coef_d)


print("\nExamples of evaluating the p3c polynomial and its Jacobian at different points:")

print("\nEvaluation of p3c polynomial and its Jacobian at xitilde = zeta0:")
print("    p3c(zeta0) = ", p3c_of_xitilde(zeta0))
print("    jac_p3c(zeta0) = ", jac_p3c_of_xitilde(zeta0))

print("\nEvaluation of p3c polynomial and its Jacobian at xitilde = zetatilde_0:")
print("    p3c(zetatilde_0) = ", p3c_of_xitilde(zetatilde_0)," (must be equal to zeta0)")
print("    jac_p3c(zetatilde_0) = ", jac_p3c_of_xitilde(zetatilde_0), "(should be the minimum of jac_p3c)" )

print("\nChecking the minimum: jac_p3c at zetatilde_0-0.001, zetatilde_0, zetatilde_0+0.001:")
print("    jac_p3c(zetatilde_0-0.001) = ", jac_p3c_of_xitilde(zetatilde_0-0.001))
print("    jac_p3c(zetatilde_0      ) = ", jac_p3c_of_xitilde(zetatilde_0))
print("    jac_p3c(zetatilde_0+0.001) = ", jac_p3c_of_xitilde(zetatilde_0+0.001))




    