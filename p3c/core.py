import math
import cmath
from typing import Callable


def polo_normalizado_plano_complejo(Xc,Yc, X1,Y1, X2,Y2, X3,Y3):
    """
    Computes the normalized pole coordinates (zeta_c, eta_c) in the complex plane
    from three interpolation points of a boundary element. Handles both quadratic
    and linear cases.
    """    

    Z1, Z2, Z3, Zc = complex(X1,Y1), complex(X2,Y2), complex(X3,Y3), complex(Xc,Yc)
    a = (Z1 - 2*Z2 + Z3)/2
    b = (Z3-Z1)/2
    c = Z2 - Zc

    if a != 0: #  Quadratic equation case
        discriminant_sqrt = cmath.sqrt(b * b - 4 * a * c)
        chi1_sol = (-b + discriminant_sqrt) / (2 * a)
        chi2_sol = (-b - discriminant_sqrt) / (2 * a)
                
        # --- Selecting the correct solution ---
        # Step 1: Compute the minimum distance from chi1_sol to the segment [-1, 1]
        if abs(chi1_sol.real)  <= 1:
            dist_chi1 = abs(chi1_sol.imag)
        else:
            dist_chi1 = abs(chi1_sol - complex(math.copysign(1, chi1_sol.real), 0))
        
        # Step 2: Compute the minimum distance from chi2_sol to the segment [-1, 1]
        if abs(chi2_sol.real)  <= 1:
            dist_chi2 = abs(chi2_sol.imag)
        else:
            dist_chi2 = abs(chi2_sol - complex(math.copysign(1, chi2_sol.real), 0))
        
        # Step 3: Choose the one farther from the segment
        if dist_chi1 < dist_chi2:
            chi1_sol_es_la_mas_cercana_al_interv_m1_1 = True
        elif dist_chi1 > dist_chi2:
            chi1_sol_es_la_mas_cercana_al_interv_m1_1 = False
        else:   # Tie: choose the one closest to (0, 0)
            if abs(chi1_sol.real) < abs(chi2_sol.real):
                chi1_sol_es_la_mas_cercana_al_interv_m1_1 = True
            else:
                chi1_sol_es_la_mas_cercana_al_interv_m1_1 = False

        # Assign the correct root
        if chi1_sol_es_la_mas_cercana_al_interv_m1_1:
            chi_sol_0 = chi1_sol
            chi_sol_1 = chi2_sol
        else:
            chi_sol_0 = chi2_sol
            chi_sol_1 = chi1_sol

    elif b != 0:  # Linear equation case, centered straight element
        chi_sol_0 = -c/b
        chi_sol_1 = math.nan
    else:
        raise ValueError("Unsolvable equation: boundary element not valid.")

    zeta_c = chi_sol_0.real
    eta_c  = chi_sol_0.imag

    return zeta_c,eta_c, chi_sol_1.real, chi_sol_1.imag
   



# BEGIN: Computing the transformed pole (with tilde)

def NR(f, fp, xini, xtol, *args): # Newton Raphson's method
    """Newton-Raphson iteration."""
    x0 = xini
    nmaxit = 500
    for number in range(nmaxit):
        fx0  = f(x0, *args)
        fpx0 = fp(x0, *args)
        x1 = x0 - fx0/fpx0
        xerr = abs(x1 - x0)
        # Stop if tolerance is reached
        if xerr < xtol: break 
        if number == nmaxit: 
            raise RuntimeError(
                f"Maximum number of iterations ({nmaxit}) reached. "
                "No convergence in Newton-Raphson method."
            )
        x0 = x1
    return x1, number+1, xerr


def ftfetc(xetilde_c,xseta_c,xeta_c):  # Function to find etilde_c
     return ( (2*xseta_c/xeta_c*xetilde_c**3)**2 - 
             ( 2/3*xetilde_c**3/xeta_c + 2*xetilde_c**2 + 8/3 )**2 *(2/3*xetilde_c**3/xeta_c - xetilde_c**2 - 1/3) )
def ftfetcp(xetilde_c,xseta_c,xeta_c): # Derivative of ftfetc
    return ( 24*xetilde_c**5 + 48*xetilde_c**3 + 64/3*xetilde_c - 40/3*xetilde_c**4/xeta_c - 
             32/3*xetilde_c**2/xeta_c - 32/3*xetilde_c**7/xeta_c**2 + 
             24*xetilde_c**5*xseta_c**2/xeta_c**2 - 40/3*xetilde_c**5/xeta_c**2 - 8/3*xetilde_c**8/xeta_c**3 )


def cubica_una_raiz_real(xeta_c):
    """Solve a cubic equation that has exactly one real root."""
    aa = 2
    bb = -3
    cc = 0
    dd = -1/xeta_c**2
    Delta_1 = bb**2-3*aa*cc
    Delta_2 = 2*bb**3-9*aa*bb*cc+27*aa**2*dd

    if (Delta_1 == 0 and Delta_2 == 0): # Multiple real solution. In this particular case:
        raise ValueError(
            "Multiple real solution encountered in cubica_una_raiz_real. "
            "This case is not implemented."
        )
    else:
        Delta_3 = Delta_2**2-4*Delta_1**3
        if Delta_3 < 0: 
            raise ValueError(
                "Three real roots found in cubica_una_raiz_real. "
                "This case is not implemented."
            )            
        else:
            Delta_4 = (Delta_2+(Delta_3)**(1/2))/2
            if Delta_4 == 0:
                Delta_4 = (Delta_2-(Delta_3)**0.5)/2
            if Delta_4 > 0:
                W = Delta_4**(1/3)
            else:    
                W = -((-Delta_4)**(1/3))
            sol = -1/(3*aa)*(bb+W+Delta_1/W)
    return sol


def radicando(xetilde_c,xeta_c):
    return 2/3*xetilde_c**3/xeta_c - xetilde_c**2 - 1/3


def stilde_c_eta_c_cero(xseta_c):
    """Special case solver for stilde when eta_c = 0."""    
    a = 1
    b = -3*xseta_c
    c = 3
    d = -xseta_c
    Delta_1 = b**2-3*a*c
    Delta_2 = 2*b**3-9*a*b*c+27*a**2*d
    if (Delta_1 == 0 and Delta_2 == 0): # Multiple real solution. In this particular case:
        sol = xseta_c
    else:
        Delta_3 = Delta_2**2-4*Delta_1**3
        if Delta_3 > 0: # One real root exists
            Delta_4 = (Delta_2+(Delta_3)**0.5)/2
            if Delta_4 == 0:
                Delta_4 = (Delta_2-(Delta_3)**0.5)/2
            if Delta_4 > 0:
                W = Delta_4**(1/3)
            else:    
                W = -((-Delta_4)**(1/3))
            sol = -1/(3*a)*(b+W+Delta_1/W)
    return sol


def fun_stilde_etilde(xchi): # Main function
    """    
    Function to compute the transformed pole parameters (stilde_c, etilde_c)
    from the input pair (xseta, xeta). Handles special cases analytically and
    uses Newton-Raphson iteration for the general case.
    """        

    xseta,xeta = xchi
    the_n_iteraciones, the_xerr = 0, 0.0
    if xeta == 0.:
        sol_stilde_c = stilde_c_eta_c_cero(xseta)
        sol_etilde_c = xeta
    else:
        # Start at the solution for the case when seta_c = 0        
        etilde_c_ini = cubica_una_raiz_real(xeta)*xeta 

        if xseta == 0.:
            sol_stilde_c = xseta
            sol_etilde_c = etilde_c_ini
        else:
            # First, adjust the initial guess:
            etilde_c_ini += math.cbrt(xeta) * math.hypot(xseta, xeta) # Esta línea es como la precedente pero optimizada # Este es mejor para valores de eta próximos al cero, por lo que tomamos este

            # Apply Newton-Raphson method
            sol_etilde_c, the_n_iteraciones, the_xerr = NR(ftfetc,ftfetcp,etilde_c_ini,1e-12,xseta,xeta)

            # Now compute stilde_0
            sol_stilde_c2 = radicando(sol_etilde_c,xeta)
            sol_stilde_c = math.copysign(math.sqrt(sol_stilde_c2), xseta)


    return sol_stilde_c,sol_etilde_c, the_n_iteraciones,the_xerr

# END: Computing the transformed pole (with tilde)






# BEGIN BLOCK: p3c polynomial and Jacobian coefficient computation

def A_B_p3c(zetatilde_c, etatilde_c):
    """Compute coefficients A_p3c and B_p3c based on transformed coordinates."""    
    A_p3c = 3/(1 +3*zetatilde_c**2 +3*etatilde_c**2)
    B_p3c = (3 +zetatilde_c**2)*zetatilde_c/(1 +3*zetatilde_c**2 +3*etatilde_c**2)
    return A_p3c, B_p3c


def coeficientes_p3c(zetatilde_c, etatilde_c, A_p3c, B_p3c):
    """
    Return the coefficients (aa, bb, cc, dd) for the cubic polynomial:
        p3c(xitilde) = xi(xitilde) = aa·xitilde³ + bb·xitilde² + cc·xitilde + dd
    """
    aa = A_p3c / 3.0
    bb = -A_p3c * zetatilde_c
    cc = A_p3c * (zetatilde_c**2 + etatilde_c**2)
    dd = -A_p3c * zetatilde_c**3 / 3.0 + B_p3c
    return aa, bb, cc, dd
def coeficientes_jactilde_p3c(A_p3c, bb_p3c, cc_p3c):
    """
    Return the coefficients (aa, bb, cc) for the Jacobian of the cubic:
        jac_p3c(xitilde) = xi'(xitilde) = aa·xitilde² + bb·xitilde + cc
    """
    aa_jac = A_p3c
    bb_jac = 2.0 * bb_p3c
    cc_jac = cc_p3c
    return aa_jac, bb_jac, cc_jac

def build_p3c(Xc: float, Yc: float,
              X1: float, Y1: float,
              X2: float, Y2: float,
              X3: float, Y3: float) -> tuple[
                  float, float, float, float,
                  float, float, int, float,
                  float, float, float, float,
                  Callable[[float], float],
                  Callable[[float], float]]:
    """
    Construct the p3c cubic transformation and its Jacobian from three control points.

    Parameters
    ----------
    Xc, Yc : float
        Center point coordinates.
    X1, Y1, X2, Y2, X3, Y3 : float
        Control point coordinates.

    Returns
    -------
    zeta0, eta0, zeta1, eta1 : float
        Normalized pole coordinates in the complex plane.
    zetatilde_0, etatilde_0 : float
        Transformed pole coordinates (with tilde).
    the_n_iteraciones : int
        Number of Newton–Raphson iterations (if used).
    the_xerr : float
        Final iteration error.
    aa, bb, cc, dd : float
        Cubic polynomial coefficients for xi(xitilde).
    p3c_of_xitilde : Callable
        Function to compute p3c(xitilde) = xi(xitilde) from xitilde.
    jac_p3c_of_xitilde : Callable
        Function to compute the Jacobian of p3c: jac_p3c(xitilde) = xi'(xitilde).
    """

    zeta0,eta0, zeta1,eta1 = polo_normalizado_plano_complejo(Xc,Yc, X1,Y1, X2,Y2, X3,Y3) 
    
    zetatilde_0,etatilde_0, the_n_iteraciones, the_xerr = fun_stilde_etilde((zeta0, eta0))

    A_p3c, B_p3c = A_B_p3c(zetatilde_0, etatilde_0)
    aa, bb, cc, dd = coeficientes_p3c(zetatilde_0, etatilde_0, A_p3c, B_p3c)
    aa_jac, bb_jac, cc_jac = coeficientes_jactilde_p3c(A_p3c, bb, cc)

    def p3c_of_xitilde(xitilde: float) -> float:
        """Evaluate the cubic polynomial xi(xitilde)."""        
        return ((aa *xitilde + bb) *xitilde + cc) *xitilde + dd

    def jac_p3c_of_xitilde(xitilde: float) -> float:
        """Evaluate the Jacobian of the cubic polynomial xi(xitilde)."""        
        return (aa_jac *xitilde + bb_jac) *xitilde + cc_jac

    return (zeta0, eta0, zeta1, eta1, 
           zetatilde_0, etatilde_0, the_n_iteraciones, the_xerr,  
           aa, bb, cc, dd, 
           p3c_of_xitilde, jac_p3c_of_xitilde
    )

# END BLOCK: p3c polynomial and Jacobian coefficient computation