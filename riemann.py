from sympy import *
from itertools import product

def christoffel(coords, g, g_inv, mu, alpha, beta):
    """ Calculates the christoffel symbol Gamma^{mu}_{alpha beta}, given metric
    g and its inverse, as well as a coordinate system
    """
    res = 0
    for nu in coords:
        res += (diff(g(nu, alpha), beta) +
                diff(g(nu, beta),  alpha) -
                diff(g(alpha, beta), nu)) * g_inv(mu, nu) / 2
    return res

def riemann_tensor(coords, g, g_inv, mu, nu, alpha, beta):
    """ Calculates the Riemann Curvature tensor, R_{mu nu alpha beta}, given
    metric g and its inverse, as well as a coordinate system
    """
    c = make_christoffel(coords, g, g_inv)
    res = diff(c(mu, beta, nu), alpha) - diff(c(mu, alpha, nu), beta)
    for sigma in coords:
        res += c(mu, alpha, sigma) * c(sigma, beta,  nu) - \
               c(mu, beta,  sigma) * c(sigma, alpha, nu)
    return res * g(mu, mu)

def make_riemann(coords, g, g_inv):
    return lambda mu, nu, alpha, beta: \
      riemann_tensor(coords, g, g_inv, mu, nu, alpha, beta)

def make_christoffel(coords, g, g_inv):
    return lambda mu, alpha, beta: \
      christoffel(coords, g, g_inv, mu, alpha, beta)

def get_christoffels(coords, g, g_inv):
    """ Prints all non-zero christoffel symbols """
    for mu, alpha, beta in product(coords, repeat=3):
        res = christoffel(coords, g, g_inv, mu, alpha, beta)
        if res != 0:
            print mu, alpha, beta, "\t: ", simplify(res)

def get_riemanns(coords, g, g_inv):
    """ Prints all non-zero components of Riemann curvature tensor """
    for mu, nu, alpha, beta in product(coords, repeat=4):
        res = riemann_tensor(coords, g, g_inv, mu, nu, alpha, beta)
        if res != 0:
            print mu, nu, alpha, beta, "\t: ", simplify(res)

# Define symbols
t, r, th, ph= symbols('t r th ph')
b, a, m = symbols('b a m')
spherical = (t, r, th, ph)
spherical_2d = (th, ph)

def g_wormhole(mu, nu):
    """ Wormhole metric """
    a = (mu, nu)
    if a == (t, t):
        return -1
    elif a == (r, r):
        return 1
    elif a == (th, th):
        return b**2 + r**2
    elif a == (ph, ph):
        return (b**2 + r**2) * sin(th)**2
    else:
        return 0

def g_schwarz(mu, nu):
    """ Schwarzchild metric """
    a = (mu, nu)
    if a == (t, t):
        return 2 * m / r - 1
    elif a == (r, r):
        return 1 / (1 - 2 * m / r)
    elif a == (th, th):
        return r**2
    elif a == (ph, ph):
        return (r * sin(th))**2
    else:
        return 0

def g_sphere(mu,nu):
    """ 2-D spherical metric """
    if (mu, nu) == (th, th):
        return a ** 2
    elif (mu, nu) == (ph, ph):
        return (a * sin(th)) ** 2
    else:
        return 0

def inv_diag(g):
    """Inverts diagonal metrics"""
    def g_inv(mu, nu):
        if mu == nu:
            return 1 / g(mu, mu)
        else:
            return 0
    return g_inv

def main():
    #get_riemanns(spherical, g_wormhole, inv_diag(g_wormhole_inv))
    #get_christoffels(spherical, g_schwarz, inv_diag(g_schwarz))
    #get_christoffels(spherical_2d, g_sphere, inv_diag(g_sphere))
    get_riemanns(spherical, g_schwarz, inv_diag(g_schwarz))

R = make_riemann(spherical, g_schwarz, inv_diag(g_schwarz))
g = g_schwarz
gi = inv_diag(g_schwarz)
