from sympy import *
from itertools import product

def christoffel(coords, g, g_inv, mu, alpha, beta):
    """Calculates the christoffel symbol Gamma^{mu}_{alpha beta}, given metric
    g and its inverse, as well as a coordinate system
    """
    res = 0
    for nu in coords:
        res += (diff(g(nu, alpha), beta) +
                diff(g(nu, beta),  alpha) -
                diff(g(alpha, beta), nu)) * g_inv(mu, nu) / 2
    return res

def riemann_tensor(coords, g, g_inv, mu, nu, alpha, beta):
    """Calculates the Riemann Curvature tensor, R_{mu nu alpha beta}, given
    metric g and its inverse, as well as a coordinate system
    """
    c = lambda mu, alpha, beta: christoffel(coords, g, g_inv, mu, alpha, beta)
    res = diff(c(mu, beta, nu), alpha) - diff(c(mu, alpha, nu), beta)
    for sigma in coords:
        res += c(mu, alpha, sigma) * c(sigma, beta,  nu) - \
               c(mu, beta,  sigma) * c(sigma, alpha, nu)
    return res * g(mu, mu)

def get_christoffels(coords, g, g_inv):
    for mu, alpha, beta in product(coords, repeat=3):
        res = christoffel(coords, g, g_inv, mu, alpha, beta)
        if res != 0:
            print mu, alpha, beta, "\t: ", res

def get_riemann(coords, g, g_inv):
    for mu, nu, alpha, beta in product(coords, repeat=4):
        res = riemann_tensor(coords, g, g_inv, mu, nu, alpha, beta)
        if res != 0:
            print mu, nu, alpha, beta, "\t: ", simplify(res)

def g_wormhole(mu, nu):
    a = (mu, nu)
    if a == (t, t):
        return -1
    elif a == (r, r):
        return 1
    elif a == (th, th):
        return b*b + r*r
    elif a == (ph, ph):
        return (b*b + r*r) * sin(th) * sin(th)
    else:
        return 0

def g_wormhole_inv(mu, nu):
    if mu == nu:
        return 1 / g_wormhole(mu, mu)
    else:
        return 0

def g_sphere(mu,nu):
    if (mu, nu) == (th, th):
        return a ** 2
    elif (mu, nu) == (ph, ph):
        return (a * sin(th)) ** 2
    else:
        return 0

def g_sphere_inv(mu, nu):
    if mu == nu:
        return 1 / g_sphere(mu, mu)
    else:
        return 0

t, r, th, ph, b, a= symbols('t r th ph b a')
spherical = (t, r, th, ph)
spherical_2d = (th, ph)
#get_riemann(spherical, g_wormhole, g_wormhole_inv)
get_christoffels(spherical_2d, g_sphere, g_sphere_inv)
