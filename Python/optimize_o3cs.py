import numpy as np
# from scipy.integrate import fixed_quad as integrate
from scipy.integrate import quadrature as integrate
# from scipy.integrate import quad as integrate
import pylab as pl


def o3cs(alpha=0.5, gamma=0.5):
    a = -2 - 4 * alpha
    b = -4 * a

    def integrand(kh):
        keq3_k3 = (
            (8 * a * np.sin(kh) + b * np.sin(2 * kh)) /
            (8 * kh ** 3 * (1 + 2 * alpha * np.cos(kh))))

        return np.abs(1 + keq3_k3)

    return integrate(integrand, 0, gamma * np.pi)


def optimize(objective_fn):
    alpha = np.arange(0.44, 0.49, 1e-4)
    gamma = np.arange(0.5, 1.0, 0.04)

    dict_obj = {}
    for g in gamma:
        obj = np.zeros_like(alpha)
        for i, a in enumerate(alpha):
            value = objective_fn(a, g)[0]
            obj[i] = value

        dict_obj[g] = obj

    return alpha, dict_obj


def plot(alpha, dict_obj):
    for k, v in dict_obj.items():
        pl.semilogy(alpha, v)

    pl.legend(dict_obj.keys())
    pl.show()


if __name__ == '__main__':
    a, d = optimize(o3cs)
    plot(a, d)
