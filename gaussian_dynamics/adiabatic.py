import numpy as np

from .potentials import avoided_crossing_diabatic


def avoided_crossing_derivative(x, k=0.02, shift=2.0, coupling=0.01, beta=1.0):
    x = float(x)
    v12 = coupling*np.exp(-beta*x*x)
    return np.array([
        [k*(x+shift), -2.0*beta*x*v12],
        [-2.0*beta*x*v12, k*(x-shift)],
    ], dtype=float)


def adiabatic_point(x, model=avoided_crossing_diabatic, derivative=avoided_crossing_derivative,
                    reference=None):
    V = model(float(x))
    dV = derivative(float(x))
    E, U = np.linalg.eigh(V)

    if reference is not None:
        U = U.copy()
        for a in range(U.shape[1]):
            if np.dot(reference[:,a], U[:,a]) < 0:
                U[:,a] *= -1.0

    G = U.T @ dV @ U
    gradients = np.diag(G).copy()
    d = np.zeros((2,2), dtype=float)
    d[0,1] = G[0,1]/(E[1]-E[0])
    d[1,0] = -d[0,1]
    return E, U, gradients, d


def adiabatic_grid(x, model=avoided_crossing_diabatic, derivative=avoided_crossing_derivative):
    E = np.zeros((len(x),2))
    U = np.zeros((len(x),2,2))
    grad = np.zeros((len(x),2))
    d = np.zeros((len(x),2,2))

    ref = None
    for i, xi in enumerate(x):
        Ei, Ui, gi, di = adiabatic_point(xi, model, derivative, reference=ref)
        E[i], U[i], grad[i], d[i] = Ei, Ui, gi, di
        ref = Ui
    return E, U, grad, d


def derivative_fd(values, dx):
    return np.gradient(values, dx, axis=0, edge_order=2)


def adiabatic_hamiltonian_action(chi, x, dx, mass, E, d):
    """Apply H_ad = -1/(2M)(partial+d)^2 + E using finite differences."""
    chi = np.asarray(chi, dtype=complex)
    chip = derivative_fd(chi, dx)
    chipp = derivative_fd(chip, dx)

    dp = derivative_fd(d, dx)
    d2 = np.einsum("xik,xkj->xij", d, d)
    tau2 = dp + d2

    first = 2.0*np.einsum("xij,xj->xi", d, chip)
    second = np.einsum("xij,xj->xi", tau2, chi)

    return -(chipp + first + second)/(2.0*mass) + E*chi


def diabatic_hamiltonian_action(psi, x, dx, mass, V):
    pp = derivative_fd(derivative_fd(psi,dx), dx)
    return -pp/(2.0*mass) + np.einsum("xij,xj->xi", V, psi)


def transform_ad_to_dia(chi, U):
    return np.einsum("xij,xj->xi", U, chi)


def transform_dia_to_ad(psi, U):
    return np.einsum("xji,xj->xi", U.conj(), psi)
