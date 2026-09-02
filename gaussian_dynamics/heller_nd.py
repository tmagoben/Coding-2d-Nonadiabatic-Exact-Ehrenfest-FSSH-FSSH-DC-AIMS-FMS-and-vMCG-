import numpy as np


def initial_heller_nd(q0, p0, sigma):
    """Initialize multidimensional Heller parameters.

    sigma may be scalar or length-D vector of coordinate standard deviations.
    The Heller convention is
        psi = exp{i[1/2 xi^T A xi + p^T xi + gamma]}
    with Im(A) positive definite.
    """
    q0 = np.asarray(q0, dtype=float)
    p0 = np.asarray(p0, dtype=float)

    if q0.shape != p0.shape:
        raise ValueError("q0 and p0 must have the same shape.")

    D = len(q0)
    sigma = np.asarray(sigma, dtype=float)
    if sigma.ndim == 0:
        sigma = np.full(D, float(sigma))
    if sigma.shape != (D,) or np.any(sigma <= 0):
        raise ValueError("sigma must be positive and scalar or length D.")

    A0 = 1j * np.diag(1.0/(2.0*sigma**2))

    # Product of 1D normalized Gaussian prefactors.
    N0 = np.prod((1.0/(2.0*np.pi*sigma**2))**0.25)
    gamma0 = -1j*np.log(N0)

    return q0, p0, A0, complex(gamma0)


def _pack(q, p, A, gamma):
    D = len(q)
    return np.concatenate([
        np.asarray(q, complex),
        np.asarray(p, complex),
        np.asarray(A, complex).reshape(D*D),
        np.array([gamma], complex),
    ])


def _unpack(z, D):
    q = z[:D].real
    p = z[D:2*D].real
    A = z[2*D:2*D+D*D].reshape(D, D)
    gamma = z[-1]
    return q, p, A, gamma


def _rhs(z, D, mass, potential, gradient, hessian):
    q, p, A, gamma = _unpack(z, D)

    dq = p/mass
    dp = -np.asarray(gradient(q), dtype=float)
    H = np.asarray(hessian(q), dtype=float)
    dA = -(A @ A)/mass - H
    dgamma = (p @ p)/(2.0*mass) - potential(q) + 0.5j*np.trace(A)/mass

    return _pack(dq, dp, dA, dgamma)


def run_thawed_gaussian_nd(
    q0,
    p0,
    sigma,
    mass,
    potential,
    gradient,
    hessian,
    dt=0.001,
    steps=1000,
    store_every=1,
):
    """Multidimensional Heller TGA with complex symmetric width matrix."""
    q, p, A, gamma = initial_heller_nd(q0, p0, sigma)
    D = len(q)
    z = _pack(q, p, A, gamma)

    times=[]; qs=[]; ps=[]; As=[]; gammas=[]

    def rhs(zz):
        return _rhs(zz, D, mass, potential, gradient, hessian)

    def record(step):
        q,p,A,gamma = _unpack(z,D)
        times.append(step*dt)
        qs.append(q.copy())
        ps.append(p.copy())
        As.append(A.copy())
        gammas.append(gamma)

    record(0)

    for step in range(1, steps+1):
        k1=rhs(z)
        k2=rhs(z+0.5*dt*k1)
        k3=rhs(z+0.5*dt*k2)
        k4=rhs(z+dt*k3)
        z = z + dt*(k1+2*k2+2*k3+k4)/6.0

        # The exact equation preserves A=A^T. Remove only roundoff antisymmetry.
        q,p,A,gamma = _unpack(z,D)
        A = 0.5*(A+A.T)
        z = _pack(q,p,A,gamma)

        if step % store_every == 0:
            record(step)

    return {
        "time":np.asarray(times),
        "q":np.asarray(qs),
        "p":np.asarray(ps),
        "A":np.asarray(As),
        "gamma":np.asarray(gammas),
    }
