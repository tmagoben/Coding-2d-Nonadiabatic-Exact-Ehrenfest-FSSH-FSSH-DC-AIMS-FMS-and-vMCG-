import numpy as np


def pointwise_matrix_exponential_2d(V, dt):
    """exp(-i V dt) for V[nx,ny,nstate,nstate], using batched eigh."""
    V = np.asarray(V, dtype=complex)
    e, U = np.linalg.eigh(V)
    phase = np.exp(-1j*dt*e)

    # U[...,i,k] phase[...,k] U*[...,j,k]
    return np.einsum("...ik,...k,...jk->...ij", U, phase, U.conj())


def apply_pointwise_matrix_2d(U, psi):
    return np.einsum("...ij,...j->...i", U, psi)


def split_operator_2d_step(psi, dx, dy, dt, mass, Uhalf):
    """One 2D multistate Strang step with precomputed half-step potential."""
    psi = apply_pointwise_matrix_2d(Uhalf, psi)

    nx, ny = psi.shape[:2]
    kx = 2.0*np.pi*np.fft.fftfreq(nx, d=dx)
    ky = 2.0*np.pi*np.fft.fftfreq(ny, d=dy)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    kinetic = np.exp(-0.5j*dt*(KX**2+KY**2)/mass)

    psi_k = np.fft.fftn(psi, axes=(0,1))
    psi = np.fft.ifftn(kinetic[...,None]*psi_k, axes=(0,1))

    psi = apply_pointwise_matrix_2d(Uhalf, psi)
    return psi


def norm_2d(psi, dx, dy):
    return float(np.sum(np.abs(psi)**2)*dx*dy)


def run_exact_2d(
    psi0,
    dx,
    dy,
    V,
    mass=20.0,
    dt=0.002,
    steps=200,
    store_every=10,
):
    """Exact-grid 2D two-state propagation for a time-independent V(x,y)."""
    psi = np.asarray(psi0, dtype=complex).copy()
    psi /= np.sqrt(norm_2d(psi,dx,dy))

    Uhalf = pointwise_matrix_exponential_2d(V, 0.5*dt)

    times=[]; states=[]; norms=[]; populations=[]

    def record(step):
        times.append(step*dt)
        states.append(psi.copy())
        norms.append(norm_2d(psi,dx,dy))
        populations.append(np.sum(np.abs(psi)**2, axis=(0,1))*dx*dy)

    record(0)
    for step in range(1,steps+1):
        psi = split_operator_2d_step(psi,dx,dy,dt,mass,Uhalf)
        if step % store_every == 0:
            record(step)

    return {
        "time":np.asarray(times),
        "psi":np.asarray(states),
        "norm":np.asarray(norms),
        "populations":np.asarray(populations),
    }
