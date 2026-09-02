import numpy as np


def matrix_exponential_hermitian(V, dt):
    V = np.asarray(V, dtype=complex)
    out = np.empty_like(V)
    for ix in range(V.shape[0]):
        e, U = np.linalg.eigh(V[ix])
        out[ix] = (U * np.exp(-1j*dt*e)) @ U.conj().T
    return out


def apply_pointwise_matrix(U, psi):
    return np.einsum("xij,xj->xi", U, psi)


def multistate_split_step(psi, x, dx, dt, mass, V):
    Uhalf = matrix_exponential_hermitian(V, 0.5*dt)
    psi = apply_pointwise_matrix(Uhalf, psi)

    k = 2*np.pi*np.fft.fftfreq(len(x), d=dx)
    phase = np.exp(-0.5j*dt*k**2/mass)
    psi = np.fft.ifft(phase[:,None]*np.fft.fft(psi, axis=0), axis=0)

    psi = apply_pointwise_matrix(Uhalf, psi)
    return psi


def multistate_norm(psi, dx):
    return float(np.sum(np.abs(psi)**2)*dx)


def run_multistate_exact(psi0, x, dx, V, mass=1.0, dt=0.002, steps=1000, store_every=10):
    psi = np.asarray(psi0, dtype=complex).copy()
    psi /= np.sqrt(multistate_norm(psi, dx))

    times=[]; states=[]; norms=[]; pops=[]
    def record(step):
        times.append(step*dt)
        states.append(psi.copy())
        norms.append(multistate_norm(psi,dx))
        pops.append(np.sum(np.abs(psi)**2, axis=0)*dx)

    record(0)
    for step in range(1,steps+1):
        psi = multistate_split_step(psi,x,dx,dt,mass,V)
        if step % store_every == 0:
            record(step)

    return {
        "time":np.asarray(times),
        "psi":np.asarray(states),
        "norm":np.asarray(norms),
        "populations":np.asarray(pops),
    }
