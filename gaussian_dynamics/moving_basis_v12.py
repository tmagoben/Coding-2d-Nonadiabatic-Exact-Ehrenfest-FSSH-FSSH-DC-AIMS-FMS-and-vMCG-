import numpy as np


def moving_basis_midpoint_cayley_step(
    C,
    S_old,
    H_old,
    S_new,
    H_new,
    T_mid,
    dt,
):
    r"""Implicit midpoint/Cayley step for the nonorthogonal moving-basis equation.

    Starting from

        i S Cdot = (H - i T) C,

    define

        K = i H + T

    so that

        S Cdot = -K C.

    Midpoint discretization gives

        [S_m + dt K_m/2] C_{n+1}
          = [S_m - dt K_m/2] C_n.

    For a fixed orthonormal basis this reduces to the ordinary Cayley/Crank-Nicolson
    unitary step.
    """
    C = np.asarray(C, dtype=complex)
    S0 = np.asarray(S_old, dtype=complex)
    S1 = np.asarray(S_new, dtype=complex)
    H0 = np.asarray(H_old, dtype=complex)
    H1 = np.asarray(H_new, dtype=complex)
    T = np.asarray(T_mid, dtype=complex)
    dt = float(dt)

    if dt <= 0.0:
        raise ValueError("dt must be positive.")
    if S0.shape != S1.shape or S0.shape != H0.shape or S0.shape != H1.shape:
        raise ValueError("S/H matrices must have equal shapes.")
    if T.shape != S0.shape or S0.ndim != 2 or S0.shape[0] != S0.shape[1]:
        raise ValueError("T and S/H matrices must be equal-size square matrices.")
    if C.shape != (S0.shape[0],):
        raise ValueError("Coefficient vector has incompatible size.")

    Sm = 0.5*(S0 + S1)
    Hm = 0.5*(H0 + H1)
    K = 1j*Hm + T

    lhs = Sm + 0.5*dt*K
    rhs = (Sm - 0.5*dt*K) @ C

    return np.linalg.solve(lhs, rhs)


def fixed_basis_cayley_operator(H, dt):
    H = np.asarray(H, dtype=complex)
    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise ValueError("H must be square.")
    if not np.allclose(H, H.conj().T, atol=1e-12):
        raise ValueError("H must be Hermitian.")

    I = np.eye(H.shape[0], dtype=complex)
    return np.linalg.solve(
        I + 0.5j*dt*H,
        I - 0.5j*dt*H,
    )


def endpoint_generalized_norm_error(C_old, S_old, C_new, S_new):
    n0 = float(np.real(np.vdot(C_old, np.asarray(S_old) @ C_old)))
    n1 = float(np.real(np.vdot(C_new, np.asarray(S_new) @ C_new)))
    return n1-n0


def phase_aligned_vector_error(candidate, reference):
    """Euclidean error after removing one physically irrelevant global phase."""
    a = np.asarray(candidate, dtype=complex)
    b = np.asarray(reference, dtype=complex)
    if a.shape != b.shape:
        raise ValueError("vectors must have equal shape.")

    overlap = np.vdot(b, a)
    if abs(overlap) > 1e-15:
        a = a * np.conj(overlap)/abs(overlap)

    return float(np.linalg.norm(a-b))
