"""Pedagogical nonadiabatic quantum-dynamics algorithms."""

from .adiabatic import (
    diagonalize_path,
    finite_difference_derivative_couplings_1d,
    hellmann_feynman_derivative_couplings,
)
from .diabatization import integrate_adt_path, transform_adiabatic_to_diabatic
from .direct import (
    build_direct_hamiltonian_1d,
    build_direct_hamiltonian_2d,
    diagonalize_hamiltonian,
    propagate_from_eigendecomposition,
)
from .grids import PeriodicGrid1D, PeriodicGrid2D
from .models import (
    linear_vibronic_coupling_2d,
    lvc_analytic_derivative_coupling,
    lvc_analytic_eigensystem,
    smooth_single_avoided_crossing,
    tully_single_avoided_crossing,
)
from .observables import (
    adiabatic_populations_1d,
    adiabatic_populations_2d,
    diabatic_populations,
    norm,
)
from .propagators import (
    potential_propagator_2x2,
    split_operator_step_1d,
    split_operator_step_2d,
)
from .wavepackets import (
    gaussian_wavepacket_1d,
    gaussian_wavepacket_2d,
    prepare_adiabatic_wavepacket_1d,
    prepare_adiabatic_wavepacket_2d,
)

__all__ = [
    "PeriodicGrid1D",
    "PeriodicGrid2D",
    "smooth_single_avoided_crossing",
    "tully_single_avoided_crossing",
    "linear_vibronic_coupling_2d",
    "lvc_analytic_derivative_coupling",
    "lvc_analytic_eigensystem",
    "diagonalize_path",
    "finite_difference_derivative_couplings_1d",
    "hellmann_feynman_derivative_couplings",
    "integrate_adt_path",
    "transform_adiabatic_to_diabatic",
    "potential_propagator_2x2",
    "split_operator_step_1d",
    "split_operator_step_2d",
    "build_direct_hamiltonian_1d",
    "build_direct_hamiltonian_2d",
    "diagonalize_hamiltonian",
    "propagate_from_eigendecomposition",
    "gaussian_wavepacket_1d",
    "gaussian_wavepacket_2d",
    "prepare_adiabatic_wavepacket_1d",
    "prepare_adiabatic_wavepacket_2d",
    "norm",
    "diabatic_populations",
    "adiabatic_populations_1d",
    "adiabatic_populations_2d",
]
