import numpy as np

from .direct_dynamics_nd import run_backend_spawned_gaussians
from .molecular_gauge_graph_v19 import (
    build_molecular_centroid_graph_v19,
)


def run_molecular_direct_dynamics_v19(
    initial_basis,
    C0,
    provider,
    *,
    dt=0.0005,
    steps=100,
    spawn_threshold=1e-4,
    overlap_block=0.85,
    max_basis=8,
    store_every=5,
):
    """Provider-neutral molecular direct-dynamics prototype.

    The propagation kernel is the inherited v0.5 local constant-electronic-quantity
    Gaussian approximation. v0.19 adds the tracked/cached molecular provider and
    records a center-centroid gauge-graph audit of the final basis.

    This is not full AIMS matrix-element theory.
    """
    before=provider.diagnostics_dict()
    out=run_backend_spawned_gaussians(
        initial_basis,C0,provider,
        dt=dt,steps=steps,
        spawn_threshold=spawn_threshold,
        overlap_block=overlap_block,
        max_basis=max_basis,
        store_every=store_every,
    )
    after=provider.diagnostics_dict()

    graph_audit=None
    try:
        graph=build_molecular_centroid_graph_v19(
            out["final_basis"],provider
        )
        Sg,Hg=graph.matrices()
        graph_audit={
            "basis_size":len(out["final_basis"]),
            "graph_nodes":len(graph.registry.graph.nodes),
            "graph_edges":len(graph.registry.graph.edges()),
            "S_hermiticity_error":float(
                np.linalg.norm(Sg-Sg.conj().T,ord="fro")
            ),
            "H_hermiticity_error":float(
                np.linalg.norm(Hg-Hg.conj().T,ord="fro")
            ),
            "condition_number":float(np.linalg.cond(Sg)),
        }
    except Exception as exc:
        graph_audit={
            "failed":True,
            "exception":type(exc).__name__,
            "message":str(exc),
        }

    out["provider_diagnostics_before"]=before
    out["provider_diagnostics_after"]=after
    out["molecular_centroid_graph_audit"]=graph_audit
    return out
