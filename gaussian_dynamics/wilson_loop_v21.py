import numpy as np

def wilson_product_v21(links):
    links = [np.asarray(x, dtype=complex) for x in links]
    n = links[0].shape[0]
    W = np.eye(n, dtype=complex)
    for L in links: W = W @ L
    return W

def gauge_transform_cycle_links_v21(links, node_gauges):
    links = [np.asarray(x, dtype=complex) for x in links]; G = [np.asarray(x, dtype=complex) for x in node_gauges]
    return [G[k].conj().T @ L @ G[(k + 1) % len(G)] for k, L in enumerate(links)]

def sorted_wilson_eigenphases_v21(links):
    return np.sort(np.angle(np.linalg.eigvals(wilson_product_v21(links))))
