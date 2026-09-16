from __future__ import annotations

import numpy as np


LOG2 = float(np.log(2.0))


def repair_isolates_deterministic(
    filtered: np.ndarray,
    adjacency: np.ndarray,
    dissimilarity: np.ndarray,
) -> np.ndarray:
    """Reconnect each isolate to its least-dissimilar original neighbor."""
    repaired = np.asarray(filtered, dtype=np.int8).copy()
    adjacency = np.asarray(adjacency)
    dissimilarity = np.asarray(dissimilarity, dtype=float)
    degree = repaired.sum(axis=1).astype(int)

    for node in range(repaired.shape[0]):
        if degree[node] != 0:
            continue
        neighbors = np.flatnonzero(adjacency[node] > 0.5)
        if neighbors.size == 0:
            continue
        neighbor = neighbors[np.argmin(dissimilarity[node, neighbors])]
        if repaired[node, neighbor] == 0:
            repaired[node, neighbor] = 1
            repaired[neighbor, node] = 1
            degree[node] += 1
            degree[neighbor] += 1
    return repaired


def effective_boundary_draws(
    eta_draws: np.ndarray,
    adjacency: np.ndarray,
    dissimilarity: np.ndarray,
    edge_i: np.ndarray | None = None,
    edge_j: np.ndarray | None = None,
    threshold: float = LOG2,
) -> np.ndarray:
    """Return post-repair boundary indicators for every draw and observed edge."""
    eta_draws = np.asarray(eta_draws, dtype=float).reshape(-1)
    adjacency = np.asarray(adjacency)
    dissimilarity = np.asarray(dissimilarity, dtype=float)
    if adjacency.shape != dissimilarity.shape or adjacency.ndim != 2:
        raise ValueError("adjacency and dissimilarity must be square matrices of equal shape.")

    if edge_i is None or edge_j is None:
        edge_i, edge_j = np.where(np.triu(adjacency > 0.5, 1))
    edge_i = np.asarray(edge_i, dtype=int).reshape(-1)
    edge_j = np.asarray(edge_j, dtype=int).reshape(-1)
    edge_z = dissimilarity[edge_i, edge_j]

    retained = edge_z[None, :] * eta_draws[:, None] <= threshold
    incidence = np.zeros((edge_i.size, adjacency.shape[0]), dtype=np.int8)
    incidence[np.arange(edge_i.size), edge_i] = 1
    incidence[np.arange(edge_j.size), edge_j] = 1
    degree = retained.astype(np.int16) @ incidence

    edge_lookup = np.full(adjacency.shape, -1, dtype=int)
    edge_lookup[edge_i, edge_j] = np.arange(edge_i.size)
    edge_lookup[edge_j, edge_i] = np.arange(edge_i.size)

    # Follow the simulator's node order: an earlier repair can de-isolate a later node.
    for node in range(adjacency.shape[0]):
        isolated_draws = degree[:, node] == 0
        if not np.any(isolated_draws):
            continue
        neighbors = np.flatnonzero(adjacency[node] > 0.5)
        if neighbors.size == 0:
            continue
        neighbor = neighbors[np.argmin(dissimilarity[node, neighbors])]
        edge_position = edge_lookup[node, neighbor]
        if edge_position < 0:
            raise ValueError("The edge list does not contain every adjacency edge.")
        retained[isolated_draws, edge_position] = True
        degree[isolated_draws, node] += 1
        degree[isolated_draws, neighbor] += 1

    return ~retained


def effective_boundary_truth(
    filtered_adjacency: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
) -> np.ndarray:
    filtered_adjacency = np.asarray(filtered_adjacency)
    edge_i = np.asarray(edge_i, dtype=int).reshape(-1)
    edge_j = np.asarray(edge_j, dtype=int).reshape(-1)
    return (filtered_adjacency[edge_i, edge_j] < 0.5).astype(np.int8)
