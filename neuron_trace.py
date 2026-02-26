"""
Neuron tracer - region grow from seed points, skeletonize, filter non-neurons.

Functions:
    trace_neuron(img, seed, channel=0, threshold_pct=0.10)
    trace_all_neurons(img, centroids, nuclei_labeled, channel=0, 
                      threshold_pct=0.10, min_reach_ratio=2.0)
    show_traces(img, traces, nuclei_labeled=None)
"""

import numpy as np
import matplotlib.pyplot as plt
from skimage import morphology, measure
from collections import deque


def region_grow(img_gray, seed, threshold_pct=0.10):
    """Grow mask from seed point based on intensity threshold."""
    h, w = img_gray.shape
    img = img_gray.astype(np.float64)
    min_intensity = img[seed[0], seed[1]] * threshold_pct

    mask = np.zeros((h, w), dtype=bool)
    visited = np.zeros((h, w), dtype=bool)
    queue = deque([seed])
    visited[seed[0], seed[1]] = True

    while queue:
        y, x = queue.popleft()
        if img[y, x] >= min_intensity:
            mask[y, x] = True
            for dy in [-1, 0, 1]:
                for dx in [-1, 0, 1]:
                    ny, nx = y + dy, x + dx
                    if 0 <= ny < h and 0 <= nx < w and not visited[ny, nx]:
                        visited[ny, nx] = True
                        queue.append((ny, nx))
    return mask


def geodesic_max_distance(mask, seed):
    """
    Quick geodesic BFS to find the max distance from seed through mask.
    Returns just the max distance value (not the full map).
    """
    h, w = mask.shape

    if not mask[seed[0], seed[1]]:
        mask_ys, mask_xs = np.where(mask)
        if len(mask_ys) == 0:
            return 0.0
        dists = (mask_ys - seed[0])**2 + (mask_xs - seed[1])**2
        nearest = np.argmin(dists)
        seed = (mask_ys[nearest], mask_xs[nearest])

    dist = np.full((h, w), np.inf)
    dist[seed[0], seed[1]] = 0
    queue = deque([seed])
    max_dist = 0.0

    neighbours = [
        (-1, 0, 1.0), (1, 0, 1.0), (0, -1, 1.0), (0, 1, 1.0),
        (-1, -1, 1.414), (-1, 1, 1.414), (1, -1, 1.414), (1, 1, 1.414),
    ]

    while queue:
        y, x = queue.popleft()
        current_dist = dist[y, x]
        for dy, dx, cost in neighbours:
            ny, nx = y + dy, x + dx
            if 0 <= ny < h and 0 <= nx < w and mask[ny, nx]:
                new_dist = current_dist + cost
                if new_dist < dist[ny, nx]:
                    dist[ny, nx] = new_dist
                    if new_dist > max_dist:
                        max_dist = new_dist
                    queue.append((ny, nx))

    return max_dist


def nucleus_radius(nucleus_mask):
    """Estimate nucleus radius from mask area (assuming circular)."""
    area = nucleus_mask.sum()
    if area == 0:
        return 15.0
    return np.sqrt(area / np.pi)


def trace_neuron(img, seed, channel=0, threshold_pct=0.10):
    """
    Trace a single neuron from a seed point.

    Parameters
    ----------
    img : array
        RGB image.
    seed : tuple (y, x)
        Seed point (e.g. nucleus centroid).
    channel : int
        Channel to trace on (0=red, 1=green, 2=blue).
    threshold_pct : float
        Intensity threshold as fraction of seed value.

    Returns
    -------
    dict with: mask, skeleton, seed
    """
    gray = img[:, :, channel].astype(np.float64)

    mask = region_grow(gray, seed, threshold_pct)
    mask = morphology.remove_small_objects(mask, min_size=50)
    mask = morphology.remove_small_holes(mask, area_threshold=100)

    labeled = measure.label(mask)
    seed_label = labeled[seed[0], seed[1]]
    if seed_label > 0:
        mask = labeled == seed_label
    else:
        return {"mask": mask, "skeleton": np.zeros_like(mask), "seed": seed}

    skeleton = morphology.skeletonize(mask)

    return {"mask": mask, "skeleton": skeleton, "seed": seed}


def trace_all_neurons(img, centroids, nuclei_labeled, channel=0,
                      threshold_pct=0.10, min_reach_ratio=2.0,
                      min_mask_size=100000):
    """
    Trace all neurons, filtering out non-neurons.

    A detection is kept only if:
    1. The tubulin mask is at least min_mask_size pixels
    2. The tubulin mask extends at least min_reach_ratio × nucleus_radius

    Parameters
    ----------
    img : array
        RGB image.
    centroids : list of (y, x)
        Seed points from nucleus_detector.
    nuclei_labeled : 2D array
        Labelled nuclei mask from nucleus_detector.
    channel : int
        Channel to trace on.
    threshold_pct : float
        Intensity threshold as fraction of seed value.
    min_reach_ratio : float
        Minimum geodesic reach as multiple of nucleus radius.
    min_mask_size : int
        Minimum neuron mask area in pixels. Default 100000.

    Returns
    -------
    traces : list of trace dicts (only neurons that passed filter)
    kept_centroids : list of (y, x) centroids that passed filter
    """
    traces = []
    kept_centroids = []

    for i, seed in enumerate(centroids):
        trace = trace_neuron(img, seed, channel, threshold_pct)
        mask_size = trace["mask"].sum()

        # Filter 1: minimum mask size
        if mask_size < min_mask_size:
            print(f"  Neuron {i+1}: FILTERED (too small) — "
                  f"mask={mask_size} px, min={min_mask_size} px")
            continue

        # Get this nucleus mask and radius
        label_at_seed = nuclei_labeled[seed[0], seed[1]]
        if label_at_seed > 0:
            nuc_mask = nuclei_labeled == label_at_seed
        else:
            h, w = nuclei_labeled.shape
            yy, xx = np.ogrid[:h, :w]
            nuc_mask = ((yy - seed[0])**2 + (xx - seed[1])**2) < 15**2

        nuc_r = nucleus_radius(nuc_mask)
        min_reach = nuc_r * min_reach_ratio

        # Filter 2: geodesic reach
        if trace["mask"].any():
            reach = geodesic_max_distance(trace["mask"], seed)
        else:
            reach = 0.0

        if reach >= min_reach:
            traces.append(trace)
            kept_centroids.append(seed)
            print(f"  Neuron {i+1}: KEPT — mask={mask_size} px, "
                  f"reach={reach:.0f} px, threshold={min_reach:.0f} px")
        else:
            print(f"  Neuron {i+1}: FILTERED (low reach) — "
                  f"mask={mask_size} px, reach={reach:.0f} px, "
                  f"threshold={min_reach:.0f} px")

    print(f"\n  Kept {len(traces)}/{len(centroids)} neurons")
    return traces, kept_centroids


def show_traces(img, traces, nuclei_labeled=None, save_path=None):
    """
    Visualise all traced neurons, optionally with nuclei overlay.
    """
    gray = img[:, :, 0].astype(np.float64)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Original with seeds and nuclei
    axes[0].imshow(img)
    if nuclei_labeled is not None:
        axes[0].contour(nuclei_labeled > 0, colors='cyan', linewidths=0.8)
    for i, t in enumerate(traces):
        y, x = t["seed"]
        axes[0].plot(x, y, 'c+', markersize=10, markeredgewidth=1.5)
        axes[0].text(x + 5, y - 5, str(i + 1), color='cyan', fontsize=7)
    axes[0].set_title(f"Seeds + nuclei ({len(traces)} neurons)")
    axes[0].axis("off")

    # All masks with nuclei
    axes[1].imshow(gray, cmap="gray")
    if nuclei_labeled is not None:
        axes[1].contour(nuclei_labeled > 0, colors='cyan', linewidths=0.8)
    colours = plt.cm.Set1(np.linspace(0, 1, max(len(traces), 1)))
    for i, t in enumerate(traces):
        if t["mask"].any():
            axes[1].contour(t["mask"], colors=[colours[i % len(colours)]], linewidths=0.5)
    axes[1].set_title("Masks + nuclei")
    axes[1].axis("off")

    # Skeletons with nuclei
    skel_overlay = np.stack([gray / max(gray.max(), 1)] * 3, axis=-1)
    if nuclei_labeled is not None:
        from skimage.segmentation import find_boundaries
        nuclei_border = find_boundaries(nuclei_labeled, mode='outer')
        skel_overlay[nuclei_border, 0] = 0
        skel_overlay[nuclei_border, 1] = 1
        skel_overlay[nuclei_border, 2] = 1
    for i, t in enumerate(traces):
        c = colours[i % len(colours)]
        skel_overlay[t["skeleton"], 0] = c[0]
        skel_overlay[t["skeleton"], 1] = c[1]
        skel_overlay[t["skeleton"], 2] = c[2]
    axes[2].imshow(skel_overlay)
    axes[2].set_title("Skeletons + nuclei")
    axes[2].axis("off")

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    plt.show()