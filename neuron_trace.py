"""
Neuron tracer - click the soma, get a trace.

Usage:
    python neuron_trace.py path/to/image.tif
    
A window opens. Click the soma. Close the window. Get results.
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from skimage import io, morphology, measure, filters
from collections import deque
from pathlib import Path


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


def trace_from_click(img_path, channel=0, threshold_pct=0.10):
    """Open image, let user click soma, run trace."""
    img = io.imread(str(img_path))
    
    # Extract channel
    if img.ndim == 3:
        gray = img[:, :, channel].astype(np.float64)
    else:
        gray = img.astype(np.float64)

    # Show image, collect click
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.imshow(img)
    ax.set_title("Click the soma, then close this window")
    clicks = []

    def on_click(event):
        if event.inaxes:
            clicks.append((int(event.ydata), int(event.xdata)))
            ax.plot(event.xdata, event.ydata, 'c+', markersize=15, markeredgewidth=2)
            fig.canvas.draw()
            print(f"Clicked: y={int(event.ydata)}, x={int(event.xdata)}")

    fig.canvas.mpl_connect('button_press_event', on_click)
    plt.show(block=True)

    if not clicks:
        print("No click detected.")
        return

    seed = clicks[0]
    print(f"\nSeed: {seed}")
    print(f"Intensity at seed: {gray[seed[0], seed[1]]:.1f}")

    # Region grow
    mask = region_grow(gray, seed, threshold_pct)
    
    # Clean up
    mask = morphology.remove_small_objects(mask, min_size=50)
    mask = morphology.remove_small_holes(mask, area_threshold=100)

    # Keep only component containing seed
    labeled = measure.label(mask)
    seed_label = labeled[seed[0], seed[1]]
    if seed_label > 0:
        mask = labeled == seed_label

    # Skeletonize
    skeleton = morphology.skeletonize(mask)

    # Show results
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    axes[0].imshow(img)
    axes[0].plot(seed[1], seed[0], 'c+', markersize=15, markeredgewidth=2)
    axes[0].set_title("Original + seed")
    axes[0].axis("off")

    axes[1].imshow(gray, cmap="gray")
    axes[1].contour(mask, colors='lime', linewidths=0.5)
    axes[1].set_title(f"Mask ({mask.sum()} px)")
    axes[1].axis("off")

    axes[2].imshow(gray, cmap="gray")
    skel_y, skel_x = np.where(skeleton)
    axes[2].scatter(skel_x, skel_y, c='lime', s=0.3)
    axes[2].set_title(f"Skeleton ({skeleton.sum()} px)")
    axes[2].axis("off")

    plt.tight_layout()
    
    # Save
    out_path = Path(img_path).with_suffix('.trace.png')
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python neuron_trace.py image.tif")
        sys.exit(1)
    
    trace_from_click(sys.argv[1])