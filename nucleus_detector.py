"""
Nucleus detector - find nuclei from DAPI (blue channel).

Usage:
    python nucleus_detector.py image.tif

Returns list of (y, x) centroids and labelled mask.
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from skimage import io, morphology, measure, filters
from scipy import ndimage
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from pathlib import Path


def detect_nuclei(img, min_size=100, blur_sigma=5, min_peak_distance=50):
    """
    Detect nuclei from the blue channel of an RGB image.
    
    Parameters
    ----------
    img : array
        RGB image.
    min_size : int
        Minimum nucleus area in pixels (removes debris).
    blur_sigma : float
        Gaussian blur sigma. Smooths internal chromatin structure
        so each nucleus becomes one blob.
    min_peak_distance : int
        Minimum distance between watershed seeds. Prevents splitting
        one nucleus into multiple detections.
    
    Returns
    -------
    centroids : list of (y, x) tuples
    labeled : 2D array, each nucleus has a unique integer label
    """
    blue = img[:, :, 2].astype(np.float64)
    
    # Blur to merge internal structure
    blurred = filters.gaussian(blue, sigma=blur_sigma)
    
    # Otsu threshold on blurred
    thresh = filters.threshold_otsu(blurred)
    mask = blurred > thresh
    
    # Clean up
    mask = morphology.remove_small_objects(mask, min_size=min_size)
    mask = morphology.remove_small_holes(mask, area_threshold=500)
    
    # Watershed with increased min_distance
    distance = ndimage.distance_transform_edt(mask)
    coords = peak_local_max(distance, min_distance=min_peak_distance, labels=mask)
    peak_mask = np.zeros(distance.shape, dtype=bool)
    peak_mask[tuple(coords.T)] = True
    markers = measure.label(peak_mask)
    labeled = watershed(-distance, markers, mask=mask)
    
    # Filter by min size
    centroids = []
    filtered_label = np.zeros_like(labeled)
    new_id = 1
    for region in measure.regionprops(labeled):
        if region.area >= min_size:
            filtered_label[labeled == region.label] = new_id
            centroids.append((int(region.centroid[0]), int(region.centroid[1])))
            new_id += 1
    
    return centroids, filtered_label


def show_nuclei(img_path):
    """Detect and display nuclei."""
    img = io.imread(str(img_path))
    centroids, labeled = detect_nuclei(img)
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    axes[0].imshow(img[:, :, 2], cmap="gray")
    axes[0].set_title("Blue channel (DAPI)")
    axes[0].axis("off")
    
    axes[1].imshow(labeled, cmap="nipy_spectral")
    axes[1].set_title(f"Labelled nuclei ({labeled.max()} found)")
    axes[1].axis("off")
    
    axes[2].imshow(img)
    for i, (cy, cx) in enumerate(centroids):
        axes[2].plot(cx, cy, 'c+', markersize=10, markeredgewidth=1.5)
        axes[2].text(cx + 5, cy - 5, str(i+1), color='cyan', fontsize=7)
    axes[2].set_title(f"Centroids ({len(centroids)} nuclei)")
    axes[2].axis("off")
    
    plt.tight_layout()
    out = Path(img_path).with_suffix(".nuclei.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Found {len(centroids)} nuclei")
    print(f"Saved: {out}")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python nucleus_detector.py image.tif")
        sys.exit(1)
    show_nuclei(sys.argv[1])