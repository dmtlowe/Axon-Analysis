# Axon Analysis Pipeline — Documentation

## Overview

A Python pipeline for automated measurement of axon length and microtubule disorganisation from fluorescence microscopy images. Built to replace manual ImageJ/Fiji measurements.

**Input:** Multi-panel TIF images with up to 3 channels (DAPI/blue, HRP/green, Tubulin/red).

**Output:** Axon length measurements in pixels (convertible to µm via scaling factor).

---

## Pipeline Flow

```
panel_splitter.py     → Split multi-panel TIF into individual panels
        ↓
show_channels.py      → (Optional) Inspect R/G/B channels
        ↓
nucleus_detector.py   → Detect neuronal nuclei from DAPI, filtered by HRP
        ↓
neuron_trace.py       → Region-grow tubulin masks from nuclei, filter non-neurons
        ↓
measure_axon.py       → Trace intensity-guided axon path, measure length
```

---

## Module: panel_splitter.py

Splits a single large TIF containing multiple panels into individual images.

### `split_panels(filepath, rows=5, cols=6, output_dir=None)`

Splits a TIF image into a grid of panels and saves each one as a separate TIF.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `filepath` | str | required | Path to the multi-panel TIF file |
| `rows` | int | 5 | Number of rows in the panel grid |
| `cols` | int | 6 | Number of columns in the panel grid |
| `output_dir` | str | None | Output directory. If None, creates a folder named `{filename}_panels` next to the input |

**Returns:** `list` of numpy arrays, one per panel.

**Output files:** Individual TIF files named `{stem}_panel{00-29}.tif` in the output directory.

**Notes:**
- Panel size is calculated as `image_height // rows` by `image_width // cols`
- Leftover pixels at the edges are trimmed
- Handles RGB images (H, W, 3) — preserves all channels in output

**Example:**
```python
from panel_splitter import split_panels
panels = split_panels("experiment.tif")
# Also works from command line:
# python panel_splitter.py experiment.tif
```

---

## Module: show_channels.py

Utility to visualise individual RGB channels.

### `show_channels(img_path)`

Opens an image and displays Composite, Red, Green, and Blue channels side by side.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `img_path` | str | required | Path to the image file |

**Output files:** Saves `{filename}.channels.png` next to the input file.

**Example:**
```python
# Command line:
# python show_channels.py panel_00.tif
```

---

## Module: nucleus_detector.py

Detects neuronal nuclei from DAPI staining, filtered by HRP colocalisation.

### `detect_nuclei(img, min_size=100, blur_sigma=5, min_peak_distance=50, require_hrp=True, hrp_channel=1, hrp_ring_width=30, hrp_threshold_factor=2.0)`

Detects nuclei from the blue channel of an RGB image. Optionally filters to keep only nuclei with nearby HRP (green) signal — removes non-neuronal DAPI signals.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `img` | numpy array | required | RGB image array (H, W, 3) |
| `min_size` | int | 100 | Minimum nucleus area in pixels. Removes small debris |
| `blur_sigma` | float | 5 | Gaussian blur sigma applied to DAPI channel before thresholding. Smooths internal chromatin structure so each nucleus is detected as one object, not split into fragments |
| `min_peak_distance` | int | 50 | Minimum distance (px) between watershed seeds. Prevents one nucleus being split into multiple detections. Increase if nuclei are still being split |
| `require_hrp` | bool | True | If True, only keep nuclei that have green (HRP) signal in a surrounding patch. Filters out non-neuronal DAPI signals |
| `hrp_channel` | int | 1 | Which channel contains HRP staining (0=R, 1=G, 2=B) |
| `hrp_ring_width` | int | 30 | Half-width of square patch around each centroid used to check for HRP signal (pixels) |
| `hrp_threshold_factor` | float | 2.0 | Nucleus is kept if mean green intensity in patch > factor × median green of whole image. Increase to be stricter about requiring HRP |

**Returns:**
- `centroids` — list of (y, x) tuples, one per detected nucleus
- `labeled` — 2D numpy array, same size as image, where each nucleus has a unique integer label (0 = background)

**Algorithm:**
1. Extract blue channel, apply Gaussian blur
2. Otsu threshold to create binary mask
3. Remove small objects and fill holes
4. Watershed segmentation to separate touching nuclei
5. For each nucleus, check mean green intensity in surrounding patch
6. Discard nuclei without sufficient HRP signal

**Tuning guide:**
- Nuclei being split into fragments → increase `blur_sigma` or `min_peak_distance`
- Missing dim nuclei → lower the Otsu threshold won't help (it's automatic); check if HRP filter is removing them (`require_hrp=False` to test)
- Rogue DAPI signals passing through → increase `hrp_threshold_factor` to 3.0 or 4.0
- HRP check too slow → it shouldn't be (uses simple array slicing), but reduce `hrp_ring_width` if needed

### `show_nuclei(img_path)`

Standalone visualisation. Runs `detect_nuclei` and displays three panels: DAPI channel, labelled nuclei, and centroids on composite.

**Example:**
```python
from nucleus_detector import detect_nuclei
from skimage import io

img = io.imread("panel_00.tif")
centroids, labeled = detect_nuclei(img)
print(f"Found {len(centroids)} nuclei")
```

---

## Module: neuron_trace.py

Region-grows tubulin masks from nucleus seed points. Filters out non-neuronal detections based on tubulin reach and mask size.

### `region_grow(img_gray, seed, threshold_pct=0.10)`

Low-level function. Flood-fill from a seed point, including all connected pixels above an intensity threshold.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `img_gray` | 2D numpy array | required | Single-channel grayscale image |
| `seed` | tuple (y, x) | required | Starting pixel coordinates |
| `threshold_pct` | float | 0.10 | Minimum intensity as fraction of seed pixel intensity. 0.10 means include pixels above 10% of the seed value |

**Returns:** 2D boolean mask.

**Notes:**
- Uses 8-connected flood fill (includes diagonal neighbours)
- If seed pixel intensity is low, the absolute threshold will be very low — may include noise

### `trace_neuron(img, seed, channel=0, threshold_pct=0.10)`

Traces a single neuron from a seed point by region growing on the specified channel.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `img` | numpy array | required | RGB image (H, W, 3) |
| `seed` | tuple (y, x) | required | Seed point, typically a nucleus centroid |
| `channel` | int | 0 | Channel to trace on. 0=red (tubulin), 1=green (HRP), 2=blue (DAPI) |
| `threshold_pct` | float | 0.10 | Intensity threshold for region growing |

**Returns:** dict with keys:
- `mask` — 2D boolean array, the neuron mask
- `skeleton` — 2D boolean array, skeletonised version of mask
- `seed` — the input seed point

**Algorithm:**
1. Extract specified channel
2. Region grow from seed
3. Remove small objects (<50 px) and fill small holes (<100 px)
4. Keep only the connected component containing the seed
5. Skeletonise the mask

### `trace_all_neurons(img, centroids, nuclei_labeled, channel=0, threshold_pct=0.10, min_reach_ratio=2.0, min_mask_size=100000)`

Traces all neurons and filters out non-neuronal detections.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `img` | numpy array | required | RGB image (H, W, 3) |
| `centroids` | list of (y, x) | required | Seed points from `detect_nuclei` |
| `nuclei_labeled` | 2D array | required | Labelled nuclei mask from `detect_nuclei` |
| `channel` | int | 0 | Channel to trace on |
| `threshold_pct` | float | 0.10 | Region grow intensity threshold |
| `min_reach_ratio` | float | 2.0 | Minimum geodesic reach as multiple of nucleus radius. A detection is kept only if tubulin extends at least this many nucleus-radii from the centroid. Increase to be stricter |
| `min_mask_size` | int | 100000 | Minimum neuron mask area in pixels. Hard floor — masks smaller than this are discarded regardless of reach |

**Returns:**
- `traces` — list of trace dicts (only neurons that passed both filters)
- `kept_centroids` — list of (y, x) centroids corresponding to kept traces

**Filtering logic (applied in order):**
1. **Mask size filter:** if `mask.sum() < min_mask_size` → discard (fast check)
2. **Geodesic reach filter:** compute max geodesic distance from centroid through the mask. If less than `min_reach_ratio × nucleus_radius` → discard. This removes DAPI+tubulin blobs that don't form neurites

**Tuning guide:**
- Small real neurons being filtered → lower `min_mask_size` or `min_reach_ratio`
- Blobs still passing → increase `min_mask_size` or `min_reach_ratio`
- Console output shows exactly why each neuron was kept or filtered

### `show_traces(img, traces, nuclei_labeled=None, save_path=None)`

Visualises all traced neurons with optional nuclei overlay.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `img` | numpy array | required | RGB image |
| `traces` | list of dicts | required | Output from `trace_all_neurons` |
| `nuclei_labeled` | 2D array | None | Labelled nuclei mask. If provided, nuclei outlines shown in cyan |
| `save_path` | str | None | Path to save figure |

**Displays three panels:** seeds + nuclei on composite, masks + nuclei on red channel, skeletons + nuclei borders.

**Example:**
```python
from skimage import io
from nucleus_detector import detect_nuclei
from neuron_trace import trace_all_neurons, show_traces

img = io.imread("panel_00.tif")
centroids, labeled = detect_nuclei(img)
traces, kept_centroids = trace_all_neurons(
    img, centroids, labeled,
    channel=0, threshold_pct=0.10,
    min_reach_ratio=2.0, min_mask_size=100000
)
show_traces(img, traces, nuclei_labeled=labeled, save_path="traces.png")
```

---

## Module: measure_axon.py

Measures axon length by tracing an intensity-guided path from the nucleus boundary to the axon tip.

### `geodesic_distance(mask, seed)`

Computes geodesic distance from a seed point to every reachable pixel in a mask.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mask` | 2D bool array | required | Binary mask defining traversable pixels |
| `seed` | tuple (y, x) | required | Starting point. If outside mask, snaps to nearest mask pixel |

**Returns:** 2D float array. 0 at seed, increasing with distance, `inf` for unreachable pixels. Cardinal steps cost 1.0, diagonal steps cost 1.414.

### `find_start_point(nucleus_mask, neuron_mask, target)`

Finds the best starting point on the nucleus boundary for axon measurement.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `nucleus_mask` | 2D bool array | required | Mask for one nucleus |
| `neuron_mask` | 2D bool array | required | Tubulin mask for the neuron |
| `target` | tuple (y, x) | required | The axon tip point. Start point is chosen as the nucleus boundary pixel closest to this target |

**Returns:** (y, x) tuple — the start point on the nucleus boundary.

### `trace_intensity_path(gray, start, end, neuron_mask, intensity_weight=2.0)`

Traces a path from start to end following bright tubulin signal, constrained to the neuron mask. Uses Dijkstra's algorithm where brighter pixels are cheaper to traverse.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gray` | 2D float array | required | Tubulin channel intensity |
| `start` | tuple (y, x) | required | Path start (nucleus boundary) |
| `end` | tuple (y, x) | required | Path end (axon tip) |
| `neuron_mask` | 2D bool array | required | Constrains path to within the mask |
| `intensity_weight` | float | 2.0 | How strongly intensity pulls the path. Higher = sticks more tightly to bright signal. Lower = more direct route |

**Returns:**
- `path` — list of (y, x) tuples from start to end
- `length_px` — float, total path length in pixels (diagonal steps count as √2)

**Cost function:** `spatial_cost / (intensity ^ weight)` — bright pixels have low cost, dark pixels have high cost. The path naturally follows the intensity ridge (centre of the axon).

**Tuning guide:**
- Path wanders off the bright signal → increase `intensity_weight` to 3.0 or 4.0
- Path is too rigid / misses curves → decrease to 1.5
- Path goes through noise → the neuron mask should constrain it, but tightening `threshold_pct` in the tracing step will help

### `measure_axon(img, neuron_mask, nucleus_centroid, nucleus_mask, channel=0, intensity_weight=2.0)`

Full axon measurement pipeline for one neuron.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `img` | numpy array | required | RGB image |
| `neuron_mask` | 2D bool array | required | Tubulin mask from `trace_neuron` |
| `nucleus_centroid` | tuple (y, x) | required | Nucleus centre from `detect_nuclei` |
| `nucleus_mask` | 2D bool array | required | Mask for this specific nucleus |
| `channel` | int | 0 | Tubulin channel (0=red) |
| `intensity_weight` | float | 2.0 | How strongly path follows bright signal |

**Returns:** dict with keys:
- `path` — list of (y, x) tuples, the axon path from nucleus boundary to tip
- `length_px` — float, axon length in pixels
- `start_point` — (y, x) where path starts on nucleus boundary
- `tip_point` — (y, x) the axon tip (geodesic furthest point)
- `distance_map` — 2D float array, full geodesic distance map

**Algorithm:**
1. Compute geodesic distance from centroid through neuron mask → find axon tip (furthest point)
2. Find start point on nucleus boundary facing the tip
3. Trace intensity-weighted shortest path from start to tip (Dijkstra)
4. Path length = axon length

**Converting pixels to µm:**
```python
scaling_factor = 29.21  # px/µm — from your image metadata (X Resolution)
length_um = result["length_px"] / scaling_factor
```

### `measure_all_axons(img, traces, centroids, nuclei_labeled, channel=0, intensity_weight=2.0)`

Measures axons for all traced neurons.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `img` | numpy array | required | RGB image |
| `traces` | list of dicts | required | From `trace_all_neurons` |
| `centroids` | list of (y, x) | required | Kept centroids from `trace_all_neurons` |
| `nuclei_labeled` | 2D array | required | Labelled nuclei mask |
| `channel` | int | 0 | Tubulin channel |
| `intensity_weight` | float | 2.0 | Path intensity weighting |

**Returns:** list of axon result dicts (same structure as `measure_axon`).

### `show_axon(img, axon_result, nucleus_mask=None, save_path=None)`

Visualises the measured axon path.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `img` | numpy array | required | RGB image |
| `axon_result` | dict | required | From `measure_axon` |
| `nucleus_mask` | 2D bool array | None | Shown as cyan outline if provided |
| `save_path` | str | None | Path to save figure |

**Displays two panels:** axon path (green line) on composite with start/tip markers, and path on tubulin channel.

**Example:**
```python
from skimage import io
from nucleus_detector import detect_nuclei
from neuron_trace import trace_all_neurons
from measure_axon import measure_all_axons, show_axon

img = io.imread("panel_00.tif")
centroids, labeled = detect_nuclei(img)
traces, kept_centroids = trace_all_neurons(img, centroids, labeled, channel=0)
results = measure_all_axons(img, traces, kept_centroids, labeled, channel=0)

# Show first neuron
nuc_mask = labeled == labeled[kept_centroids[0][0], kept_centroids[0][1]]
show_axon(img, results[0], nucleus_mask=nuc_mask, save_path="axon.png")

# Print all lengths
for i, r in enumerate(results):
    print(f"Neuron {i+1}: {r['length_px']:.1f} px")
```

---

## Module: furthest_tubulin_finder.py

Earlier iteration — geodesic distance and path tracing. Superseded by `measure_axon.py` for measurement, but `geodesic_distance` and `show_geodesic` remain useful for debugging.

---

## Full Pipeline Example

```python
from skimage import io
from panel_splitter import split_panels
from nucleus_detector import detect_nuclei
from neuron_trace import trace_all_neurons, show_traces
from measure_axon import measure_all_axons, show_axon

# Step 1: Split panels
panels = split_panels("experiment.tif")

# Step 2: Process one panel
img = io.imread("experiment_panels/experiment_panel00.tif")

# Step 3: Detect nuclei
centroids, labeled = detect_nuclei(img)

# Step 4: Trace neurons (with filtering)
traces, kept_centroids = trace_all_neurons(
    img, centroids, labeled,
    channel=0,
    threshold_pct=0.10,
    min_reach_ratio=2.0,
    min_mask_size=100000
)

# Step 5: Measure axons
results = measure_all_axons(img, traces, kept_centroids, labeled, channel=0)

# Step 6: Output
scaling_factor = 29.21  # px/µm
for i, r in enumerate(results):
    length_um = r["length_px"] / scaling_factor
    print(f"Neuron {i+1}: {r['length_px']:.1f} px = {length_um:.1f} µm")
```

---

## Dependencies

```
pip install numpy matplotlib scikit-image scipy tifffile
```

## Scaling Factor

The X Resolution in your TIF metadata is (29210000, 1000000) = 29.21 px/µm. To convert pixel measurements to micrometres, divide by this value. This factor depends on your microscope objective and camera settings (typically: 100x oil objective, optovar 1, binning 1x1).
