import os
import pandas as pd
from pathlib import Path
from skimage import io 
from nucleus_detector import detect_nuclei
from neuron_trace import trace_all_neurons, show_traces
from measure_axon import measure_axon, show_axon

# Test script for batch processing multiple images in a directory 

img_folder_path = r"C:\Uni\Projects\Axon_Analysis\data\PAVNCD1_1-5-_8bit_panel_panels"
save_folder_path = r"C:\Uni\Projects\Axon_Analysis\data\results\PAVNCD1_1-5-_8bit_panel"


os.makedirs(save_folder_path, exist_ok=True)

files = sorted([f for f in os.listdir(img_folder_path) if f.endswith(".tif")])
print(f"Found {len(files)} tif files")

rows = []

for file_idx, filename in enumerate(files, start=1):
    print(f"Processing file {file_idx}/{len(files)}: {filename}")

    img_path = os.path.join(img_folder_path, filename)
    img = io.imread(img_path)

    centroids, labeled = detect_nuclei(img)

    traces, kept_centroids = trace_all_neurons(
        img, centroids, labeled,
        channel=0,
        threshold_pct=0.10,
        min_reach_ratio=4.0
    )

    for neuron_idx in range(len(traces)):
        nuc_mask = labeled == labeled[
            kept_centroids[neuron_idx][0],
            kept_centroids[neuron_idx][1]
        ]

        result = measure_axon(
            img,
            traces[neuron_idx]["mask"],
            kept_centroids[neuron_idx],
            nuc_mask,
            channel=0
        )

        print(f"  Neuron {neuron_idx}: Axon length {result['length_px']:.1f} px")

        save_name = f"{Path(filename).stem}_neuron_{neuron_idx}.png"

        show_axon(
            img,
            result,
            nucleus_mask=nuc_mask,
            save_path=os.path.join(save_folder_path, save_name)
        )

        rows.append({
            "filename": filename,
            "neuron_index": neuron_idx,
            "axon_length_px": result["length_px"]
        })

axon_results = pd.DataFrame(rows)
axon_results.to_csv(
    os.path.join(save_folder_path, "axon_measurements.csv"),
    index=False
)

print(f"Saved {len(axon_results)} neuron measurements")