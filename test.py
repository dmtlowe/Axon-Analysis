from skimage import io
from nucleus_detector import detect_nuclei
from neuron_trace import trace_all_neurons, show_traces
from measure_axon import measure_axon, show_axon

img_path = r"C:\Uni\Projects\Axon_Analysis\data\PAVNCD1_1-5-_8bit_panel_panels\PAVNCD1_1-5-_8bit_panel_R3_C5.tif"

img = io.imread(img_path)
centroids, labeled = detect_nuclei(img)

# Now returns (traces, kept_centroids) — filtered list
traces, kept_centroids = trace_all_neurons(
    img, centroids, labeled, 
    channel=0, 
    threshold_pct=0.10,
    min_reach_ratio=6.0   # tubulin must extend 4× nucleus radius
)


nuc_mask = labeled == labeled[kept_centroids[0][0], kept_centroids[0][1]]
result = measure_axon(img, traces[0]["mask"], kept_centroids[0], nuc_mask, channel=0)
print(f"Axon length: {result['length_px']:.1f} px")

show_axon(img, result, nucleus_mask=nuc_mask)


# img = io.imread(img_path)
# centroids, labeled = detect_nuclei(img)
# traces = trace_all_neurons(img, centroids, channel=0, threshold_pct=0.10)

# # Measure one neuron
# i = 0
# nuc_mask = labeled == labeled[centroids[i][0], centroids[i][1]]
# result = measure_axon(img, traces[i]["mask"], centroids[i], nuc_mask, channel=0)
# print(f"Axon length: {result['length_px']:.1f} px")

# show_axon(img, result, nucleus_mask=nuc_mask, save_path="axon.png")

# img = io.imread(img_path)
# centroids, labeled = detect_nuclei(img)

# # Now returns (traces, kept_centroids) — filtered list
# traces, kept_centroids = trace_all_neurons(
#     img, centroids, labeled, 
#     channel=0, 
#     threshold_pct=0.10,
#     min_reach_ratio=4.0   # tubulin must extend 4× nucleus radius
# )

# show_traces(img, traces, nuclei_labeled=labeled, save_path="traces.png")