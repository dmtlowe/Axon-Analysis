import sys
from skimage import io
import matplotlib.pyplot as plt
import matplotlib.patches as patches

image_path = r"C:\Uni\Projects\Axon_Analysis\data\PAVNCD1_1-5-_8bit_panel_panels\PAVNCD1_1-5-_8bit_panel_R2_C4.tif"
inset_px = 350

def show_border(image_path, inset_px):
    # Load image
    img = io.imread(image_path)

    h, w = img.shape[:2]

    if inset_px * 2 >= min(h, w):
        raise ValueError("Inset too large for image dimensions.")

    fig, ax = plt.subplots()
    ax.imshow(img)
    ax.axis("off")

    # Rectangle: (x, y), width, height
    rect = patches.Rectangle(
        (inset_px, inset_px),
        w - 2 * inset_px,
        h - 2 * inset_px,
        linewidth=2,
        edgecolor='red',
        facecolor='none'
    )

    ax.add_patch(rect)
    plt.show()


if __name__ == "__main__":

    show_border(image_path, inset_px)