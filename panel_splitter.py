"""
Split a multi-panel TIF into individual panels (5 rows × 6 cols).
"""

import numpy as np
import tifffile
from pathlib import Path


def split_panels(filepath, rows=5, cols=6, output_dir=None):
    """
    Split a TIF image into a grid of panels and save each one.
    
    Parameters
    ----------
    filepath : str
        Path to the TIF file.
    rows, cols : int
        Grid layout (default: 5×6 = 30 panels).
    output_dir : str, optional
        Where to save. Defaults to a folder next to the input file.
    """
    filepath = Path(filepath)
    img = tifffile.imread(str(filepath))
    
    h, w = img.shape[0], img.shape[1]
    ph = h // rows
    pw = w // cols
    
    if output_dir is None:
        output_dir = filepath.parent / f"{filepath.stem}_panels"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Image: {img.shape} → {rows}×{cols} grid, panel size {ph}×{pw}")
    
    panels = []
    for r in range(rows):
        for c in range(cols):
            idx = r * cols + c
            panel = img[r*ph:(r+1)*ph, c*pw:(c+1)*pw]
            panels.append(panel)
            
            out_path = output_dir / f"{filepath.stem}_panel{idx:02d}.tif"
            tifffile.imwrite(str(out_path), panel)
    
    print(f"Saved {len(panels)} panels to {output_dir}")
    return panels


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        split_panels(sys.argv[1])
    else:
        print("Usage: python panel_splitter.py image.tif")