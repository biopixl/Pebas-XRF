#!/usr/bin/env python3
"""
Process Itrax optical images: brighten and add position scale
"""

import os
from pathlib import Path
from PIL import Image, ImageEnhance, ImageDraw, ImageFont

# Configuration
BASE_PATH = Path(__file__).parent.parent
DATA_PATH = BASE_PATH / "TAM-SC-IsaacA"
OUTPUT_PATH = BASE_PATH / "output" / "figures" / "optical_review"
OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

BRIGHTNESS_FACTOR = 1.8
CONTRAST_FACTOR = 1.3

def get_scan_positions(optical_path):
    """Read start/stop positions from document.txt"""
    doc_path = optical_path.parent / "document.txt"

    if not doc_path.exists():
        return None, None

    start, stop = None, None

    with open(doc_path, 'r', encoding='latin-1') as f:
        for line in f:
            if "Start coordinate" in line:
                parts = line.strip().split('\t')
                try:
                    start = float(parts[1])
                    stop = float(parts[3])
                except (IndexError, ValueError):
                    pass

    return start, stop

def process_optical_image(optical_path, brightness=1.8, contrast=1.3):
    """Brighten image and save as PNG"""

    section_name = optical_path.parent.name
    print(f"Processing: {section_name}")

    try:
        img = Image.open(optical_path)
    except Exception as e:
        print(f"  Error reading image: {e}")
        return None

    # Convert to RGB if needed
    if img.mode != 'RGB':
        img = img.convert('RGB')

    # Brighten
    enhancer = ImageEnhance.Brightness(img)
    img = enhancer.enhance(brightness)

    # Increase contrast
    enhancer = ImageEnhance.Contrast(img)
    img = enhancer.enhance(contrast)

    # Get positions
    start_mm, stop_mm = get_scan_positions(optical_path)

    # Add scale markers
    draw = ImageDraw.Draw(img)
    width, height = img.size

    # Try to load a font, fall back to default
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 24)
        font_small = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 18)
    except:
        font = ImageFont.load_default()
        font_small = font

    # Add title
    title = f"{section_name}"
    if start_mm and stop_mm:
        title += f" | {start_mm:.0f} - {stop_mm:.0f} mm"
    draw.text((10, 10), title, fill="red", font=font)

    # Add scale markers every 100mm if we have position info
    if start_mm and stop_mm:
        mm_range = stop_mm - start_mm
        mm_per_px = mm_range / width

        for pos_mm in range(0, int(mm_range) + 1, 100):
            px = int(pos_mm / mm_per_px)
            if px < width:
                # Draw tick mark
                draw.line([(px, height - 30), (px, height)], fill="red", width=2)
                # Draw label
                label = f"{int(start_mm + pos_mm)}"
                draw.text((px - 15, height - 50), label, fill="red", font=font_small)

    # Save
    out_file = OUTPUT_PATH / f"optical_{section_name}.png"

    # Resize if too large (max 2000px width for easier viewing)
    if width > 2000:
        ratio = 2000 / width
        new_size = (2000, int(height * ratio))
        img = img.resize(new_size, Image.Resampling.LANCZOS)

    img.save(out_file, "PNG", quality=95)
    print(f"  Saved: {out_file.name}")

    return out_file

def main():
    # Find all optical images (excluding 'copied' directories)
    optical_files = []
    for f in DATA_PATH.rglob("optical.tif"):
        if "/copied/" not in str(f):
            optical_files.append(f)

    print(f"Found {len(optical_files)} optical images\n")

    # Process each
    processed = 0
    for optical_path in sorted(optical_files):
        result = process_optical_image(optical_path, BRIGHTNESS_FACTOR, CONTRAST_FACTOR)
        if result:
            processed += 1

    print(f"\n{processed} images processed")
    print(f"Output directory: {OUTPUT_PATH}")

    # Create summary
    print("\n=== SECTIONS PROCESSED ===")
    for f in sorted(OUTPUT_PATH.glob("optical_*.png")):
        print(f"  {f.name}")

if __name__ == "__main__":
    main()
