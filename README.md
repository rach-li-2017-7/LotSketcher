# Lot Sketcher

A Python program that generates lot sketches within a given polygon, optimizing for maximum lot numbers and minimum road length.

## Requirements

- Python 3.7+
- Required packages are listed in `requirements.txt`

## Installation

1. Clone this repository
2. Install the required packages:
```bash
pip install -r requirements.txt
```

## Usage

Run the program:
```bash
python lot_sketcher.py
```

The program will:
1. Create a sample polygon (currently using New York City coordinates as an example)
2. Ask for input parameters:
   - Minimum lot size (in square meters)
   - Minimum front line length (in meters)
   - Road width (in meters)
3. Generate a lot sketch
4. Save the result as a GeoJSON file named `lot_sketch.geojson`

## Output

The program generates a GeoJSON file containing:
- Lot polygons
- Road networks
- All coordinates are in WGS84 (latitude/longitude) format

## Note

This is a work in progress. The current version includes the basic structure and coordinate transformation, but the actual lot generation algorithm needs to be implemented.