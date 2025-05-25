import numpy as np
from shapely.geometry import Polygon, Point, LineString, shape
import geojson
from pyproj import Transformer
import json
from typing import Tuple, List
import folium
import os

def square_feet_to_square_meters(sq_ft: float) -> float:
    """Convert square feet to square meters."""
    return sq_ft * 0.09290304

def feet_to_meters(ft: float) -> float:
    """Convert feet to meters."""
    return ft * 0.3048

def load_polygon_from_geojson(file_path: str) -> Tuple[Polygon, Point]:
    """
    Load polygon and entrance point from a GeoJSON file.
    The file should contain a FeatureCollection with:
    - A polygon feature for the boundary
    - A point feature for the entrance
    """
    with open(file_path, 'r') as f:
        data = json.load(f)
    
    polygon = None
    entrance_point = None
    
    for feature in data['features']:
        if feature['geometry']['type'] == 'Polygon':
            polygon = shape(feature['geometry'])
        elif feature['geometry']['type'] == 'Point':
            entrance_point = Point(feature['geometry']['coordinates'])
    
    if polygon is None or entrance_point is None:
        raise ValueError("GeoJSON file must contain both a polygon and a point feature")
    
    return polygon, entrance_point

def create_sample_polygon() -> Tuple[Polygon, Point]:
    """
    Creates a sample polygon and entrance point for testing.
    Returns a tuple of (polygon, entrance_point)
    """
    # Sample coordinates (in lat/lon)
    coords = [
        (40.7128, -74.0060),  # New York City coordinates as example
        (40.7128, -74.0050),
        (40.7138, -74.0050),
        (40.7138, -74.0060)
    ]
    
    # Create polygon
    polygon = Polygon(coords)
    
    # Create entrance point on the first edge
    entrance_point = Point(40.7128, -74.0055)
    
    return polygon, entrance_point

def generate_lot_sketch(
    polygon: Polygon,
    entrance_point: Point,
    min_lot_size: float,
    min_front_line: float,
    road_width: float
) -> dict:
    """
    Generate a lot sketch within the given polygon.
    
    Args:
        polygon: Shapely Polygon object
        entrance_point: Shapely Point object representing the entrance
        min_lot_size: Minimum lot size in square feet
        min_front_line: Minimum front line length in feet
        road_width: Road width in feet
    
    Returns:
        GeoJSON dictionary containing the lots and roads
    """
    # Convert feet to meters for calculations
    min_lot_size_meters = square_feet_to_square_meters(min_lot_size)
    min_front_line_meters = feet_to_meters(min_front_line)
    road_width_meters = feet_to_meters(road_width)
    
    # Convert lat/lon to meters for calculations
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    
    # Transform polygon coordinates to meters
    x, y = polygon.exterior.xy
    coords_meters = [transformer.transform(lon, lat) for lon, lat in zip(x, y)]
    polygon_meters = Polygon(coords_meters)
    
    # Transform entrance point to meters
    entrance_x, entrance_y = transformer.transform(entrance_point.x, entrance_point.y)
    entrance_point_meters = Point(entrance_x, entrance_y)
    
    # TODO: Implement the actual lot generation algorithm
    # This is a placeholder that returns a simple grid of lots
    
    # Create a simple grid of lots (this is just a placeholder)
    lots = []
    roads = []
    
    # Convert back to lat/lon for GeoJSON
    transformer_reverse = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)
    
    # Create GeoJSON structure
    feature_collection = {
        "type": "FeatureCollection",
        "features": []
    }
    
    return feature_collection

def create_visualization(
    input_polygon: Polygon,
    entrance_point: Point,
    output_geojson: dict,
    output_file: str = "visualization.html"
) -> None:
    """
    Create an interactive map visualization of the input polygon and output sketch.
    
    Args:
        input_polygon: The input boundary polygon
        entrance_point: The entrance point
        output_geojson: The output GeoJSON containing lots and roads
        output_file: Name of the output HTML file
    """
    # Create a map centered on the polygon's centroid
    center_lat, center_lon = input_polygon.centroid.y, input_polygon.centroid.x
    m = folium.Map(location=[center_lat, center_lon], zoom_start=15)
    
    # Add the input polygon
    folium.GeoJson(
        geojson.Feature(
            geometry=geojson.Polygon([list(input_polygon.exterior.coords)]),
            properties={"name": "Input Boundary"}
        ),
        name="Input Boundary",
        style_function=lambda x: {
            "fillColor": "blue",
            "color": "blue",
            "fillOpacity": 0.1
        }
    ).add_to(m)
    
    # Add the entrance point
    folium.CircleMarker(
        location=[entrance_point.y, entrance_point.x],
        radius=5,
        color="red",
        fill=True,
        popup="Entrance Point"
    ).add_to(m)
    
    # Add the output lots and roads
    folium.GeoJson(
        output_geojson,
        name="Lots and Roads",
        style_function=lambda x: {
            "fillColor": "green",
            "color": "black",
            "fillOpacity": 0.3
        }
    ).add_to(m)
    
    # Add layer control
    folium.LayerControl().add_to(m)
    
    # Save the map
    m.save(output_file)

def get_float_input(prompt: str) -> float:
    """
    Get a float input from the user with error handling.
    
    Args:
        prompt: The input prompt to display
        
    Returns:
        float: The validated float input
    """
    while True:
        try:
            value = float(input(prompt))
            if value <= 0:
                print("Please enter a positive number.")
                continue
            return value
        except ValueError:
            print("Please enter a valid number.")

def print_polygon_area(polygon: Polygon):
    """
    Print the area of the polygon in square meters and square feet.
    """
    from pyproj import Transformer
    # Convert lat/lon to meters
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    x, y = polygon.exterior.xy
    coords_meters = [transformer.transform(lon, lat) for lon, lat in zip(x, y)]
    polygon_meters = Polygon(coords_meters)
    area_m2 = polygon_meters.area
    area_ft2 = area_m2 / 0.09290304
    print(f"Polygon area: {area_m2:.2f} square meters ({area_ft2:.2f} square feet)")

def main():
    # Check if input file exists, otherwise use sample polygon
    input_file = "input_polygon.geojson"
    if os.path.exists(input_file):
        print(f"Loading polygon from {input_file}")
        polygon, entrance_point = load_polygon_from_geojson(input_file)
    else:
        print("No input file found. Using sample polygon.")
        polygon, entrance_point = create_sample_polygon()
    
    # Print polygon area
    print_polygon_area(polygon)
    
    # Get user input
    print("\nPlease enter the following parameters:")
    min_lot_size = get_float_input("Minimum lot size (in square feet): ")
    min_front_line = get_float_input("Minimum front line length (in feet): ")
    road_width = get_float_input("Road width (in feet): ")
    
    # Generate lot sketch
    result = generate_lot_sketch(
        polygon,
        entrance_point,
        min_lot_size,
        min_front_line,
        road_width
    )
    
    # Save to GeoJSON file
    output_file = "lot_sketch.geojson"
    with open(output_file, 'w') as f:
        json.dump(result, f, indent=2)
    
    print(f"\nLot sketch has been saved to {output_file}")
    
    # Create visualization
    create_visualization(polygon, entrance_point, result)
    print("Visualization has been saved to visualization.html")

if __name__ == "__main__":
    main() 