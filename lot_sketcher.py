from typing import List, Tuple
from shapely.geometry import Polygon, Point, LineString, box
from shapely.ops import unary_union
from pyproj import Transformer
import math
import json
import os
from scipy.spatial import Voronoi
import numpy as np
import random
from shapely.ops import nearest_points

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
            polygon = Polygon(feature['geometry']['coordinates'][0])
        elif feature['geometry']['type'] == 'Point':
            entrance_point = Point(feature['geometry']['coordinates'])
    
    if polygon is None or entrance_point is None:
        raise ValueError("GeoJSON file must contain both a polygon and a point feature")
    
    return polygon, entrance_point

def create_sample_polygon() -> Tuple[Polygon, Point]:
    """
    Creates a sample 5-sided polygon and entrance point for testing.
    Returns a tuple of (polygon, entrance_point)
    """
    # Sample coordinates (in lat/lon)
    coords = [
        (40.7128, -74.0060),  # New York City coordinates as example
        (40.7128, -74.0050),
        (40.7138, -74.0040),  # Added point to make it 5-sided
        (40.7148, -74.0050),
        (40.7138, -74.0060)
    ]
    
    # Create polygon
    polygon = Polygon(coords)
    
    # Create entrance point on the first edge
    entrance_point = Point(40.7128, -74.0055)
    
    return polygon, entrance_point

def create_rectangular_lots(
    polygon: Polygon,
    min_lot_size: float,
    min_front_line: float
) -> List[Polygon]:
    """
    Create lots using a Voronoi-based approach with minimum lot size constraints.
    Ensures each lot has at least one edge meeting the minimum front line requirement.
    
    Args:
        polygon: The boundary polygon
        min_lot_size: Minimum lot size in square feet
        min_front_line: Minimum front line length in feet
    
    Returns:
        List of lot polygons
    """
    # Get the polygon's bounds
    minx, miny, maxx, maxy = polygon.bounds
    
    # Calculate lot dimensions
    lot_size = math.sqrt(min_lot_size * 1.5)  # Add buffer for spacing
    print(f"Target lot size: {lot_size:.2f} feet")
    
    # Calculate number of lots based on area
    num_lots = max(1, int(polygon.area / (min_lot_size * 1.5)))
    print(f"Target number of lots: {num_lots}")
    
    # Generate initial points within the polygon
    points = []
    remaining_polygon = polygon
    attempts = 0
    max_attempts = 1000
    
    while len(points) < num_lots and attempts < max_attempts:
        # Try to place a point in the remaining polygon
        if isinstance(remaining_polygon, Polygon):
            pieces = [remaining_polygon]
        else:
            pieces = list(remaining_polygon.geoms)
        
        if not pieces:
            break
            
        # Find the largest piece
        largest_piece = max(pieces, key=lambda p: p.area)
        
        # Generate a random point within the piece
        minx, miny, maxx, maxy = largest_piece.bounds
        x = minx + (maxx - minx) * random.random()
        y = miny + (maxy - miny) * random.random()
        point = Point(x, y)
        
        # Check if point is far enough from existing points
        if all(p.distance(point) >= lot_size for p in points):
            points.append(point)
            # Remove a small circle around the point from the remaining polygon
            buffer = point.buffer(lot_size / 2)
            remaining_polygon = remaining_polygon.difference(buffer)
            attempts = 0
        else:
            attempts += 1
    
    if not points:
        print("Could not generate initial points!")
        return []
    
    # Create Voronoi diagram
    points_array = np.array([[p.x, p.y] for p in points])
    vor = Voronoi(points_array)
    
    # Convert Voronoi regions to polygons
    lots = []
    for i, region in enumerate(vor.regions):
        if not region or -1 in region:  # Skip infinite regions
            continue
            
        # Get vertices of the region
        vertices = [vor.vertices[i] for i in region]
        if not vertices:
            continue
            
        # Create polygon from vertices
        lot = Polygon(vertices)
        
        # Intersect with original polygon
        lot = lot.intersection(polygon)
        
        if isinstance(lot, Polygon) and not lot.is_empty:
            # Check if lot meets size requirements
            if lot.area >= min_lot_size:
                # Check if any edge meets minimum front line requirement
                coords = list(lot.exterior.coords)
                has_valid_edge = False
                for j in range(len(coords) - 1):
                    edge = LineString([coords[j], coords[j + 1]])
                    if edge.length >= min_front_line:
                        has_valid_edge = True
                        break
                
                if has_valid_edge:
                    lots.append(lot)
                    print(f"Created lot with area {lot.area:.2f} sq ft")
    
    # If we didn't create enough lots, try to fill remaining spaces
    if len(lots) < num_lots:
        print("Filling remaining spaces...")
        remaining_polygon = polygon
        for lot in lots:
            remaining_polygon = remaining_polygon.difference(lot)
        
        if isinstance(remaining_polygon, Polygon):
            pieces = [remaining_polygon]
        else:
            pieces = list(remaining_polygon.geoms)
        
        for piece in pieces:
            if piece.area >= min_lot_size:
                # Try to create a lot from the center
                center = piece.centroid
                # Create a rectangular lot with minimum front line
                half_width = min_front_line / 2
                half_depth = min_lot_size / min_front_line / 2
                lot = box(
                    center.x - half_width,
                    center.y - half_depth,
                    center.x + half_width,
                    center.y + half_depth
                )
                
                lot = lot.intersection(piece)
                if isinstance(lot, Polygon) and not lot.is_empty:
                    if lot.area >= min_lot_size:
                        lots.append(lot)
                        print(f"Created additional lot with area {lot.area:.2f} sq ft")
    
    print(f"Created {len(lots)} lots")
    return lots

def create_road_network(
    polygon: Polygon,
    entrance_point: Point,
    lots: List[Polygon],
    road_width: float
) -> List[LineString]:
    """
    Create a road network using lot edges.
    Ensures each lot has at least one edge connected to the road network.
    
    Args:
        polygon: The boundary polygon
        entrance_point: The entrance point
        lots: List of lot polygons
        road_width: Road width in feet
    
    Returns:
        List of road centerlines
    """
    if not lots:
        print("No lots to connect!")
        return []
    
    print(f"Creating road network for {len(lots)} lots")
    
    # Get lot edges and track which lots are connected
    lot_edges = []
    edge_to_lot = {}  # Map edges to their lot
    lot_connected = {i: False for i in range(len(lots))}  # Track if lot has any edge connected
    
    for i, lot in enumerate(lots):
        coords = list(lot.exterior.coords)
        for j in range(len(coords) - 1):
            edge = LineString([coords[j], coords[j + 1]])
            if edge.length >= road_width * 2:  # Only consider edges long enough for roads
                lot_edges.append((edge, i))
                edge_to_lot[edge] = i
    
    # Sort edges by length
    lot_edges.sort(key=lambda x: x[0].length)
    
    # Use Kruskal's algorithm to find MST
    parent = {}
    for edge, _ in lot_edges:
        parent[edge.coords[0]] = edge.coords[0]
        parent[edge.coords[1]] = edge.coords[1]
    
    def find(v):
        if parent[v] != v:
            parent[v] = find(parent[v])
        return parent[v]
    
    def union(v1, v2):
        parent[find(v1)] = find(v2)
    
    roads = []
    
    # First pass: Ensure each lot has at least one edge connected
    for edge, lot_idx in lot_edges:
        if not lot_connected[lot_idx]:
            if polygon.contains(edge) or polygon.touches(edge):
                roads.append(edge)
                lot_connected[lot_idx] = True
                union(edge.coords[0], edge.coords[1])
    
    # Second pass: Add remaining edges to complete MST
    for edge, lot_idx in lot_edges:
        if find(edge.coords[0]) != find(edge.coords[1]):
            if polygon.contains(edge) or polygon.touches(edge):
                roads.append(edge)
                union(edge.coords[0], edge.coords[1])
    
    # Add entrance connection
    if roads:
        entrance_road = None
        min_dist = float('inf')
        for road in roads:
            dist = road.distance(entrance_point)
            if dist < min_dist:
                min_dist = dist
                entrance_road = road
        
        if entrance_road:
            nearest_point = nearest_points(entrance_point, entrance_road)[1]
            entrance_connection = LineString([entrance_point, nearest_point])
            roads.append(entrance_connection)
    
    # Add additional roads to ensure all lots have access
    for i, lot in enumerate(lots):
        if not lot_connected[i]:
            # Find the nearest road
            min_dist = float('inf')
            nearest_road = None
            for road in roads:
                dist = lot.distance(road)
                if dist < min_dist:
                    min_dist = dist
                    nearest_road = road
            
            if nearest_road is not None:
                # Find the nearest edge of the lot to the nearest road
                coords = list(lot.exterior.coords)
                lot_edges = [
                    LineString([coords[j], coords[j + 1]]) 
                    for j in range(len(coords) - 1)
                ]
                nearest_edge = min(lot_edges, key=lambda e: e.distance(nearest_road))
                nearest_point = nearest_points(nearest_edge, nearest_road)[0]
                connecting_road = LineString([
                    nearest_point,
                    nearest_points(nearest_point, nearest_road)[1]
                ])
                roads.append(connecting_road)
                lot_connected[i] = True
    
    # Merge overlapping or connected roads
    merged_roads = []
    for road in roads:
        # Check if this road connects to any existing road
        connected = False
        for i, existing_road in enumerate(merged_roads):
            if (road.touches(existing_road) or 
                    road.distance(existing_road) < road_width):
                # Merge the roads
                merged_roads[i] = unary_union([road, existing_road])
                connected = True
                break
        if not connected:
            merged_roads.append(road)
    
    # Convert merged roads back to LineStrings
    final_roads = []
    for road in merged_roads:
        if isinstance(road, LineString):
            final_roads.append(road)
        else:
            # If the road is a MultiLineString, split it into LineStrings
            for line in road.geoms:
                final_roads.append(line)
    
    # Ensure all roads are within the polygon
    final_roads = [
        road for road in final_roads 
        if polygon.contains(road) or polygon.touches(road)
    ]
    
    total_length = sum(road.length for road in final_roads)
    print(f"Created {len(final_roads)} road segments with total length {total_length:.2f} feet")
    return final_roads

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
    print("\nUsing measurements in feet:")
    print(f"Minimum lot size: {min_lot_size:.2f} sq ft")
    print(f"Minimum front line: {min_front_line:.2f} ft")
    print(f"Road width: {road_width:.2f} ft")
    
    # Check if coordinates are in lat/lon
    x, y = polygon.exterior.xy
    is_latlon = abs(x[0]) <= 180 and abs(y[0]) <= 90
    
    if is_latlon:
        print("\nConverting coordinates from lat/lon to feet")
        # Convert lat/lon to meters first
        transformer = Transformer.from_crs(
            "EPSG:4326",
            "EPSG:3857",
            always_xy=True
        )
        
        # Transform polygon coordinates to meters
        coords_meters = [
            transformer.transform(lon, lat)
            for lon, lat in zip(x, y)
        ]
        polygon_meters = Polygon(coords_meters)
        
        # Transform entrance point to meters
        entrance_x, entrance_y = transformer.transform(
            entrance_point.x,
            entrance_point.y
        )
        entrance_point_meters = Point(entrance_x, entrance_y)
        
        # Convert meters to feet
        x_feet = [x * 3.28084 for x in polygon_meters.exterior.xy[0]]
        y_feet = [y * 3.28084 for y in polygon_meters.exterior.xy[1]]
        polygon_feet = Polygon(zip(x_feet, y_feet))
        entrance_point_feet = Point(
            entrance_point_meters.x * 3.28084,
            entrance_point_meters.y * 3.28084
        )
    else:
        print("\nCoordinates appear to be already in feet")
        polygon_feet = polygon
        entrance_point_feet = entrance_point
    
    print("\nPolygon bounds:", polygon_feet.bounds)
    print(
        "Entrance point:",
        (entrance_point_feet.x, entrance_point_feet.y)
    )
    
    # Convert to local coordinates (0,0 based)
    minx, miny, maxx, maxy = polygon_feet.bounds
    local_coords = []
    for x, y in polygon_feet.exterior.coords:
        local_coords.append((x - minx, y - miny))
    polygon_local = Polygon(local_coords)
    entrance_point_local = Point(
        entrance_point_feet.x - minx,
        entrance_point_feet.y - miny
    )
    
    # First create lots
    lots = create_rectangular_lots(
        polygon_local,
        min_lot_size,
        min_front_line
    )
    
    if not lots:
        print("No lots were created!")
        return {"type": "FeatureCollection", "features": []}
    
    # Then create road network
    roads = create_road_network(
        polygon_local,
        entrance_point_local,
        lots,
        road_width
    )
    
    print(f"Created {len(roads)} roads")
    
    # Create GeoJSON structure
    feature_collection = {
        "type": "FeatureCollection",
        "features": []
    }
    
    # Check if we need to convert back to lat/lon
    if is_latlon:
        print("\nConverting coordinates back to lat/lon")
        transformer_reverse = Transformer.from_crs(
            "EPSG:3857",
            "EPSG:4326",
            always_xy=True
        )
        
        # Add roads to GeoJSON
        for road in roads:
            coords = []
            for x, y in road.coords:
                # Convert back to original coordinates
                x += minx
                y += miny
                # Convert feet back to meters
                x_meters = x / 3.28084
                y_meters = y / 3.28084
                # Convert meters back to lat/lon
                lon, lat = transformer_reverse.transform(x_meters, y_meters)
                coords.append([lon, lat])
            
            feature = {
                "type": "Feature",
                "properties": {
                    "type": "road",
                    "width": road_width / 3.28084,  # Convert to meters
                    "stroke": "#FF0000",
                    "stroke-width": 3,
                    "stroke-opacity": 1
                },
                "geometry": {
                    "type": "LineString",
                    "coordinates": coords
                }
            }
            feature_collection["features"].append(feature)
        
        # Add lots to GeoJSON
        for lot in lots:
            coords = []
            for x, y in lot.exterior.coords:
                # Convert back to original coordinates
                x += minx
                y += miny
                # Convert feet back to meters
                x_meters = x / 3.28084
                y_meters = y / 3.28084
                # Convert meters back to lat/lon
                lon, lat = transformer_reverse.transform(x_meters, y_meters)
                coords.append([lon, lat])
            
            feature = {
                "type": "Feature",
                "properties": {
                    "type": "lot",
                    "area": lot.area,  # Already in square feet
                    "stroke": "#000000",
                    "stroke-width": 1,
                    "stroke-opacity": 0.5,
                    "fill": "#FFFFFF",
                    "fill-opacity": 0.1
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [coords]
                }
            }
            feature_collection["features"].append(feature)
    else:
        # Coordinates are already in feet
        # Add roads to GeoJSON
        for road in roads:
            coords = []
            for x, y in road.coords:
                # Convert back to original coordinates
                x += minx
                y += miny
                coords.append([x, y])
            
            feature = {
                "type": "Feature",
                "properties": {
                    "type": "road",
                    "width": road_width,
                    "stroke": "#FF0000",
                    "stroke-width": 3,
                    "stroke-opacity": 1
                },
                "geometry": {
                    "type": "LineString",
                    "coordinates": coords
                }
            }
            feature_collection["features"].append(feature)
        
        # Add lots to GeoJSON
        for lot in lots:
            coords = []
            for x, y in lot.exterior.coords:
                # Convert back to original coordinates
                x += minx
                y += miny
                coords.append([x, y])
            
            feature = {
                "type": "Feature",
                "properties": {
                    "type": "lot",
                    "area": lot.area,  # Already in square feet
                    "stroke": "#000000",
                    "stroke-width": 1,
                    "stroke-opacity": 0.5,
                    "fill": "#FFFFFF",
                    "fill-opacity": 0.1
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [coords]
                }
            }
            feature_collection["features"].append(feature)
    
    return feature_collection

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
    Print the area of the polygon in square feet.
    """
    # Check if coordinates are in lat/lon
    x, y = polygon.exterior.xy
    if abs(x[0]) <= 180 and abs(y[0]) <= 90:
        print("\nConverting coordinates from lat/lon to feet")
        # Convert lat/lon to meters first
        transformer = Transformer.from_crs(
            "EPSG:4326",
            "EPSG:3857",
            always_xy=True
        )
        
        # Transform polygon coordinates to meters
        coords_meters = [
            transformer.transform(lon, lat)
            for lon, lat in zip(x, y)
        ]
        polygon_meters = Polygon(coords_meters)
        
        # Convert meters to feet
        x_feet = [x * 3.28084 for x in polygon_meters.exterior.xy[0]]
        y_feet = [y * 3.28084 for y in polygon_meters.exterior.xy[1]]
        polygon_feet = Polygon(zip(x_feet, y_feet))
    else:
        print("\nCoordinates appear to be already in feet")
        polygon_feet = polygon
    
    area_ft2 = polygon_feet.area
    print(f"Polygon area: {area_ft2:.2f} square feet")

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

if __name__ == "__main__":
    main() 