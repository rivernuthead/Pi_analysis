import rasterio
import numpy as np
import os
from math import cos, sin, radians
import matplotlib.pyplot as plt

def generate_cross_section(dem_path, start_point, length, angle_deg, output_txt):
    """
    Generates a cross-section from a DEM file.

    Parameters:
        dem_path (str): Path to the DEM `.tif` file.
        start_point (tuple): (x, y) coordinates of the starting point of the cross-section.
        length (float): Length of the cross-section in map units.
        angle_deg (float): Angle of the cross-section in degrees (from horizontal).
        output_txt (str): Path to the output text file.

    Returns:
        list of tuples: Filtered coordinates of the cross-section line.
    """
    # Load the DEM
    with rasterio.open(dem_path) as src:
        dem = src.read(1)  # Read the first band
        transform = src.transform

    # Convert angle to radians
    angle_rad = radians(angle_deg)

    # Calculate the x and y components of the cross-section line
    dx = cos(angle_rad)
    dy = sin(angle_rad)

    # Prepare output data
    distances = []
    elevations = []
    cross_section_coords = []

    # Iterate along the cross-section
    for dist in np.arange(0, length, step=0.5):  # 100 points along the cross-section
        x = start_point[0] + dist * dx
        y = start_point[1] + dist * dy

        # Convert map coordinates to row, col indices
        row, col = rasterio.transform.rowcol(transform, x, y)

        try:
            elevation = dem[row, col]
        except IndexError:
            # If the coordinates fall outside the DEM bounds
            elevation = np.nan

        # Calculate the distance in real units considering grid resolution
        real_distance = dist * 0.5  # Convert grid distance to map distance by multiplying by the resolution

        distances.append(real_distance)
        elevations.append(elevation)
        cross_section_coords.append((x, y))

    # Filter out NaN values and values less than 0
    filtered_data = [(dist, elev, coord) for dist, elev, coord in zip(distances, elevations, cross_section_coords) if not np.isnan(elev) and elev >= 0]

    # Recalculate distances starting from 0 with a step of 0.5
    filtered_distances = np.arange(0, len(filtered_data) * 0.5, 0.5)
    filtered_elevations = [elev for _, elev, _ in filtered_data]
    filtered_coords = [coord for _, _, coord in filtered_data]

    # Set elevation at the start and end to 360
    if len(filtered_elevations) > 0:
        filtered_elevations[0] = 360
        filtered_elevations[-1] = 360

    # Write to output file
    with open(output_txt, 'w') as f:
        f.write("Cross-sections Rees\n")
        f.write(f"1\t{len(filtered_elevations)}\t1\n")

        for dist, elevation in zip(filtered_distances, filtered_elevations):
            f.write(f"{dist:.2f}\t{elevation:.2f}\t30\n")

    return filtered_coords, filtered_distances, filtered_elevations

def plot_dems_with_cross_sections(dem_paths, all_cross_sections_coords):
    """
    Plots the DEMs with multiple cross-section lines.

    Parameters:
        dem_paths (list of str): Paths to the DEM `.tif` files.
        all_cross_sections_coords (list of list of tuples): List of cross-section coordinates for each cross-section.

    Returns:
        None
    """
    fig, ax = plt.subplots(1, len(dem_paths), figsize=(10 * len(dem_paths), 10))

    if len(dem_paths) == 1:
        ax = [ax]  # Ensure ax is iterable if there's only one DEM

    for idx, dem_path in enumerate(dem_paths):
        with rasterio.open(dem_path) as src:
            dem = src.read(1)
            transform = src.transform

            # Plot the DEM
            dem = np.where(dem<0, np.min(dem[dem>0])-2, dem)
            ax[idx].imshow(dem, cmap='terrain', extent=(
                transform[2],
                transform[2] + dem.shape[1] * transform[0],
                transform[5] + dem.shape[0] * transform[4],
                transform[5]
            ))

            ax[idx].set_title(f"DEM {idx + 1}")
            ax[idx].set_xlabel("X Coordinate")
            ax[idx].set_ylabel("Y Coordinate")

            # Plot all cross-section lines
            for cross_section_coords in all_cross_sections_coords[idx]:
                cross_x, cross_y = zip(*cross_section_coords)
                ax[idx].plot(cross_x, cross_y, linewidth=2, label='Cross-Section')

            ax[idx].legend()

    plt.tight_layout()
    plt.show()

def plot_all_cross_section_profiles(dem_paths, all_cross_section_profiles):
    """
    Plots elevation profiles of all cross-sections with unique colors for each DEM.

    Parameters:
        dem_paths (list of str): Paths to the DEM `.tif` files.
        all_cross_section_profiles (list of list of tuples): List of distances and elevations for each DEM's cross-sections.

    Returns:
        None
    """
    colors = plt.cm.viridis(np.linspace(0, 1, len(dem_paths)))  # Generate unique colors for each DEM

    plt.figure(figsize=(12, 8))

    for idx, (dem_path, cross_section_profiles) in enumerate(zip(dem_paths, all_cross_section_profiles)):
        for distances, elevations in cross_section_profiles:
            plt.plot(distances, elevations, label=f'DEM {idx + 1}', color=colors[idx])

    plt.title("Cross-Section Elevation Profiles")
    plt.xlabel("Distance (m)")
    plt.ylabel("Elevation (m)")
    plt.legend()
    plt.grid()
    plt.show()

def generate_multiple_cross_sections(dem_path, start_point, length, angle_deg, separation, num_sections, output_dir):
    """
    Generates multiple cross-sections from a DEM file.

    Parameters:
        dem_path (str): Path to the DEM `.tif` file.
        start_point (tuple): (x, y) coordinates of the starting point of the first cross-section.
        length (float): Length of each cross-section in map units.
        angle_deg (float): Angle of the cross-sections in degrees (from horizontal).
        separation (float): Distance between consecutive cross-sections.
        num_sections (int): Number of cross-sections to generate.
        output_dir (str): Directory to save output text files.

    Returns:
        list of list of tuples: Filtered coordinates and profiles of all cross-sections.
    """
    cross_sections_coords = []
    cross_sections_profiles = []
    dem_filename = dem_path[-13:-4]
    for i in range(num_sections):
        # Offset the starting point for each cross-section
        offset_x = i * separation * cos(radians(angle_deg + 90))
        offset_y = i * separation * sin(radians(angle_deg + 90))
        current_start_point = (start_point[0] + offset_x, start_point[1] + offset_y)

        output_txt = os.path.join(output_dir, f"{dem_filename}_cross_section_{i + 1}.txt")

        coords, distances, elevations = generate_cross_section(
            dem_path=dem_path,
            start_point=current_start_point,
            length=length,
            angle_deg=angle_deg,
            output_txt=output_txt
        )

        cross_sections_coords.append(coords)
        cross_sections_profiles.append((distances, elevations))

    return cross_sections_coords, cross_sections_profiles

#Script
home_folder = os.getcwd()
path_DEMs = os.path.join(home_folder, 'Data', 'DEMS')
dem_paths = [os.path.join(path_DEMs, "e00_clean.tif")]#,os.path.join(path_DEMs, "e10_clean.tif")] 
start_point = (1236544, 5032389)  # Replace with your starting point coordinates
length = 1000  # Length of each cross-section in map units
angle_deg = -35  # Angle in degrees
separation = -200  # Separation distance between cross-sections in map units
num_sections = 11  # Number of cross-sections to generate
output_dir = os.path.join(home_folder,'Data', "Cross_sections")  # Output directory

os.makedirs(output_dir, exist_ok=True)

all_cross_sections_coords = []
all_cross_section_profiles = []

for dem_path in dem_paths:
    cross_sections_coords, cross_sections_profiles = generate_multiple_cross_sections(
        dem_path=dem_path,
        start_point=start_point,
        length=length,
        angle_deg=angle_deg,
        separation= separation, num_sections=num_sections, output_dir=output_dir)
    all_cross_sections_coords.append(cross_sections_coords)
    all_cross_section_profiles.append(cross_sections_profiles)

# Plot DEMs with cross-section lines
plot_dems_with_cross_sections(dem_paths, all_cross_sections_coords)

# Plot elevation profiles of all cross-sections
plot_all_cross_section_profiles(dem_paths, all_cross_section_profiles)
