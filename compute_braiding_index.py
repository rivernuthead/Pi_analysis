import os
import numpy as np
import matplotlib.pyplot as plt
from skimage import morphology

####### Functions ##################################################################################
def read_data(file_path):
    """
    Reads the DEM data from a text file.

    Parameters:
    file_path (str): Path to the DEM file.

    Returns:
    tuple: A tuple containing the headers and the DEM data array.
    """
    with open(file_path, 'r') as f:
        # Read the first eight lines as headers
        headers = [next(f) for _ in range(8)]
        # Read the rest of the data into a numpy array
        data = np.loadtxt(f)
        # Replace -999 values with NaN
        data = np.where(data == -999, np.nan, data)
    return headers, data

def update_headers(headers, east_coord):
    # Update rows and cols in headers
    headers[2] = f"east: {east_coord}\n"
    return headers
    
def generate_elevation_array(shape, max_elevation, slope, x_res):
    """
    Generates a synthetic 2D array of elevation with a specified maximum elevation and slope.

    Parameters:
    shape (tuple): Shape of the array (y, x).
    max_elevation (float): Maximum elevation value at the first x column.
    slope (float): Slope to be applied along the x-axis.
    x_res (float): Resolution along the x-axis.

    Returns:
    np.ndarray: 2D array of synthetic elevation values.
    """
    y, x = shape
    elevation_array = np.zeros((y, x))
    
    # Create a linear gradient along the x-axis with the given resolution
    for i in range(x):
        elevation_array[:, i] = max_elevation - slope * i * x_res
    
    return elevation_array

def compute_above_length(dem, synthetic_elevation, y_res):
    """
    Computes the length of segments where the synthetic elevation values are above the DEM values.

    Parameters:
    dem (np.ndarray): Original DEM array.
    synthetic_elevation (np.ndarray): Synthetic elevation array.
    y_res (float): Resolution along the y-axis.

    Returns:
    np.ndarray: Array containing lengths of segments where synthetic values are above DEM values.
    """
    # Count the number of times synthetic elevation is greater than DEM along each column
    above_length = np.nansum(synthetic_elevation > dem, axis=0) * y_res
    return above_length

def plot_profiles_and_above_length(dem, synthetic_elevation, above_length, x_res, y_res):
    """
    Plots the elevation profiles for the DEM, synthetic elevation, and the lengths of segments 
    where synthetic values are above DEM values.

    Parameters:
    dem (np.ndarray): Original DEM array.
    synthetic_elevation (np.ndarray): Synthetic elevation array.
    above_length (np.ndarray): Lengths of segments where synthetic values are above DEM values.
    x_res (float): Resolution along the x-axis.
    y_res (float): Resolution along the y-axis.
    """
    y = np.arange(dem.shape[0]) * y_res
    x = np.arange(dem.shape[1]) * x_res
    num_profil = 200
    y_dem = dem[:, num_profil]  # Taking the profile from the 200th column (adjust as needed)
    y_synthetic = synthetic_elevation[:, num_profil]  # Taking the profile from the 200th column (adjust as needed)

    # Calculate the difference between synthetic elevation and DEM
    difference = synthetic_elevation - dem
    
    # Identify where the sign changes along each column
    sign_changes = np.diff(np.sign(difference), axis=0)
    
    # Create a mask for the sign change points
    sign_change_points_mask = np.abs(sign_changes) == 2
    y_mask, x_mask = np.where(sign_change_points_mask)
    y_mask = y_mask[x_mask==num_profil]
    x_points = y_mask 
    y_points = []
    for i in x_points:
        elevation = y_dem[i]
        y_points = np.append(y_points,elevation)
    
    max_y = np.nanmax(dem)
    min_y = np.nanmin(dem)
    
    plt.figure(figsize=(12, 12))
    
    plt.subplot(3, 1, 1)
    plt.plot(y, y_dem,'k.', label='DEM', color='black',ls='-',lw=1)
    plt.plot(y, y_synthetic, 'b.', label='Synthetic Elevation', color='blue',ls='-',lw=1)
    plt.scatter(x_points* y_res, y_points, color='red', s=20, label='Sign Change Points')
    plt.xlabel('Y Axis (mm)')
    plt.ylabel('Elevation')
    plt.title('Elevation Profile along the Y Axis')
    plt.ylim(min_y,max_y)
    plt.legend()
    plt.grid(True)
    
    plt.subplot(3, 1, 2)
    plt.plot(x, above_length, label='Length of Synthetic Values Above DEM', color='green')
    plt.xlabel('X Axis (mm)')
    plt.ylabel('Length (mm)')
    plt.title('Length of Segments where Synthetic Values are Above DEM')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

def find_best_max_elevation(dem, border_mask, slope, x_res, y_res, mean_wetted_width, max_elevation_start, max_elevation_end, step):
    best_max_elevation = None
    best_mean_above_length = None
    min_diff = float('inf')
    for max_elevation in np.arange(max_elevation_start, max_elevation_end, step):
        synthetic_elevation = generate_elevation_array(dem.shape, max_elevation, slope, x_res)
        synthetic_elevation = np.where(border_mask == 0, np.nan, synthetic_elevation)
        above_length = compute_above_length(dem, synthetic_elevation, y_res)
        mean_above_length = np.nanmean(above_length)
        
        diff = abs(mean_above_length - mean_wetted_width)
        
        if diff < min_diff:
            min_diff = diff
            best_max_elevation = max_elevation
            best_mean_above_length = mean_above_length
        
        if min_diff == 0:
            break
    
    return best_max_elevation, best_mean_above_length

def compute_water_depth(dem, synthetic_elevation):
    """
    Computes the water depth as the difference between synthetic elevation and DEM.

    Parameters:
    dem (np.ndarray): Original DEM array.
    synthetic_elevation (np.ndarray): Synthetic elevation array.

    Returns:
    np.ndarray: Array containing the water depth values.
    """
    water_depth = synthetic_elevation - dem
    water_depth = np.where(water_depth < 0, 0, water_depth)
    return water_depth

def save_water_depth(file_path, headers, water_depth):
    """
    Saves the water depth array to a text file with the given headers.

    Parameters:
    file_path (str): Path to the output file.
    headers (list): List of header lines to be written at the top of the file.
    water_depth (np.ndarray): Water depth array to be saved.
    """
    with open(file_path, 'w') as f:
        for header in headers:
            f.write(header)
        np.savetxt(f, water_depth, fmt='%.6f')

def plot_dem_with_sign_change_points(dem, synthetic_elevation):
    """
    Plots the DEM with points highlighted in red where the sign of the difference
    between DEM and synthetic elevation changes.

    Parameters:
    dem (np.ndarray): Original DEM array.
    synthetic_elevation (np.ndarray): Synthetic elevation array.
    """
    # Calculate the difference between synthetic elevation and DEM
    difference = synthetic_elevation - dem
    
    # Identify where the sign changes along each column
    sign_changes = np.diff(np.sign(difference), axis=0)
    
    # Create a mask for the sign change points
    sign_change_points_mask = np.abs(sign_changes) == 2
    
    plt.figure(figsize=(10, 8))
    plt.imshow(dem, cmap='terrain', interpolation='none')
    plt.colorbar(label='Elevation')
    
    y, x = np.where(sign_change_points_mask)
    plt.scatter(x, y, color='red', s=10, label='Sign Change Points')
    
    plt.xlabel('X Axis')
    plt.ylabel('Y Axis')
    plt.title('DEM with Sign Change Points for Water Surface Intersection')
    plt.legend()
    plt.show()


def filter_water_depth_mask(array, y_res, min_width_ones=3, min_width_zeros=2):
    """
    Filter the water depth mask to keep only channels with wetted width > 3 mm.

    Parameters:
    mask (np.ndarray): Water depth mask.
    """
    y, x = array.shape
    min_width_ones = min_width_ones // y_res # transform mm in pixels
    min_width_zeros = min_width_zeros // y_res # transform mm in pixels
    for j in range(x):  # iterate over each column
        consecutive_ones = 0
        for i in range(y):  # iterate over each element in the column
            if array[i, j] == 1:
                consecutive_ones += 1
            else:
                if consecutive_ones < min_width_ones:  # if consecutive 1s is less than 3
                    for k in range(i - consecutive_ones, i):
                        array[k, j] = 0
                consecutive_ones = 0
        if consecutive_ones < min_width_ones:  # handle case where consecutive 1s continue till the end of column
            for k in range(y - consecutive_ones, y):
                array[k, j] = 0

    # Deals with zeros          
    for j in range(x):  # iterate over each column
        consecutive_zeros = 0
        for i in range(y):  # iterate over each element in the column
            if array[i, j] == 0:
                consecutive_zeros +=1
            else:
                if consecutive_zeros < min_width_zeros:  # if consecutive 1s is less than 3
                    for k in range(i - consecutive_zeros, i):
                        array[k, j] = 1
                consecutive_zeros = 0
                
        if consecutive_zeros < min_width_zeros:  # handle case where consecutive 1s continue till the end of column
            for k in range(y - consecutive_zeros, y):
                array[k, j] = 1
    return array


def compute_braiding_index(array, step=1):
    # Ensure the array is a NumPy array
    array = np.array(array)
    
    # Get the shape of the array
    y, x = array.shape
    
    # Initialize a list to store the count of changes for each column
    nb_measurement = 0
    for col in range(0,x,step):
        nb_measurement = nb_measurement + 1

    nb_channels = [0] * nb_measurement

    # Iterate through each column
    for i,col in enumerate(range(0,x,step)):
        # Start with the initial condition check
        if array[0, col] == 1:
            nb_channels[i] += 1
            
        # Check for changes from 0 to 1 in the column
        for row in range(1, y):
            if array[row-1, col] == 0 and array[row, col] == 1:
                nb_channels[i] += 1
                
    braiding_index = np.mean(nb_channels)  
    
    return braiding_index

####### Script ##################################################################################

# Options 
Active_braiding_index = 1
Total_braiding_index = 0

# Parameters
folder_home = os.getcwd()
slope = 0.00037  # longitudinal Slope applied along the x-axis (only for Total_braiding index)
data_name = 'TPiQs'



if Active_braiding_index == 1:
    # Parameters
    runs = ['q05rgm5','q07rgm','q10rgm2','q15rgm2','q20rgm2']
    set_names = ['q05_1', 'q07_1', 'q10_2', 'q15_2', 'q20_2']
    x_res = 5  # Resolution along the x-axis in mm
    y_res = 5  # Resolution along the y-axis in mm
    cs_steps = [100] # step between two cross sections in mm. Minimum is 5.
    thrs_rso = 1600 # number of pixels ; consider a square of n*n pixels ; 400= 10 cm * 10 cm ; 1600 = 20 cm * 20 cm
    path_out_fig =  os.path.join(folder_home, 'bedload_analysis', 'plots', 'Active_braiding_index')
    if not(os.path.exists(path_out_fig)):
        os.makedirs(path_out_fig)
        
    for ii,run in enumerate(runs):
        print('*****************')
        print(run)
        print('*****************')
        # Load images
        path_images = os.path.join(folder_home, 'bedload_analysis' ,'output_data','3_PiQs_BAW_2Dt_filter', set_names[ii])
        stack = np.load(os.path.join(path_images,run+'_BAA_stack_LR5_cld.npy'),mmap_mode='r+')
        time, dy, dx = np.shape(stack)
        stack = np.where(stack>0,1,0) # create boolean
    
        # Compute active braiding index
        active_braiding_index_steps = []
        std = []
        for s in cs_steps:
            cs_step_res = int(s//x_res)
            active_braiding_index_list = []
            for i in range(0,time): # for each image
                image = stack[i,:,:]
                
                # 4.d REMOVE SMALL OBJECTS
                image_rso = morphology.remove_small_objects(image>0, min_size=thrs_rso, connectivity=1)
                image = np.where(np.isnan(image)==False,image_rso,image)
                # Compute braiding index
                active_braiding_index_im = compute_braiding_index(image,step=cs_step_res)
                mean_ABI = np.mean(active_braiding_index_im)
                #if mean_ABI<1:
                    #mean_ABI = 1
                active_braiding_index_list = np.append(active_braiding_index_list,mean_ABI)
                
            active_braiding_index_steps = np.append(active_braiding_index_steps,np.round(np.mean(active_braiding_index_list),decimals=2))
            std = np.append(std,np.round(np.std(active_braiding_index_list),decimals=2))
            #Plot
            fig, ax = plt.subplots(figsize=(10, 8))
            x = np.arange(0,time,1)
            ax.plot(x, active_braiding_index_list,'k-')
            ax.set_xlabel('time')
            ax.set_ylabel('ABI')
            ax.set_ylim(0,2.5)
            ax.set_title(f'Reach_scale mean ABI through time ({run})')
            # Add text outside the plot using figure coordinates (0.1 is bottom-left, 0.9 is top-right)
            fig.text(0.1, 0.10, f" For each image, the Nb of channel is computed along cross-sections every {s} mm.\n Each image has been filtered by filtering out small objects (thr_rso = {thrs_rso}) to keep only long channels \n Mean ABI = {active_braiding_index_steps}; Std ABI = {std}", fontsize=11)
            plt.subplots_adjust(bottom=0.25)
            # Save the plot with 300 dpi
            plt.savefig(os.path.join(path_out_fig,f'Reach_scale_mean_ABI_through_time_{run}.pdf'), dpi=300, bbox_inches='tight')
            #plt.tight_layout()
            plt.show()
            
            
        #active_braiding_index = np.mean(active_braiding_index_steps)
        print(f"Active braiding index is {active_braiding_index_steps}")
            
if Total_braiding_index == 1:
    # Parameters
    dem_name = 'matrix_bed_W06_Q05_s9.txt'
    path_dem = os.path.join(folder_home, data_name ,'Laser_survey','W06_Q05', 'DEMs', 'survey_9', dem_name)
    mask_name = 'mask_borders_W06_Q05.txt'
    path_border_mask = os.path.join(folder_home,data_name ,'Laser_survey', 'W06_Q05', mask_name)
    path_water_depth = os.path.join(folder_home,data_name ,'Laser_survey', 'W06_Q05', 'Water_depth')
    x_res = 50  # Resolution along the x-axis in mm
    y_res = 5  # Resolution along the y-axis in mm
    water_depth_filename = dem_name[:-4] + '_water_depth.txt'
    east_coord = 13807 # to plot the DoD with the actual extent
    mean_wetted_width = 510 # User-defined mean wetted width in mm

    # Load DEM
    headers, DEM = read_data(path_dem)
    shape = np.shape(DEM)

    # Load border mask
    headers, border_mask = read_data(path_border_mask)

    # Apply mask
    DEM = np.where(border_mask == 0, np.nan, DEM)

    # Find the best max_elevation
    max_elevation_start = np.nanmax(DEM[:, 0])  # Starting max elevation
    max_elevation_end = np.nanmin(DEM[:, 0])  # Ending max elevation
    step = -0.1  # Step size for max elevation
    
    best_max_elevation, best_mean_above_length = find_best_max_elevation(DEM, border_mask, slope, x_res, y_res, mean_wetted_width, max_elevation_start, max_elevation_end, step)
    
    # Generate synthetic elevation array with the best max elevation
    synthetic_elevation = generate_elevation_array(shape, best_max_elevation, slope, x_res)
    synthetic_elevation = np.where(border_mask == 0, np.nan, synthetic_elevation)
    
    # Compute lengths of segments where synthetic values are above DEM values
    above_length = compute_above_length(DEM, synthetic_elevation, y_res)
    
    # Plot the profiles and above length
    plot_profiles_and_above_length(DEM, synthetic_elevation, above_length, x_res, y_res)
    
    # Plot DEM with closest points
    plot_dem_with_sign_change_points(DEM, synthetic_elevation)
    
    # Compute water depth
    water_depth = compute_water_depth(DEM, synthetic_elevation)
    water_depth = np.where(np.isnan(water_depth), -999, water_depth)
    
    # Save water depth to a .txt file
    if not os.path.exists(path_water_depth):
        os.mkdir(path_water_depth)
    #save dem
    update_headers(headers, east_coord)
    save_water_depth(os.path.join(path_water_depth, dem_name), headers, DEM)
    
    save_water_depth(os.path.join(path_water_depth, water_depth_filename), headers, water_depth)
    
    # Water depth mask
    water_depth_mask = np.where(water_depth>0,1,0) # value is in mm
    water_depth_mask_filename = dem_name[:-4] + '_water_depth_mask.txt'
    
    save_water_depth(os.path.join(path_water_depth, water_depth_mask_filename), headers, water_depth_mask)
    
    # Filter water_depth mask
    water_depth_filtered = filter_water_depth_mask(water_depth_mask, y_res,min_width_ones=50,min_width_zeros=10)
    water_depth_mask_filename = dem_name[:-4] + '_water_depth_mask_filtered.txt'
    save_water_depth(os.path.join(path_water_depth, water_depth_mask_filename), headers, water_depth_filtered)
    
    
    # Compute total braided index
    total_braiding_index = compute_braiding_index(water_depth_filtered)
    print(f"Best Max Elevation: {best_max_elevation}")
    print(f"Mean Above Length for Best Max Elevation: {best_mean_above_length}")
    print(f"Water depth data saved to {os.path.join(path_water_depth, water_depth_filename)}")
    print(f"Total braiding index is {total_braiding_index}")
