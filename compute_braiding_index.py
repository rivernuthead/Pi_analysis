import os
import rasterio
from skimage import morphology
import numpy as np
import matplotlib.pyplot as plt
from skimage import morphology
import cv2

####### Functions ##################################################################################

def compute_braiding_index(array, step=100):
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
    
def fill_small_holes(matrix, avg_target_kernel, area_threshold, connectivity, value=None):
    
    # Convert matrix in bool matrix
    matrix_bool = np.where(matrix>0, 1, 0)
    matrix_bool = np.where(matrix<0, 1, matrix_bool)
    matrix_bool = np.array(matrix_bool, dtype='bool') # Convert in array of bool
    
    # Set target and average
    ker=np.ones((avg_target_kernel,avg_target_kernel), np.float32)/(avg_target_kernel**2)
    matrix_target = np.where(np.isnan(matrix), 0, matrix)
    
    if value == None:
        matrix_target = cv2.filter2D(src=matrix_target,ddepth=-1, kernel=ker)
    elif value != None:
        matrix_target = np.ones((matrix.shape))*value
    

    # Perform morphological analysis
    matrix_out = morphology.remove_small_holes(matrix_bool, area_threshold=area_threshold, connectivity=connectivity)
    
    # Apply the target where holes have to be filled
    matrix_out = np.where(matrix_bool!=matrix_out, matrix_target, matrix)
    
    # matrix_out = matrix_target*matrix_out
    
    return matrix_out, matrix_target


def filter_water_depth_mask(array, y_res, min_width_ones=3, min_width_zeros=2):
    """
    Filter the water depth mask to keep only channels with wetted width > 3 mm.

    Parameters:
    mask (np.ndarray): Water depth mask.
    """
    y, x = array.shape
    y_res_mm = y_res * 1000
    min_width_ones = min_width_ones // y_res_mm  # transform mm in pixels
    min_width_zeros = min_width_zeros // y_res_mm  # transform mm in pixels
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
    
def read_tif(file_path):
    with rasterio.open(file_path) as src:
        array = src.read(1)  # Read the first band
    return array

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

####### Script ##################################################################################

# Options 
Active_braiding_index = 0
Total_braiding_index = 1

# Parameters
folder_home = os.getcwd()

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
    filename = 'q05_1_DEM9_hw.tif'
    print('*****************')
    print(filename)
    print('*****************')
    path = os.path.join(folder_home,  'hydraulic_analysis', 'output_data','water_depth_maps')
    target_wetted_width = 0.51 # in m
    x_res = 0.050  # Resolution along the x-axis in m
    y_res = 0.005  # Resolution along the y-axis in m
    cs_step = 50  #number of pixels between each cross section to compute braiding index,  10px=5cm (along dx) and 1px=0.5cm (along dy) 
    
    # Load water depth map computed from the hydraulic model
    water_depth = read_tif(os.path.join(path,filename))
    dy, dx = np.shape(water_depth)

    ##### 1. Select the optimal threshold procedure ######################################################
    # Loop over water depth threshold to fit the targetted wetted width (target_wetted_width)
    WD_threshold_range = np.arange(np.min(water_depth), np.max(water_depth), 2e-4)

    wetted_width_thrs_list = []
    for t in WD_threshold_range:
        water_depth_thrs = np.where(water_depth>=t, water_depth, 0)
        water_depth_thrs_area = len(water_depth_thrs[water_depth_thrs>0])
        wetted_width_thrs =  water_depth_thrs_area/dx * y_res
        wetted_width_thrs_list = np.append(wetted_width_thrs_list, wetted_width_thrs)
        
    # Find the index of the element closest to target_wetted_width
    closest_index = min(range(len(wetted_width_thrs_list)), key=lambda v: abs(wetted_width_thrs_list[v] - target_wetted_width))
    wetted_width = wetted_width_thrs_list[closest_index]
    WD_thrs = WD_threshold_range[closest_index]
    print(f'The optimal threshold is: {WD_thrs} m')
    print(f'The wetted width is: {wetted_width} m')
    ######################################################################################################

    ##### 2. Filtering procedure ######################################################
    # 2.1 Filter the water_depth mask with the optimal threshold
    water_depth_thrs = np.where(water_depth >= WD_thrs, water_depth, 0) 
    water_depth_mask = np.where(water_depth_thrs >0, 1, 0)
    water_depth_mask_copy = np.copy(water_depth_mask)

    water_depth_rsm_mask = morphology.remove_small_objects(water_depth_mask_copy>0, min_size=200, connectivity=1)
    water_depth_rsm = water_depth_mask*water_depth_rsm_mask
    
    # 2.2 FILL SMALL HOLES
    water_depth_rsm_hls= morphology.remove_small_holes(water_depth_rsm, area_threshold=100, connectivity=1)
    water_depth_rsm_hls_copy = np.copy(water_depth_rsm_hls)

    # 2.3 clean edges between channels and bars
    # Step 1: Convert the image to 0 and 255 (required for cv2.erode)
    water_depth_rsm_hls_binary = water_depth_rsm_hls_copy * 255
    water_depth_rsm_hls_binary = water_depth_rsm_hls_binary.astype(np.uint8)
    # Step 2: Define the kernel (structuring element) for erosion
    kernel = np.ones((3, 3), np.uint8)  # A 5x5 square kernel
    # Step 3: Perform opening operation 
    water_depth_opened = cv2.morphologyEx(water_depth_rsm_hls_binary, cv2.MORPH_OPEN, kernel)#cv2.erode(water_depth_rsm_hls_binary, kernel, iterations=1)

    #2.4 apply connected component to label the objects and keep only those that represents channels
    # Step 1:
    num_labels, labels_im, stats, centroids = cv2.connectedComponentsWithStats(water_depth_opened, connectivity=4)
    # Step 3: Set a minimum area threshold for filtering
    min_area = 2000  # Define your area threshold in pixels, 10px=5cm (along dx) and 1px=0.5cm (along dy) 

    # Step 4: Create a new binary image to hold only the large components
    water_depth_filtered = np.zeros_like(water_depth_opened)

    # Step 5: Loop through each component and filter based on area
    for i in range(1, num_labels):  # Start from 1 to skip the background (label 0)
        area = stats[i, cv2.CC_STAT_AREA]
        if area >= min_area:
            # Keep the component if it meets the area threshold
            water_depth_filtered[labels_im == i] = 255
            
    # Step 6: Convert the eroded image back to 0 and 1 if needed
    water_depth_filtered = water_depth_filtered // 255  # Dividing by 255 to get 0 and 1 back
    water_depth_filtered_copy = np.copy(water_depth_filtered)
    
    # 2.5 Filter out channels width < min_width_ones along each cross section
    water_depth_map_final = filter_water_depth_mask(water_depth_filtered_copy, y_res,min_width_ones=40,min_width_zeros=20)

    ######################################################################################################

    
    # Option: plot the filtered water depth map
    fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7,1,figsize=(14,14))
    aspect= 1
    ax1.imshow(water_depth_mask, cmap = 'Blues', aspect=aspect)
    ax2.imshow(water_depth_rsm, cmap = 'Blues', aspect=aspect)
    ax3.imshow(water_depth_rsm_hls, cmap = 'Blues', aspect=aspect)
    ax4.imshow(water_depth_opened, cmap = 'Blues', aspect=aspect)
    # Step 3: Normalize the labels image for visualization
    # Normalize the labels to the range 0-255 for visualization purposes
    labels_normalized = (labels_im.astype(np.float32) / num_labels) * 255
    labels_normalized = labels_normalized.astype(np.uint8)

    # Step 4: Apply a color map to the normalized labels image
    colored_labels = cv2.applyColorMap(labels_normalized, cv2.COLORMAP_TURBO)
    ax5.imshow(colored_labels, aspect=aspect)
    ax6.imshow(water_depth_filtered, cmap = 'Blues', aspect=aspect)
    ax7.imshow(water_depth_map_final, cmap = 'Blues', aspect=aspect)
    # Vertical cross-sections (for columns)
    for col in range(0, dx, cs_step):
        ax7.axvline(x=col, color='k', linestyle='--', linewidth=1)
    ax1.set_title('Water depth map thresholded (mask)')
    ax2.set_title('1. Water depth map rsm (mask)')
    ax3.set_title('2. Water depth map remove holes (mask)')
    ax4.set_title('3. Water depth map cleaned (opening operation)')
    ax5.set_title('4a. Water depth map labeled')
    ax6.set_title('4b. Water depth map filtered')
    ax7.set_title('5. Water depth map (final)')
    plt.tight_layout()
    plt.show()
    
    # Compute total braided index
    total_braiding_index = compute_braiding_index(water_depth_filtered, cs_step)
    print(f"Total braiding index is {total_braiding_index}")
