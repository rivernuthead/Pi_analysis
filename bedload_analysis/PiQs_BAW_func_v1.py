#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 17:52:03 2023

@author: erri
"""
import os
import cv2
import numpy as np
import matplotlib.pyplot as plt
# import math as m
import PyPDF2
from PIL import Image
from skimage import morphology
from scipy.signal import convolve
from scipy.interpolate import interp1d
from scipy import ndimage

def detect_blurry_regions(image_path, threshold=100):
    # Load the image
    image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)

    # Apply Laplacian operator to calculate gradient magnitude
    gradient = cv2.Laplacian(image, cv2.CV_64F)
    gradient_magnitude = np.absolute(gradient)

    # Normalize the gradient magnitude
    gradient_magnitude = (gradient_magnitude - np.min(gradient_magnitude)) / (np.max(gradient_magnitude) - np.min(gradient_magnitude))
    gradient_magnitude = (gradient_magnitude * 255).astype(np.uint8)

    # # Apply threshold to identify blurry regions
    # blurry_regions = gradient_magnitude < threshold
    
    return gradient_magnitude





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

def neigh_filter(matrix, k_dim, param):
    matrix_out=np.zeros((matrix.shape))
    # domain = matrix_dom = np.pad(matrix, 1, mode='edge') # Create neighbourhood analysis domain
    
    for i in range(0,matrix.shape[0]):
        for j in range(0,matrix.shape[1]):
            y_inf, y_sup = int(i - (k_dim-1)/2), int(i + (k_dim-1)/2)
            x_inf, x_sup = int(j - (k_dim-1)/2), int(j + (k_dim-1)/2)
            # print('rows: ', y_inf, y_sup)
            # print('cols: ', x_inf, x_sup)
            
            if y_inf<0:
                y_inf=0
            if y_sup>matrix.shape[0]:
                y_sup=matrix.shape[0]
            
            if x_inf<0:
                x_inf=0
            if x_sup>matrix.shape[1]:
                x_sup=matrix.shape[1]

            kernel = matrix[y_inf:y_sup+1, x_inf:x_sup+1]       
            
            # print(kernel)

            if param=='variance':
                matrix_out[i,j] = np.var(kernel)
            elif param=='stdev':
                matrix_out[i,j] = np.std(kernel)
    return matrix_out, kernel

def moving_average_filter(x, window_size):
    """
    Applies a moving average filter to the input array x with the given window size.
    """
    # Define the filter coefficients
    coefficients = np.ones(window_size) / window_size
    
    # Apply the filter using convolution
    y = convolve(x, coefficients, mode='valid')
    
    return y

def remove_grouped_cells(matrix, kernel_dim, threshold):
    '''
    INPUT:
        matrix
        kernel_dim
        threshold
    '''
    padding = ((1, 1), (1, 1)) # Padding amount for each side (top, bottom, left, right)
    domain = np.pad(matrix, padding, mode='constant', constant_values=0) # Create the domain
    output_matrix = np.copy(domain) # Create the output matrix as a copy of the domain
    m, n = domain.shape # Set domain dimension
    k_dim_x, k_dim_y = kernel_dim

    for i in range(m - k_dim_x + 1):
        for j in range(n - k_dim_y + 1):
            region = domain[i:i+k_dim_x, j:j+k_dim_y]
            # print(region)
            
            # Compute the number of zeros along the edge
            edge_zeros = 0
            edge_zeros += np.sum(region[0, :] == 0) + np.sum(region[-1, :] == 0)  # Top and bottom edges
            edge_zeros += np.sum(region[1:-1, 0] == 0) + np.sum(region[1:-1, -1] == 0) # Left and right edges (excluding corners to avoid double counting)

            # So, if the region is surrounded by zeros and the number of ones inside is less than the threshold
            if (edge_zeros == 2*region.shape[0]+2*(region.shape[1]-2)) and (np.nansum(region) <= threshold):
                output_matrix[i:i+k_dim_x, j:j+k_dim_y] = 0
                # print('test')
            else:
                pass

    output_matrix = output_matrix[1:-1, 1:-1] # remove the pad

    return output_matrix



    
    
def resample_array(x, y):
    """
    Resamples the input array y to fit the same size as the input array x using linear interpolation.
    """
    # Create a linear interpolation function
    f = interp1d(np.arange(len(y)), y, kind='linear')
    
    # Resample the function at the x values
    y_resampled = f(np.linspace(0, len(y)-1, len(x)))
    
    return y_resampled

def interpolate_nans(arr):
    """
    Fill NaN values in a 1-dimensional numpy array with linear interpolation.
    """
    nan_mask = np.isnan(arr)
    x = np.arange(len(arr))
    arr[nan_mask] = np.interp(x[nan_mask], x[~nan_mask], arr[~nan_mask])
    
    return arr


def downsample_matrix_interpolation(matrix, factor):
    # Get the shape of the original matrix
    height, width = matrix.shape

    # Calculate the new dimensions after downsampling
    new_height = height // factor
    new_width = width // factor

    # Create a grid of coordinates for the new downsampling points
    row_indices = np.linspace(0, height - 1, new_height)
    col_indices = np.linspace(0, width - 1, new_width)

    # Perform bilinear interpolation to estimate the values at new points
    downsampled_matrix = ndimage.map_coordinates(matrix, 
                                                 np.meshgrid(row_indices, col_indices),
                                                 order=1,
                                                 mode='nearest')

    # Reshape the downsampled matrix to the new dimensions
    downsampled_matrix = downsampled_matrix.reshape(new_height, new_width)

    return downsampled_matrix

def non_overlapping_average(image, kernel_size):
    # Get the shape of the input image
    height, width = image.shape

    if isinstance(height / kernel_size, int) or isinstance(width / kernel_size, int):
        print("Warning: kernel size does not fit the image size: the function will ignore the remaining pixels at the right and bottom edges")
    # else:
    #     print(f"{number} is not an integer.")

    # Calculate the new dimensions for the non-overlapping blocks
    new_height = height // kernel_size
    new_width = width // kernel_size

    # Reshape the image into non-overlapping blocks
    blocks = image[:new_height * kernel_size, :new_width * kernel_size].reshape(new_height, kernel_size, new_width, kernel_size)

    # Calculate the average within each block
    block_averages = blocks.mean(axis=(1, 3))

    return block_averages


def update_matrix_dimensions(matrix1, matrix2):
    # Calculate the maximum dimensions
    max_rows = max(matrix1.shape[0], matrix2.shape[0])
    max_cols = max(matrix1.shape[1], matrix2.shape[1])

    # Create new matrices with NaN values and the new dimensions
    updated_matrix1 = np.full((max_rows, max_cols), np.nan)
    updated_matrix2 = np.full((max_rows, max_cols), np.nan)

    # Copy the original data to the new matrices
    updated_matrix1[:matrix1.shape[0], :matrix1.shape[1]] = matrix1
    updated_matrix2[:matrix2.shape[0], :matrix2.shape[1]] = matrix2

    return updated_matrix1, updated_matrix2


def cut_matrices_to_minimum_dimension(matrix1, matrix2):
    # Calculate the minimum dimensions (number of columns and rows)
    min_cols = min(matrix1.shape[1], matrix2.shape[1])
    min_rows = min(matrix1.shape[0], matrix2.shape[0])
    
    # Cut both matrices to match the new dimensions
    matrix1_cut = matrix1[:min_rows, :min_cols]
    matrix2_cut = matrix2[:min_rows, :min_cols]

    
    return matrix1_cut, matrix2_cut
# ---------------------------------------------------------------------------- #




# ---------------------------------------------------------------------------- #
#                           ACTIVITY DURATION PERIODS                          #
# ---------------------------------------------------------------------------- #

# UNCOMMENT FOR TEST:-------------------------------------------------------#

# stack_bool = np.array([0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0,1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0])  # Test array
# # stack_bool = np.array([0,0,0,0,0,-1,-1,0,0,0,1,0,0,0,1,1,1,0,0,0,-1,-1,-1,1,1,1,0,1,0,0,0]) # Test array
# # stack_bool = np.array([1, -1,  0,  0,  0,  1,  0,  0,  0]) # Test array
# # stack_bool = np.array([0, 0,  1, -1,  0,  0,  0,  -1,  0,  0]) # Test array
# # stack_bool = np.array([1, 1,  1, 1,  1,  1,  1,  1,  1,  1]) # Test array

# stack_bool = np.expand_dims(stack_bool, axis=1)
# stack_bool = np.expand_dims(stack_bool, axis=1)

# analysis_list = ['Consecutive_numbers', 'PiQs_bool_maps', 'Trim_1_switch_only_pixel']
# analysis_list = ['Consecutive_numbers', 'PiQs_bool_maps']

#-----------------------------------------------------------------------------#

def activity_duration_periods(stack_bool, analysis_list):
    '''

    Parameters
    ----------
    stack_bool : TYPE
        DESCRIPTION.
    analysis_list : TYPE
        DESCRIPTION.

    Returns
    -------
    if PiQs_bool_maps in analysis list: the input file is a boolean active list
    where 1 means activity and 0 means inactivity.
    
    consecutive_ones_stack: returns the length of the activity periods as the
        number of time a pixel shows activity.
        Data are interposed by a number of zeros that corresponds to the number
        of time a cell remains inactive.
    consecutive_zeros_stack: returns the length of the inactivity periods as
        the number of time a pixel shows inactivity.
    single_ones_stack: returns an array of zeros. Where a 1 appear, this
        identify the position of a single (isolated) one - isolated periods of
        activity (flashes)
    

    '''
    def find_consecutive_number_lengths(array, number):
        consecutive_number_lengths = []
        current_length = 0

        for num in array:
            if num != number:
                consecutive_number_lengths.append(0)
            if num == number:
                current_length += 1
            elif current_length > 0:
                consecutive_number_lengths.append(current_length)
                current_length = 0

        if current_length > 0:
            consecutive_number_lengths.append(current_length)

        # consecutive_number_lengths = consecutive_number_lengths[0:]
        return consecutive_number_lengths

    # INITIALIZE STACK AND ARRAY
    # activation time stack contains the time between switches. The first layer of this stack contains the first sctivation time that is a lower limit in time because we ignore how long the pixel has keept the same nature in the past.
    act_time_stack = np.zeros(stack_bool.shape)
    # This is the 2D matrix that collect the number of switch over time
    switch_matrix = np.zeros(stack_bool.shape[1:])

    consecutive_ones_stack = np.zeros(stack_bool.shape)
    consecutive_zeros_stack = np.zeros(stack_bool.shape)
    consecutive_minus_ones_stack = np.zeros(stack_bool.shape)
    single_ones_stack = np.zeros(stack_bool.shape)

    # This is a matrix that collect all the pixel locations that have never been active
    never_active_matrix = np.zeros(stack_bool.shape[1:])

    dim_t, dim_y, dim_x = stack_bool.shape

    # analysis_list = ['switch_number', 'consecutive_numbers']

    for x in range(0, dim_x):
        for y in range(0, dim_y):
            # Slice the stack in a single pixel array where data is collected over time
            slice_array = stack_bool[:, y, x]

            # PRINT FOR TEST
            # print('Trimmed slice: ', slice_array)

            if 'Consecutive_numbers' in analysis_list:
                raw_slice_array = stack_bool[:, y, x]
                consecutive_ones_array = find_consecutive_number_lengths(
                    raw_slice_array, 1)
                consecutive_zeros_array = find_consecutive_number_lengths(
                    raw_slice_array, 0)
                consecutive_minus_ones_array = find_consecutive_number_lengths(
                    raw_slice_array, -1)

                # Compute where single activity instats are within zeros (activity flash)
                # Find the indices where ones are surrounded by zeros
                indices = np.where((raw_slice_array == 1) & (
                    np.roll(raw_slice_array, 1) == 0) & (np.roll(raw_slice_array, -1) == 0))

                # Create a new array with only the isolated ones
                isolated_ones_array = np.zeros_like(raw_slice_array)
                isolated_ones_array[indices] = 1

            # This is the number of zero values before the first non-zero value in the sliced array. It is initialized to zero.
            n_zero = 0
            time_array = []
            '''Check if the sliced array has np.nan. If so, fill this array with
            np.nan and then fill the act_time_stack. This avoids to include in the
            computation pixel that are near the domain edge or pixels that during
            the runs fail in signal detection'''
            if np.isnan(slice_array).any():  # check if a has np.nan value, if so fill matrix with np.nan
                switch_matrix[y, x] = np.nan
                act_time_stack[:, y, x] = np.nan
                pass

            ''' This part treats the cases in which the first entry of the sliced
                array is zero.
            1. If the entire array is full of zeros means that in that position
                nothing occurs and the active periods array will be filled with
                np.nan
            2. If the first entry is zero and there are non-zero value in the array
                the script count the number of zeros before the first non-zero
                value and then trim from the slice array the zeros before the first
                non-zero value. Now the resized slice array is ready to go on to
                the next section.
            '''
            if slice_array[0] == 0:
                if np.all(slice_array == 0):  # Check if the sliced array is full of zeros
                    # Fill the array with np.nan to keep them transparent
                    time_array = np.array([np.nan])
                    never_active_matrix[y, x] = 1
                else:
                    # Number of zero values before the first non-zero value
                    n_zero = np.array(np.where(slice_array != 0))[0, 0]
                    if 'Trim_zero_values_start_period' in analysis_list:
                        # Trim the zero values before the first non-zero value.
                        slice_array = slice_array[n_zero:]

            ''' This part treats the cases in which the first entry of the sliced
                array is a non-zero value.
                1. The counter is then initialize to one and the first entry sign
                    is detected as the target sign
                2. The script takes all the adjacent elements of the slice array
                    and detect the sign of each one.
                    a. If the two adjacent elements have the same sign, the counter
                        is updated with +1
                    b. Elif if both the second element is zero and the first or the
                    second element have the target sign the counter is +1
                    c. Elif a switch occurs, so there is a change in sign in the
                        two adjacent elements, or the second element shows a sign
                        that is different from the target sign the counter values
                        is stored in the time array, the target sign change and the
                        count is updated to 1.
                3. The script create the time_array trimming the last period
                    because it is a lover boundary as well as the zeros before the
                    first non-zero entry.
                    
            '''
            if slice_array[0] != 0:  # If the first entry of the sliced array is non-zero
                period_count = 1  # Initialize the count variable. This variable will count the number of activation instants
                switch_count = 0
                # This variable collects the sign of the first element of each same-nature period
                target_sign = np.sign(slice_array[0])
                for i in range(0, len(slice_array)-1):  # Loop over the sliced array
                    # a1 and a2 are the two adjacent element in the sliced array
                    a1, a2 = slice_array[i], slice_array[i+1]

                    # If two consecutive elements have the same naure
                    if np.sign(a1) == np.sign(a2):
                        period_count += 1  # If two consecutive elements have the same nature the count increases

                    elif np.sign(a1)*np.sign(a2) == 0 and (np.sign(a2) == target_sign or np.sign(a1) == target_sign):
                        period_count += 1  # The count increases also if one or both elements are zero but the non-zero value has the target sign

                    # The count stops when a switch occours or when one of the two elements shows a sign different from the target sign
                    elif np.sign(a1) != np.sign(a2) and (np.sign(a2) != target_sign or np.sign(a1) != target_sign):
                        # This operation append to time_array the count value with his sign. This could be useful to keep trace of the nature of the period.
                        time_array = np.append(
                            time_array, period_count*target_sign)
                        target_sign = -1*target_sign  # Update the target sign
                        period_count = 1  # Update the count variable that will starts again from zero
                        switch_count += 1
                        pass
                # FILL THE SWITCH MATRIX
                switch_matrix[y, x] = switch_count

                # By now the last period is not calculated (actually because, as the first one, it is only a lower boundary of time because it doesn't appear within two switches) so this operation appends this value manually
                time_array = np.append(
                    time_array, (len(slice_array)-np.sum(np.abs(time_array)))*target_sign)
                # time_array[0] = time_array[0] + np.sign(time_array[0])*n_zero # Ths operation append, if present, the number of zeroes before the first non-zero value calculated on the very first sliced array (n_zero variable)
            if 'Trim_1st_period' in analysis_list:
                # TRIM THE LAST PERIOD BECAUSE IT IS NOT RELIABLE
                # This number correspond to the index of the last period in the time_array that is not reliable (as the first one)
                ind = np.max(np.where(time_array != 0))
                # So in the filling process I want to exclude the last period:
                # This operation fills the stack with time_array
                act_time_stack[:ind, y, x] = time_array[:ind]

            # SLICE ARRAY WITH NO SWITCH LEADS TO A TIME ARRAY WITH NO PERIODS
            if len(time_array) == 0 or len(time_array) == 1:  # If sliced array does not contain any switch (so if the length of the time array is 0 in the case we do not consider the last period - see above - or if the length is 1 so only one period is considered - the last one - this operation fills the output matrix with np.nan)
                # switch_matrix[y,x] = np.nan # Fill switch matrix with np.nan
                # Fill activation time stack with np.nan
                for t in range(0, act_time_stack.shape[0]):
                    act_time_stack[t, y, x] = np.nan

                ''' This section is introduced to eliminate to the dataset pixels that
                show only one switch and so where only the first period is considered'''
            if 'Trim_1_switch_only_pixel' in analysis_list:
                if len(time_array) == 2:  # If only one switch occours, two periods were identified. This two periods are not completely reliable since they were not observed between two switches
                    # switch_matrix[y,x] = np.nan # Fill switch matrix with np.nan
                    # Fill activation time stack with np.nan
                    for t in range(0, act_time_stack.shape[0]):
                        act_time_stack[t, y, x] = np.nan

            else:  # If at least one (or two) switches occours
                # Fill switch matrix
                # switch_matrix[y,x] = len(time_array[time_array!=0]) # To provide the number of switch we decide to trim the zero values to keep only the number of switch

                # Fill activation time stack
                # act_time_stack[:len(time_array),y,x]=time_array # To provide the time between each detected switch
                # This number correspond to the index of the last period in the time_array that is not reliable (as the first one)
                ind = np.max(np.where(time_array != 0))
                # So in the filling process I want to exclude the last period:
                # To provide the time between each detected switch
                act_time_stack[:ind, y, x] = time_array[:ind]

            # FILL consecutive_ones_stack
            consecutive_ones_stack[:len(
                consecutive_ones_array), y, x] = consecutive_ones_array

            # FILL consecutive_zeros_stack
            consecutive_zeros_stack[:len(
                consecutive_zeros_array), y, x] = consecutive_zeros_array

            # FILL consecutive_minus_ones_stack
            consecutive_minus_ones_stack[:len(
                consecutive_minus_ones_array), y, x] = consecutive_minus_ones_array

            # FILL stack that identify where single ones appear
            single_ones_stack[:len(isolated_ones_array),
                              y, x] = isolated_ones_array

            # # FILL THE SWITCH MATRIX
            # switch_matrix[y,x] = switch_count

    switch_dist_stack = act_time_stack[:, :, :]
    consecutive_ones_stack = consecutive_ones_stack[:, :, :]
    consecutive_zeros_stack = consecutive_zeros_stack[:, :, :]

    if 'PiQs_bool_maps' in analysis_list:
        return consecutive_ones_stack, consecutive_zeros_stack, single_ones_stack
    else:
        return act_time_stack, switch_dist_stack, consecutive_ones_stack, consecutive_zeros_stack, consecutive_minus_ones_stack


# UNCOMMENT FOR TEST:-------------------------------------------------------#
# consecutive_ones_stack, consecutive_zeros_stack, single_ones_stack = activity_duration_periods(stack_bool, analysis_list)

# print(consecutive_ones_stack) # Uncomment for test
# # print(consecutive_zeros_stack) # Uncomment for test
# # print(single_ones_stack) # Uncomment for test

# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
#                                 REMOVE FILES                                 #
# ---------------------------------------------------------------------------- #

def remove_file(file_path):
    try:
        os.remove(file_path)
        print("File "+file_path+" has been removed.")
    except FileNotFoundError:
        print("File "+file_path+" not found.")
    except PermissionError:
        print("Permission denied: Unable to remove file" +file_path+".")
    except Exception as e:
        print("An error occurred while removing the file"+file_path+": {e}")
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
#                             PDF MATPLOTLIB REPORT                            #
# ---------------------------------------------------------------------------- #

def pdf_matplotlib_report(file_path, report_path, i):
    '''
    This function allows PDF report creation in loop of plots
    
    INPUT:
        file_path: is the path in which the plot save the current plot in PDF
        report_path: is the path in which the function save the plot report in PDF
        i: the loop iteration (starts from 0)
    Note: To create different report for different runs, insert the run variable in the
        file_path and report_path definition
    '''

    merger = PyPDF2.PdfMerger()

    # At the very first iteration, if a report file already exists, remove it.
    if i==0 and os.path.exists(report_path):
        remove_file(report_path)
    
    # If a report file do not already exists, print the first page of the pdf report
    if not(os.path.exists(report_path)):
        plt.savefig(report_path, dpi=1200)

    # Then, go on with the append procedure:
    else:
        # Open and append the existing PDF
        with open(report_path, "rb") as existing_file:
            merger.append(existing_file)

        # Open and append the new PDF chart
        with open(file_path, "rb") as chart_file:
            merger.append(chart_file)

        # Save the merged PDF
        with open(report_path, "wb") as merged_file:
            merger.write(merged_file)

    '''EXAMPLE
    home_dir = os.getcwd() # Get the home directory
    report_path = os.path.join(home_dir, 'report_TEST.pdf') # Define the report path


    for i in range(10):
        
        # Create a random matrix
        matrix = np.random.rand(10, 10)
        
        # Create the plot
        plt.imshow(matrix, cmap='viridis', origin='upper', interpolation='nearest')
        plt.colorbar()  # Optional: Adds a colorbar to the plot.
        plt.title('Iteration: '+ str(i))
        plt.xlabel('Columns')
        plt.ylabel('Rows')
        
        file_path = os.path.join(home_dir, 'img_TEST.pdf') # Define the file path
        plt.savefig(file_path, dpi=1200)
        
        # Create PDF report
        pdf_matplotlib_report(file_path, report_path, i)
        
        
        plt.show()
    '''
# ---------------------------------------------------------------------------- #




# ---------------------------------------------------------------------------- #
#                            PDF REPORT FROM IMAGES                            #
# ---------------------------------------------------------------------------- #
import os
from PIL import Image
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas

def create_pdf_report(folder_path, common_name, output_pdf):
    # Get a list of image files in the folder with the common name
    image_files = [file for file in os.listdir(folder_path) if common_name in file]
    
    # Sort the image files by name
    image_files.sort()
    
    # Create a PDF report
    c = canvas.Canvas(output_pdf, pagesize=A4)
    
    for image_file in image_files:
        # Add each image to the PDF
        image_path = os.path.join(folder_path, image_file)
        img = Image.open(image_path)
        c.drawImage(image_path, 50, 50, width=img.width, height=img.height)
        c.showPage()
    
    c.save()

# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
#                     SORT NAMES WITH FIRST PART AS NUMBER                     #
# ---------------------------------------------------------------------------- #
# Define a custom sorting key function
def custom_sort_key(item):
    # Split the string at the first underscore ('_') character
    parts = item.split('_')

    try:
        # Attempt to convert the part before the underscore to an integer
        num = int(parts[0])
        return num
    except ValueError:
        # If the conversion fails, return a large number to push it to the end
        return float('inf')



# ---------------------------------------------------------------------------- #
#                            MORPHOLOGICAL ANALYSIS                            #
# ---------------------------------------------------------------------------- #
def morphological_three_class_analysis(matrix, min_size=4, conn=0):
    '''

    '''
    mask_remove_isolated_ones      = morphology.remove_small_objects(np.array(matrix == +1, dtype='bool'), min_size=min_size, connectivity=conn)
    mask_remove_isolated_minusones = morphology.remove_small_objects(np.array(matrix == -1, dtype='bool'), min_size=min_size, connectivity=conn)
    matrix_cld                     = matrix*mask_remove_isolated_ones + matrix*mask_remove_isolated_minusones

    # REMOVE LONGITUDINAL SLICES OF ACTIVE GROUPED PIXELS
    kernel_dim = np.array([7, 11, 15, 19, 21, 25, 29]) # Array of different kernel longitudinal length

    matrix_cld2 = remove_grouped_cells(matrix_cld, (3, 3), 1) # First iteration (length = 3, threshold = 1)

    for r, s in zip(kernel_dim, kernel_dim - 1): # Loop with different lengths and thresholds
        matrix_cld2 = remove_grouped_cells(matrix_cld2, (3, r), s)






def spatial_temporal_activation_filtering(stack, ker_dim, thrs):
    
    '''
    '''
    
    stack_out = np.copy(stack)
    
    
    
    # Define the index that identifys the kernel centered cell
    t_centered_index = int(np.floor(ker_dim[0]/2))
    x_centered_index = int(np.floor(ker_dim[1]/2))
    y_centered_index = int(np.floor(ker_dim[1]/2))
    
    # Pad the stack with a one-thick layer
    stack_padded = np.pad(stack, ((t_centered_index, t_centered_index), (y_centered_index, y_centered_index), (x_centered_index, x_centered_index)), mode='edge')
    
    for i in range(0, stack.shape[1]): # y axes
        for j in range(0, stack.shape[2]): # x axes
            for t in range(0, stack.shape[0]): # time
            
                neighborhood = stack_padded[t:t+ker_dim[0], i:i+ker_dim[1], j:j+ker_dim[1]] # Define the kernel
                
                neigh_dim = neighborhood.shape[0]*neighborhood.shape[1]*neighborhood.shape[2]
                
                
                if neighborhood[t_centered_index,y_centered_index,x_centered_index]==1 and np.nansum(neighborhood[:,y_centered_index,x_centered_index])==1: # 0,1,0 over time
                    # print('Slice: ', neighborhood[:,y_centered_index,x_centered_index])
                    neigh_sum = np.nansum(neighborhood==1) # Compute the number of ones
                    
                    # print('# of ones: ', neigh_sum)
                    # print(t,i,j)
                    
                    if neigh_sum < thrs*(neigh_dim):
                        stack_out[t,i,j] = 0
                        # print(neighborhood)
                        
                        # print('Out value', stack_out[t,i,j])
                        
                        # print()
                        
                    else:
                        pass
                    
                elif neighborhood[t_centered_index,y_centered_index,x_centered_index]==0 and np.nansum(neighborhood[:,y_centered_index,x_centered_index])==2: # 1,0,1 over time
                    # print('Slice: ', neighborhood[:,y_centered_index,x_centered_index])
                    neigh_sum = np.nansum(neighborhood==0) # Compute the number of zeros
                    
                    # print('# of zeros: ', neigh_sum)
                    # print(t,i,j)
                    
                    if neigh_sum < thrs*(neigh_dim):
                        stack_out[t,i,j] = 1
                        # print(neighborhood)
                        
                        # print(stack_out[t,i,j])
                        
                        # print()
                    else:
                        pass
                    
                    
    return stack_out