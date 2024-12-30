#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 11:38:09 2023

@author: erri

This script compute the DoD envelope.

INPUT:
    the script takes in input a stack where all the DoDs are stored in a structure as shown below:
        
    DoD input stack structure:
        
        DoD_stack[time,y,x,delta]
        DoD_stack_bool[time,y,x,delta]
        
         - - - 0 - - - - 1 - - - - 2 - - - - 3 - - - - 4 - - - - 5 - - - - 6 - - - - 7 - - - - 8 - -  >    delta
      0  |  DoD 1-0   DoD 2-0   DoD 3-0   DoD 4-0   DoD 5-0   DoD 6-0   DoD 7-0   DoD 8-0   DoD 9-0
      1  |  DoD 2-1   DoD 3-1   DoD 4-1   DoD 5-1   DoD 6-1   DoD 7-1   DoD 8-1   DoD 9-1
      2  |  DoD 3-2   DoD 4-2   DoD 5-2   DoD 6-2   DoD 7-2   DoD 8-2   DoD 9-2
      3  |  DoD 4-3   DoD 5-3   DoD 6-3   DoD 7-3   DoD 8-3   DoD 9-3
      4  |  DoD 5-4   DoD 6-4   DoD 7-4   DoD 8-4   DoD 9-4
      5  |  DoD 6-5   DoD 7-5   DoD 8-5   DoD 9-5
      6  |  DoD 7-6   DoD 8-6   DoD 9-6
      7  |  DoD 8-7   DoD 9-7
      8  |  DoD 9-8
         |
         v
        
         time

OUTPUT:
    
"""

# =============================================================================
# IMPORT LIBRERIES
# =============================================================================
import os
import numpy as np

# =============================================================================
# FUNCTIONS
# =============================================================================

def keep_every_j_elements(arr, j):
    if j <= 0:
        raise ValueError("J should be a positive integer")

    result = arr[::j,:,:]
    return result

def compute_stats_over_stack_groups(stack, t_scale):
    """
    Compute the mean value and standard deviation of the number of active pixels (value > 0)
    for each map, for the sum of all two adjacent maps without overlapping, for the sum of all
    three adjacent maps without overlapping, and so on until the sum of all maps.

    Parameters:
    stack (np.ndarray): The input stack with shape (40, 140, 1259).

    Returns:
    np.ndarray: A matrix where the first column is the group size, the second column is the mean
                number of active pixels, and the third column is the standard deviation.
    """
    num_maps = stack.shape[0]
    results = []

    # Iterate over different group sizes
    for group_size in range(1, num_maps + 1):
        active_pixels = []
        
        # Iterate over the stack with steps equal to group_size
        for start_idx in range(0, num_maps, group_size):
            end_idx = start_idx + group_size
            if end_idx > num_maps:
                break
            
            # Sum the maps in the current group
            group_sum = np.sum(stack[start_idx:end_idx], axis=0)

            # Compute the active width
            num_active_pixels = np.count_nonzero(group_sum > 0)
            active_pixels.append(num_active_pixels/120/stack.shape[2])
        
        # Calculate active width mean and standard deviation
        mean_active_pixels = np.mean(active_pixels)
        stdev_active_pixels = np.std(active_pixels)

        # Append the results as a row to the list
        results.append([group_size*t_scale, mean_active_pixels, stdev_active_pixels])
    
    # Convert the results list to a numpy array
    results_matrix = np.array(results)
    
    return results_matrix

def trim_nan_matrices(stack):
    """
    Trim the input stack by removing matrices that are completely full of np.nan values.

    Parameters:
    stack (np.ndarray): The input stack with shape (t, x, y).

    Returns:
    np.ndarray: A trimmed stack with matrices full of np.nan removed.
    """
    # Collect indices of matrices that are not completely full of np.nan
    valid_indices = [i for i in range(stack.shape[0]) if not np.all(np.isnan(stack[i]))]
    
    # Create a new stack with only the valid matrices
    trimmed_stack = stack[valid_indices]
    
    return trimmed_stack
# =============================================================================
# RUN PARAMETERS
# =============================================================================

set_names = ['q05_1', 'q07_1','q10_2', 'q10_3', 'q10_4', 'q15_2', 'q15_3', 'q20_2']
# set_names = ['q05_1']
# set_names = ['q05_1', 'q07_1','q10_2','q15_3', 'q20_2']


# =============================================================================
# LOOP OVER RUNS
# =============================================================================
home_dir = os.path.join(os.getcwd(), 'morphological_analysis') # Home directory

for set_name in set_names:
    print()
    print()
    print(set_name, ' is running...')
    # IMPORT DoD STACK AND DoD BOOL STACK
    output_dir = os.path.join(home_dir,'output_data', 'envelopes')
    if not(os.path.exists(output_dir)):
        os.mkdir(output_dir)
    DoDs_folder = os.path.join(home_dir, 'output_data', 'DoDs', 'DoDs_stack') # Input folder
    stack_name = 'DoD_stack' + '_' + set_name + '.npy' # Define stack name
    stack_bool_name = 'DoD_stack_bool_' + set_name + '.npy' # Define stack bool name
    stack_path = os.path.join(DoDs_folder,stack_name) # Define stack path
    stack_bool_path = os.path.join(DoDs_folder,stack_bool_name) # Define stack bool path
    
    stack = np.load(stack_path) # Load DoDs stack
    stack_bool = np.load(stack_bool_path) # Load DoDs boolean stack
    
    dim_t, dim_y, dim_x, dim_delta = stack.shape # Define time dimension, crosswise dimension and longitudinal dimension
    
    # Create activity map from the stack (1=active, 0=inactive)
    act_stack = np.abs(stack_bool)
    
    # =============================================================================
    # LOOP OVER DIFFERENT DoD TIMESPAN
    # =============================================================================
    
    for d in range(0,stack.shape[3]):
        act_stack_d = act_stack[:,:,:,d]
        act_stack_trimmed = trim_nan_matrices(act_stack_d)
        
        envMAA_report =  compute_stats_over_stack_groups(act_stack_trimmed, 1)
        
        # SAVE
        fmt = ['%d', '%.4f', '%.4f']
        fmt_string = ', '.join(fmt)
        header='Generated by: ' + os.path.basename(__file__) + '\nGroup size, mean, stdev'
        np.savetxt(os.path.join(output_dir, set_name + '_timespan'+ str(d) + '_envMAW_report.txt'), envMAA_report, delimiter=',', header=header, fmt=fmt_string)