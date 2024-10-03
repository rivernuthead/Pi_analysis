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

    
"""

import os
import numpy as np

def keep_every_j_elements(arr, j):
    if j <= 0:
        raise ValueError("J should be a positive integer")

    result = arr[::j,:,:]
    return result


runs = ['q07_1','q10_2', 'q10_3', 'q10_4', 'q15_2', 'q15_3', 'q20_2']
# runs = ['q07_1']
# runs = ['q05_1', 'q07_1','q10_2','q15_3', 'q20_2']

home_dir = os.path.join(os.getcwd(), 'morphological_analysis') # Home directory
output_folder = os.path.join(home_dir,'output_data', 'envelopes')
if not(os.path.exists(output_folder)):
    os.mkdir(output_folder)
for run in runs:
    print()
    print()
    print(run, ' is running...')
    
    # IMPORT DoD STACK AND DoD BOOL STACK
    DoDs_folder = os.path.join(home_dir, 'output_data', 'DoDs', 'DoDs_stack') # Input folder
    stack_name = 'DoD_stack' + '_' + run + '.npy' # Define stack name
    stack_bool_name = 'DoD_stack_bool_' + run + '.npy' # Define stack bool name
    stack_path = os.path.join(DoDs_folder,stack_name) # Define stack path
    stack_bool_path = os.path.join(DoDs_folder,stack_bool_name) # Define stack bool path
    
    MAV_env_0_timespan_list = []
    
    
    stack = np.load(stack_path) # Load DoDs stack
    stack_bool = np.load(stack_bool_path) # Load DoDs boolean stack
    
    dim_t, dim_y, dim_x, dim_delta = stack.shape # Define time dimension, crosswise dimension and longitudinal dimension
    
    # Create activity map from the stack (1=active, 0=inactive)
    act_stack = np.abs(stack_bool)
    
    # Open a file in write mode
    with open(os.path.join(output_folder, run + '_envelope_timespan0_starting_points.txt'), 'w') as file:
        header0 = run+'\nMAW envelopes for the finest timeaspan starting from different consecutive DoD\nGenerated by ' + str(os.path.basename(__file__))
        # Write the header to the file
        file.write(header0 + "\n")
    
        # PERFORM A LOOP FOR DIFFERENT DoDs TIMESPAN
        for timespan in range(0,int(stack_bool.shape[3]//2)):
            print('timespan: ' + str(timespan))
            # timespan = 1
                
            if timespan==0:
                sliced_stack = act_stack[:,:,:,timespan]
            else:
                sliced_stack = act_stack[:-timespan,:,:,timespan]
            
            
            if timespan==0:
                m=0
                for s in range(0, dim_t): # S is the starting point in envelope computing
                    
                    sliced_stack = act_stack[s:,:,:,timespan]
                
                    # Compute the array stack that contains the array on which perform the envelope:
                    env_sliced_stack = keep_every_j_elements(sliced_stack, timespan+1)
                    
                    act_stack_envelope = np.cumsum(env_sliced_stack, axis=0)
                    
                    if s==0:
                        # Save the envelope stack as numpy binary file
                        np.save(os.path.join(output_folder, run + '_envelope_timespan'+str(timespan)+'.npy'), act_stack_envelope)
                    
                    # Compute the morphological active width (MAW)
                    act_stack_envelope_bool = np.where(act_stack_envelope>0,1,act_stack_envelope)
                    MAW_envelope = np.nansum(act_stack_envelope_bool,axis=2)
                    MAW_envelope = np.nansum(MAW_envelope,axis=1)
                    report_MAW_envelope = MAW_envelope/dim_x/120
                    
                    # Create the row string with the starting point counter
                    row_string = f"Starting point {s} ," + ", ".join(f"{num:.4f}" for num in report_MAW_envelope)
                    # row_string = ", ".join(f"{num:.4f}" for num in report_MAW_envelope)
                    
                    # # Save report_MAW_envelope as row table
                    MAV_env_0_timespan_list.append(report_MAW_envelope)
                    
                    # Write the row to the file
                    file.write(row_string + "\n")
                    
                    
            elif timespan>0: # Perform all the possible combination give the same timespan
                for m in range(0,timespan+1):
                    sliced_stack = sliced_stack[m:,:,:]
                    
                    # Compute the array stack that contains the array on which perform the envelope:
                    env_sliced_stack = keep_every_j_elements(sliced_stack, timespan+1)
                    
                    act_stack_envelope = np.cumsum(env_sliced_stack, axis=0)
                    
                    # Save the envelope stack as numpy binary file
                    np.save(os.path.join(output_folder, run + '_envelope_timespan'+str(timespan)+'_rep'+str(m)+'.npy'), act_stack_envelope)
                    
                    # Compute the morphological active width (MAW)
                    act_stack_envelope_bool = np.where(act_stack_envelope>0,1,act_stack_envelope)
                    MAW_envelope = np.nansum(act_stack_envelope_bool,axis=2)
                    MAW_envelope = np.nansum(MAW_envelope,axis=1)
                    report_MAW_envelope = MAW_envelope/dim_x/120
                    
                    # Save envMAW array
                    header1 = 'All possible MAW envelopes (repetition) at a given timespan\nGenerated by ' + str(os.path.basename(__file__))
                    np.savetxt(os.path.join(output_folder, run + '_envelope_timespan'+str(timespan)+'_rep'+str(m)+'.txt'), report_MAW_envelope, fmt='%.4f', header=header1)
            
                
    
    # Find the length of the longest array
    max_length = max(len(arr) for arr in MAV_env_0_timespan_list)
    
    # Create an empty matrix filled with the desired value (e.g., 0)
    # Using the dtype of the first array for consistency
    padded_matrix = np.zeros((len(MAV_env_0_timespan_list), max_length), dtype=MAV_env_0_timespan_list[0].dtype)
    
    # Fill the matrix with the original arrays, padding the end of shorter arrays with the desired value
    for i, arr in enumerate(MAV_env_0_timespan_list):
        padded_matrix[i, :len(arr)] = arr
    
    padded_matrix = padded_matrix.T
    header2='Machine readable version\nGenerated by ' + str(os.path.basename(__file__))
    np.savetxt(os.path.join(output_folder, run + '_envelope_timespan0_starting_points_raw.txt'), padded_matrix, delimiter=',', header=header2, fmt='%0.4f')
            
            