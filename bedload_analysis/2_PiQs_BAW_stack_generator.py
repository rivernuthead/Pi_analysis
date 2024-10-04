#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 17:54:30 2023

@author: erri

Description:
-------------

INPUT:
-------
    BAA_map
    TYPE: NUMPY 2D ARRAY
    NOTE: BEDLOAD ACTIVE AREA MAP AS A 2D NUMPY ARRY WHERE EVERY CELL CONTAINS
        BEDLOAD INTENSITY VALUE
        
OUTPUT:
-------
    The script outputs come fron two differet sections:
        1. Section 1 creates the stack of the BAA maps for single runs
        2. Section 2 creates the stack of the BAA maps of all consecutive runs
        with the same discharge
"""
# Import packages
import numpy as np
import os
from PIL import Image
from PiQs_BAW_func_v1 import *

# SELECT RUN ------------------------------------------------------------------
# runs = ['q05r1','q05r2','q05r3','q05r4','q05r5','q05r6','q05r7','q05r8','q05r9'
#         ,'q07r1','q07r2','q07r3','q07r4','q07r5','q07r6','q07r7','q07r8','q07r9'
#         ,'q10r1','q10r2','q10r3','q10r4','q10r5','q10r6','q10r7','q10r8','q10r9'
#         ,'q15r1','q15r2','q15r3','q15r4','q15r5','q15r6','q15r7','q15r8','q15r9'
#         ,'q20r1','q20r2','q20r3','q20r4','q20r5','q20r6','q20r7','q20r8','q20r9']

# runs = ['q07r1', 'q10r1', 'q15r1', 'q20r9']
# runs = ['q07r1', 'q10r1', 'q15r1']

# runs = ['q07r1']
# runs = ['q10r1']
# runs = ['q15r1']
# runs = ['q20r9']
# runs = ['q20r1']
runs = ['q05rgm5']
set_names = ['q05_1']

#runs = ['q05r1','q05r2','q05r3','q05r4','q05r5','q05r6','q05r7','q05r8','q05r9']
# runs = ['q07r1','q07r2','q07r3','q07r4','q07r5','q07r6','q07r7','q07r8','q07r9']
# runs = ['q10r1','q10r2','q10r3','q10r4','q10r5','q10r6','q10r7','q10r8','q10r9']
# runs = ['q15r1','q15r2','q15r3','q15r4','q15r5','q15r6','q15r7','q15r8','q15r9']
# runs = ['q20r1','q20r2','q20r3','q20r4','q20r5','q20r6','q20r7','q20r8','q20r9']
# runs = ['q07rgm', 'q10rgm2', 'q15rgm2', 'q20rgm2']
# runs = ['q10rgm2', 'q15rgm2', 'q20rgm2']
# runs = ['q10rgm2', 'q20rgm2']
# runs = ['q07rgm', 'q10rgm2', 'q15rgm2']

# runs = ['q20rgm2']
# runs = ['q15rgm2']
# runs = ['q07rgm']

# SCRIPT PARAMETERS------------------------------------------------------------
downsampling_dim = 5
save_single_run_merged=1 # if 1: create a stack of all image of all runs

# FOLDERS SETUP----------------------------------------------------------------
folder_home = os.getcwd() # Setup home folder 

for run_index,run in enumerate(runs):
    print('****************')
    print(run)
    print('****************')
    
    path_output_data = os.path.join(folder_home, 'output_data')
    path_report = os.path.join(path_output_data, 'output_report',set_names[run_index], run)
    diff_path_out = os.path.join(path_output_data, '1_PiQs_BAW_maps',set_names[run_index],run)
    path_output_stack = os.path.join(path_output_data,'2_PiQs_BAW_stacks',set_names[run_index], run)

     # Check if the folders already exist and create them
    if not(os.path.exists(path_output_stack)):
        os.makedirs(path_output_stack)
    # CREATE A LIST OF THE NUMPY ARRAY OF THE BAA MAPS-------------------------
    BAA_maps_path = os.path.join(diff_path_out)
    BAA_maps_folder_names = sorted(os.listdir(BAA_maps_path))
    BAA_maps_names = []
    for name in BAA_maps_folder_names:
        if name.endswith('_ultimate_map_LR' + str(downsampling_dim) + '.npy'):
            BAA_maps_names = np.append(BAA_maps_names, name)
        else:
            pass
    
    # CREATE THE LIST OF BEDLOAD ACTIVITY MAPS --------------------------------
    BAA_map_list = []    
    for name in BAA_maps_names: # For each images in the folder...
        # LOAD THE BAA MAP (map of the bedload intensity values)
        BAA_map_path = os.path.join(diff_path_out, name) # Define the BAA map path
        BAA_map_LR = np.array(np.load(BAA_map_path)) # Load the .npy file
        BAA_map_list.append(BAA_map_LR)
    
    # CREATE THE BAA MAP STACK WITH INTENSITY VALUES
    BAA_map_stack = np.stack((BAA_map_list))
    
    # CREATE THE BOOL STACK
    BAA_map_stack_bool = (BAA_map_stack>0)*1
    np.save(os.path.join(path_output_stack, run + '_BAA_stack_LR' +
            str(downsampling_dim) + '.npy'), np.around(BAA_map_stack, decimals=0))
    
    # CREATE THE BOOL ENVELOPE STACK
    BAA_map_stack_env = np.nansum(BAA_map_stack_bool, axis =0)
    BAA_map_stack_env_bool = (BAA_map_stack_env>0)*1
    np.save(os.path.join(path_output_stack, run + '_BAA_stack_envelope_LR' +
            str(downsampling_dim) + '.npy'), np.around(BAA_map_stack_env, decimals=2))
    
    # CREATE THE ACTIVE TIMES MAP (Number of time a certain pixel has been active)
    BAA_active_times_map = np.nansum(BAA_map_stack_bool, axis=0)
    np.save(os.path.join(path_output_stack, run + '_envBAA_act_period_LR' +
            str(downsampling_dim) + '.npy'), np.around(BAA_active_times_map, decimals=2))
    
    # CREATE THE BAA_MEAN INTENSITY  (mean of the intensity without zeros)
    BAA_mean_intensity_map = np.nanmean(np.where(BAA_map_stack==0, np.nan, BAA_map_stack) , axis=0)
    np.save(os.path.join(path_output_stack, run + '_envBAA_act_mean_intensity_LR' +
            str(downsampling_dim) + '.npy'), np.around(BAA_mean_intensity_map, decimals=2))

    envBAA_act_cumulative_img = Image.fromarray(np.array(BAA_mean_intensity_map*1).astype(np.uint16))
    envBAA_act_cumulative_img.save(os.path.join(diff_path_out, run + '_envBAA_act_cumulative.tiff'))
    # envBAA_act_cumulative_img = Image.fromarray(np.array(envBAA_act_cumulative*1))
    # envBAA_act_cumulative_img.save(os.path.join(diff_path_out, run, run + '_envBAA_act_cumulative.tiff'))

#%%
#-----------------------------------------------------------------------------#
# MERGE SINGLE RUNS IN A UNIQUE STACK
#-----------------------------------------------------------------------------#
"""

This script merge all the single run PiQs maps (q07r1, q07r2, ...) and create
the '_envBAA_single_runs_merged_LR5.npy' file that is the stack of all the
maps. This was necessary due to the slight differences between BAW and env BAW
from regime runs and the single runs (in particular for the q07_1 set)
"""

# FUNCTIONS--------------------------------------------------------------------
if save_single_run_merged==1:
    def collect_and_stack_files(base_dir, subfolder_filter, file_filter, exclude_filter):
        stack = []
    
        for root, dirs, files in os.walk(base_dir):
            # Check if the subfolder contains "07" in the name and does not contain "rgm"
            if subfolder_filter in os.path.basename(root) and exclude_filter not in os.path.basename(root):
                for file in files:
                    # Check if the file contains "ultimate_map_LR5.npy" in the name
                    if file_filter in file:
                        full_path = os.path.join(root, file)
                        # print(full_path)
                        stack.append(full_path)
    
        # Sort the stack
        stack.sort()
    
        npy_stack = []
        for file_path in stack:
            
            data = np.load(file_path)
            npy_stack.append(data)
    
        # Convert list of arrays to a single stacked array
        npy_stack = np.stack(npy_stack, axis=0)
    
        return npy_stack
    
    # =============================================================================
    # LOOP OVER RUNS
    # =============================================================================
    # SEARCH PARAMETERS -----------------------------------------------------------
    base_dir =  'output_data/1_PiQs_BAW_maps'
    file_filter = "ultimate_map_LR5.npy"
    exclude_filter = "rgm"

    for subfolder_filter, set_name in zip(['q05', 'q07', 'q10', 'q15', 'q20'], ['q05_1','q07_1', 'q10_2', 'q15_2', 'q20_2']):
        # Collect and stack the files
        npy_stack = collect_and_stack_files(
            base_dir, subfolder_filter, file_filter, exclude_filter)
        
        np.save(os.path.join(path_output_data,'2_PiQs_BAW_stacks',run[0:3], set_name +
                '_single_run_merged_BAA_stack_LR' + str(downsampling_dim) + '.npy'), npy_stack)

        # Print the shape of the stacked array
        print(f"Stacked single runs merged array shape: {npy_stack.shape}")
        