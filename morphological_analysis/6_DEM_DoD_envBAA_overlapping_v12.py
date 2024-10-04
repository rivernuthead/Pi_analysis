#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 15:06:13 2023

@author: erri

The aim of this script is to overlap camera photo over topographic surveys
information.

Outputs:
 - 3D stack: [arrays, dy, dx] ; arrays =  [DEM_rsz_rsc, DoD_rsz_rsc, envBAA_act_period_rsh, envBAA_intensity_rsh_filt, envBAA_act_period_4th]
"""


# =============================================================================
# IMPORT PACKAGES
# =============================================================================
import os
import cv2
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import tifffile as tiff 

#==============================================================================
# OPTION
#==============================================================================
DoS_map_overlap = 0 # if 1 also overlap the envelop map of the mean Difference of saturation.

# =============================================================================
# FUNCTIONS
# =============================================================================
def array_to_bool(arr, mode):
    if mode == 'semi-bool':
        # Create a copy of the array to avoid modifying the original
        result = np.copy(arr)
        # Apply the transformation for mode 1
        result[~np.isnan(result) & (result < 0)] = -1
        result[~np.isnan(result) & (result > 0)] = 1
    elif mode == 'bool':
        # Create a copy of the array to avoid modifying the original
        result = np.copy(arr)
        # Apply the transformation for mode 2
        result[~np.isnan(result) & (result != 0)] = 1
        result[~np.isnan(result) & (result == 0)] = 0
    else:
        raise ValueError("Invalid mode. Choose either mode 'semi_bool' or mode 'bool'.")
    return result

def img_scaling_to_DEM(image, scale, dx, dy, rot_angle):

    # Create the transformation matrix
    M1 = np.float32([[scale, 0, dx], [0, scale, dy]])
    
    # Apply the transformation to img1 and store the result in img2
    rows, cols = image.shape
    image_rsh = cv2.warpAffine(image, M1, (cols, rows))
    
    # Rotate the image
    M2 = cv2.getRotationMatrix2D((image_rsh.shape[1]/2, image_rsh.shape[0]/2), rot_angle, 1)
    image_rsh = cv2.warpAffine(image_rsh, M2, (image_rsh.shape[1], image_rsh.shape[0]))
    
    # Trim zeros rows and columns due to shifting
    x_lim, y_lim = dx+int(cols*scale), dy+int(rows*scale)
    image_rsh = image_rsh[:y_lim, :x_lim]

    return image_rsh

def assign_closest_values(matrix, target_values=[0, 0.25, 0.5, 0.75, 1]):
    """
    This function replaces values in the input matrix with the closest values
    from the target_values list.
    
    Parameters:
    matrix (np.ndarray): A matrix with values ranging from 0 to 1.
    target_values (list): A list of target values to assign based on proximity. 
                          Defaults to [0, 0.25, 0.5, 0.75, 1].
    
    Returns:
    np.ndarray: A matrix with values replaced by the closest target values.
    """
    target_values = np.array(target_values)
    
    # Find the closest target value for each element in the matrix
    closest_values = target_values[np.abs(matrix[:, :, None] - target_values).argmin(axis=2)]
    
    return closest_values


run_mode = [
    'full_BAA_stack',
    # 'partial_envBAA_stack'
    ]
'''
run_mode.
        - full_BAA_stack: full BAA map stack
        - partial_envBAA_stack: stack computed as the stack of partial BAA envelopes (see /PiQs_analysis/4_PiQs_BAW_envelopes_generator.py for details)
'''

home_dir = os.getcwd() # Home directory
report_dir = os.path.join(home_dir, 'output_data')
output_folder = os.path.join(report_dir,'DEM_DoD_envBAA_overlapping')
run_dir = os.path.join(home_dir, 'input_data','surveys')
PiQs_dir = os.path.join(home_dir, '..','bedload_analysis')
stack_dir = os.path.join(home_dir, 'output_data', 'DoDs', 'DoDs_stack')

if not(os.path.exists(output_folder)):
    os.mkdir(output_folder)


#run_names = ['q05r1', 'q05r2', 'q05r3', 'q05r4', 'q05r5', 'q05r6', 'q05r7', 'q05r8', 'q05r9']

run_names = ['q07r1', 'q07r2', 'q07r3', 'q07r4', 'q07r5', 'q07r6', 'q07r7', 'q07r8', 'q07r9']
# run_names = ['q10r1', 'q10r2', 'q10r3', 'q10r4', 'q10r5', 'q10r6', 'q10r7', 'q10r8', 'q10r9']

# # run_names = ['q15r1', 'q15r2', 'q15r3', 'q15r4', 'q15r5', 'q15r6', 'q15r7', 'q15r8', 'q15r9']

# run_names = ['q15r1', 'q15r3', 'q15r4', 'q15r5', 'q15r6', 'q15r8', 'q15r9']

# run_names = ['q20r1', 'q20r2', 'q20r3', 'q20r4', 'q20r5', 'q20r6', 'q20r7', 'q20r8', 'q20r9']

# run_names = ['q05r1']

# run_names = ['q07r1']

# run_names = ['q10r1']

# run_names = ['q15r1']

# run_names = ['q20r1']

mask = tiff.imread(os.path.join(home_dir, 'input_data','masks','mask_DEM_DoD_envBAA.tiff'))

# ENVELOPES TIMESCALE FOR THE envBAA COMPUTATION-------------------------------
#env_tscale_array = [2,3,4,9,20]       #1/4 of the frame number
# env_tscale_array = [3,4,6,9,21]        # 1/8 Txnr
# env_tscale_array = [5,7,12,18,42]      # 1/4 Txnr
env_tscale_array = [10,15,24,37,84]    # 1/2 Txnr
# env_tscale_array = [15,22,36,55, 126]    # 3/4 Txnr
# env_tscale_array = [20,30,48,74,168]    # 1 Txnr

# envBAA LIST, ONE FOR EACH RUN------------------------------------------------
envBAA_intensity_list = [] # Collect all the intensity maps
envBAA_act_period_list = [] # Collect all the activity period
envBAA_act_period_4th_list = [] # Collect all the activity periods as a fraction of the run duratio (1/4)

for run_name in run_names:
    
    print(run_name, ' is running...')
    run_index = int(run_name[-1])
    
    # DEFINE RUN PARAMETERS
    if '05' in run_name:
        set_name = 'q05_1'
        env_tscale = env_tscale_array[4]
        # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
        thrs_actDS = 6.5 
        thrs_actUS = 6.5
        if 'rgm' in run_name:
            thrs_actDS = 6.5
            thrs_actUS = 6.5
    if '07' in run_name:
        set_name = 'q07_1'
        env_tscale = env_tscale_array[3]
        thrs_actDS = 6.5
        thrs_actUS = 6.4
        if 'rgm' in run_name:
            thrs_actDS = 6.5
            thrs_actUS = 6.4
    if '10' in run_name:
        set_name = 'q10_2'
        env_tscale = env_tscale_array[2]
        thrs_actDS = 5.3
        thrs_actUS = 5.8
        if 'rgm' in run_name:
            thrs_actDS = 6.4
            thrs_actUS = 6.0
    if '15' in run_name:
        set_name = 'q15_2'
        env_tscale = env_tscale_array[1]
        thrs_actDS = 7.0
        thrs_actUS = 7.0
        if 'rgm' in run_name:
            thrs_actDS = 6.0
            thrs_actUS = 6.0
    if '20' in run_name:
        set_name = 'q20_2'
        env_tscale = env_tscale_array[0]
        thrs_actDS = 7.0
        thrs_actUS = 5.7
        if 'rgm' in run_name:
            thrs_actDS = 7.0
            thrs_actUS = 5.7
    
    # =========================================================================
    #   IMPORT DATA 
    # =========================================================================
    
    # IMPORT BAA MAPS ---------------------------------------------------------
    if 'full_BAA_stack' in run_mode:
        ''' Values range from 0 (never active) to the number of frames for each run (always active) - file from 2_PiQs_BAW_stack_generator_v7.py''' 
        print('Script mode: Full stack')
        #                                        '/home/erri/Documents/PhD/Research/5_research_repos/PiQs_analysis/activity_stack/q05r1_envBAA_act_period_LR5.npy'
        #                                        '/home/erri/Documents/PhD/Research/5_research_repos/PiQs_analysis/activity_stack/q05r1_envBAA_act_period_LR5.npy'
        envBAA_act_period = np.load(os.path.join(PiQs_dir,'output_data','2_PiQs_BAW_stacks',run_name[0:3], run_name, run_name + '_envBAA_act_period_LR5.npy'))
        envBAA_act_period = np.array(envBAA_act_period, dtype=np.uint16)
        
        test = np.load(os.path.join(PiQs_dir,'output_data','2_PiQs_BAW_stacks',run_name[0:3], run_name, run_name + '_envBAA_act_period_LR5.npy'))
        
        
    elif 'partial_envBAA_stack' in run_mode:
        ''' Given the temporal resolution at which the envelope has been taken, the values range from 0 to 4 in T envelope is 1/8 of the Txnr'''
        print('Script mode: Partial stack')
        partial_envBAA_path = os.path.join('/home/erri/Documents/PhD/Research/5_research_repos/PiQs_analysis/activity_stack/activity_stack_cleaned/envelopes_cleaned/space_discr_5')
        envBAA_act_period = np.load(os.path.join(partial_envBAA_path, run_name, run_name + '_envTscale' + str(env_tscale),  run_name + '_envT' + str(env_tscale) + '_partial_envBAA_sum_stack.npy')) # Output from /PiQs_analysis/4_PiQs_BAW_envelopes_generator.py
        envBAA_act_period = np.array(envBAA_act_period, dtype=np.uint16)
    
    
    
    # IMPORT envBAA INTENSITY MAP ---------------------------------------------
    envBAA_intensity = np.load(os.path.join(PiQs_dir,'output_data','2_PiQs_BAW_stacks',run_name[0:3], run_name, run_name + '_envBAA_act_mean_intensity_LR5.npy'))
    envBAA_intensity = np.array(envBAA_intensity, dtype=np.uint16)
    
    # IMPORT DEM MAPS ---------------------------------------------------------
    # DEM = np.loadtxt(os.path.join(run_dir, set_name, 'matrix_bed_norm_' + set_name +'s'+ str(run_names.index(run_name))+'.txt'), skiprows=8)
    DEM_stack = np.load(os.path.join(home_dir, 'output_data', 'DEMs','DEMs_stack', set_name + '_DEM_stack.npy'))
    DEM = DEM_stack[run_index-1,:,:]
    
    # IMPORT DoD MAP ----------------------------------------------------------
    # DoD = np.loadtxt(os.path.join(report_dir, 'DoDs', 'DoDs_'+ set_name,'DoD_' + str(run_names.index(run_name)+1) + '-' + str(run_names.index(run_name)) + '_filt_ult.txt'))
    DoD_stack = np.load(os.path.join(home_dir, 'output_data','DoDs', 'DoDs_stack','DoD_stack_' + set_name + '.npy'))
    DoD_stack = DoD_stack[:,:,:,0] # Extract timespan 1 DoDs
    DoD_stack_bool = array_to_bool(DoD_stack, 'semi-bool')
    DoD = DoD_stack[run_index-1,:,:]

    # OPTIONAL: IMPORT ENV DOS MAP
    if DoS_map_overlap == 1 :
        envDoS_intensity = np.load(os.path.join(PiQs_dir,'output_data','3_Image_filtering_BAW_map','DoS_maps',run_name,'envelop', run_name + '_DoS_mean_map_LR5.npy'))
        envDoS_intensity = np.array(envDoS_intensity, dtype=np.uint16)
    # -------------------------------------------------------------------------
    
    # TEMPORARY RESHAPING FOR THE q05_1 RUN -----------------------------------
    # TODO REMOVE THIS PART AFTER PiQs q05_1 RESHAPING
    # Desired dimensions
    desired_shape = (140, 1252)
    
    # Check if the current matrix shape is (132, 1252)
    if envBAA_act_period.shape == (132, 1252):
        # Padding calculation
        rows_to_pad = desired_shape[0] - envBAA_act_period.shape[0]  # 140 - 132 = 8
        pad_top = rows_to_pad // 2  # 4 rows on top
        pad_bottom = rows_to_pad - pad_top  # 4 rows on bottom
        
        # Create the padding with np.nan
        padding_top = np.full((pad_top, envBAA_act_period.shape[1]), 0)
        padding_bottom = np.full((pad_bottom, envBAA_act_period.shape[1]), 0)
        envBAA_act_period = np.vstack([padding_top, envBAA_act_period, padding_bottom])
        
        padding_top = np.full((pad_top, envBAA_intensity.shape[1]), 0)
        padding_bottom = np.full((pad_bottom, envBAA_intensity.shape[1]), 0)
        envBAA_intensity = np.vstack([padding_top, envBAA_intensity, padding_bottom])

        if DoS_map_overlap == 1 :
            padding_top = np.full((pad_top, envDoS_intensity.shape[1]), 0)
            padding_bottom = np.full((pad_bottom, envDoS_intensity.shape[1]), 0)
            envDoS_intensity = np.vstack([padding_top, envDoS_intensity, padding_bottom])
            
    # -------------------------------------------------------------------------
    
    
    # DEM RESIZE AND RESCALE---------------------------------------------------
    # Number of rows and columns to pad
    pad_rows = 0
    pad_cols = 16
    DEM = np.where(DEM==-999, np.nan, DEM) # Set NaN as np.nan
    # if '05' in run_name:
    #     DEM_rsz = np.copy(DEM)[14:146,157:]  # Resize DEM to fit the envBAA map domain
    # else:
    DEM_rsz = np.copy(DEM)[4:144,157:]  # Resize DEM to fit the envBAA map domain
        
    DEM_rsz_rsc = np.repeat(DEM_rsz, 10, axis=1) # Rescale the resized DEM to fit the envBAA map scale
    DEM_rsz_rsc = np.pad(DEM_rsz_rsc, ((pad_rows, 0), (pad_cols, 0)), constant_values=np.nan)
    
    # DoD RESIZE AND RESCALE---------------------------------------------------
    # Number of rows and columns to pad
    # if '05' in run_name:
    #     DoD_rsz = np.copy(DoD)[14:146,157:] # Resize DEM to fit the envBAA map domain
    # else:
    DoD_rsz = np.copy(DoD)[4:144,157:] # Resize DEM to fit the envBAA map domain
    DoD_rsz_rsc = np.repeat(DoD_rsz, 10, axis=1) # Rescale the resized DEM to fit the envBAA map scale
    DoD_rsz_rsc = np.pad(DoD_rsz_rsc, ((pad_rows, 0), (pad_cols, 0)), constant_values=np.nan)
    
    # RESIZE envBAA------------------------------------------------------------
    envBAA_act_period = envBAA_act_period[:,:DoD_rsz_rsc.shape[1]]
    envBAA_intensity = envBAA_intensity[:,:DoD_rsz_rsc.shape[1]]
    if DoS_map_overlap == 1 :
        envDoS_intensity = envDoS_intensity[:,:DoD_rsz_rsc.shape[1]]
        
    ''' Sometimes runs have different numbers of frames. This set of tscale is
    designed to investigate close to the 1/8 of Txnr. Whith this setup it may
    occurs that the envelope calculation has the number of frame that allows
    another step so isted of finish with 4 it finish with 5. Insted of delating
    raw photos, and considering the two classes (active 1/8 of the run and
    active the full run) it doesn't matter if the full runn is 4 or 5 then we
    decided to put averything exceedes 4 at 4.'''
    if env_tscale_array==[2,3,4,9,20] and 'partial_envBAA_stack' in run_mode:
        envBAA_act_period = np.where(envBAA_act_period>4, 4, envBAA_act_period)
    
    if set_name == 'q05_1':
        # Define the transformation parameters
        scale = 1.0 # Enlargement scale
        dx = 2 # Shift in x direction
        dy = 4 # Shift in y direction
        rot_angle = -0.42
    if set_name == 'q07_1':
        # Define the transformation parameters
        scale = 1.01 # Enlargement scale
        dx = 0 # Shift in x direction
        dy = 5 # Shift in y direction
        rot_angle = -0.421
    if set_name == 'q10_2':
        # Define the transformation parameters
        scale = 1.0 # Enlargement scale
        dx = 2 # Shift in x direction
        dy = 4 # Shift in y direction
        rot_angle = -0.42
    if set_name == 'q15_2':
        # Define the transformation parameters
        scale = 1.0 # Enlargement scale
        dx = 2 # Shift in x direction
        dy = 4 # Shift in y direction
        rot_angle = -0.42
    if set_name == 'q20_2':
        # Define the transformation parameters
        scale = 1.0 # Enlargement scale
        dx = 2 # Shift in x direction
        dy = 4 # Shift in y direction
        rot_angle = -0.42
    
    
    
    
    # ADJUST AND RESCALE envBAA MAP
    envBAA_intensity = np.array(envBAA_intensity, dtype=np.uint16)
    envBAA_intensity_rsh = img_scaling_to_DEM(envBAA_intensity, scale, dx, dy, rot_angle)
    envBAA_act_period = np.array(envBAA_act_period, dtype=np.uint16)
    envBAA_act_period_rsh = img_scaling_to_DEM(envBAA_act_period, scale, dx, dy, rot_angle)
    if DoS_map_overlap == 1 :
        envDoS_intensity = np.array(envDoS_intensity, dtype=np.uint16)
        envDoS_intensity = img_scaling_to_DEM(envDoS_intensity, scale, dx, dy, rot_angle)
    
    # =============================================================================
    # APPLY ACTIVITY THRESHOLD
    # =============================================================================
    envBAA_intensity_rsh_filt = np.copy(envBAA_intensity_rsh)
    
    for y in range(1,envBAA_intensity_rsh.shape[0]-1):
        for x in range(1,envBAA_intensity_rsh.shape[1]-1):
            ker = [[envBAA_intensity_rsh[y-1,x-1], envBAA_intensity_rsh[y-1,x], envBAA_intensity_rsh[y-1,x+1]],
                   [envBAA_intensity_rsh[y,x-1], envBAA_intensity_rsh[y,x], envBAA_intensity_rsh[y,x+1]],
                   [envBAA_intensity_rsh[y+1,x-1], envBAA_intensity_rsh[y+1,x], envBAA_intensity_rsh[y+1,x+1]]]
            
            ker=ker[ker!=0]
            
            # Define activity threshold
            if x<envBAA_intensity_rsh.shape[1]//2:
                thrs_act = thrs_actUS
            else:
                thrs_act = thrs_actDS
              
            if envBAA_intensity_rsh_filt[y,x]>=thrs_act:
                envBAA_intensity_rsh_filt[y,x] = envBAA_intensity_rsh[y,x]
            elif np.logical_and(envBAA_intensity_rsh_filt[y,x]<thrs_act, np.nanmean(ker)>0.8*thrs_act):
                envBAA_intensity_rsh_filt[y,x] = thrs_act
            else:
                envBAA_intensity_rsh_filt[y,x] = 0
            
    envBAA_intensity_rsh_filt = envBAA_intensity_rsh_filt*(envBAA_intensity_rsh_filt>np.nanmin((thrs_actUS,thrs_actDS))) # remove sparse points
    
    envBAA_rsh_mask = envBAA_intensity_rsh_filt>0
    envBAA_act_period_rsh = envBAA_act_period_rsh*envBAA_rsh_mask
    
    
    # STACK DEM, DoD AND envBAA MAP -------------------------------------------
    # DEM_DoD_envBAA_stack[layer,dim_y, dim_x] where layer0=DEM, layer1=DoD, layer2=envBAA, layer3=envBAA intensity
    
    # APPEND envBAA INTENSITY MAP
    envBAA_intensity_list.append(envBAA_intensity_rsh_filt)
    
    # APPEND envBAA ACTIVE PERIOD MAP
    envBAA_act_period_list.append(envBAA_act_period_rsh)
    
    # APPEND envBAA ACTIVE PERIODS MAP AS A 1/4 OF THE TUN DURATIUON
    envBAA_act_period_4th = (assign_closest_values(envBAA_act_period_rsh/np.max(envBAA_act_period_rsh)))*4
    envBAA_act_period_4th_list.append(envBAA_act_period_4th)
    
    # PREPARE FILE FOR SAVING
    envBAA_plot = np.where(envBAA_act_period_rsh>0,100,0)
    DEM_plot = np.where(DEM_rsz_rsc ==0, np.nan, DEM_rsz_rsc)
    DoD_plot = np.where(DoD_rsz_rsc ==0, np.nan, DoD_rsz_rsc)
    
    # SAVE IMAGES
    output_path = os.path.join(report_dir, 'surveys_images_overlapping_test', set_name)
    if not(os.path.exists(output_path)):
        os.makedirs(output_path)
    DoD_plot_img = Image.fromarray(np.array(DoD_plot*1).astype(np.uint16))
    DoD_plot_img.save(os.path.join(output_path, run_name + '_DoD_resized_image.tiff'))
    envBAA_plot_img = Image.fromarray(np.array(envBAA_plot*1).astype(np.uint16))
    envBAA_plot_img.save(os.path.join(output_path, run_name + '_envBAA_resized_image.tiff'))
    
    # =============================================================================
    # CREATE THE PLOT
    # =============================================================================
    fig, axs = plt.subplots(figsize=(15, 6))
    # img1 = plt.imshow(DEM_plot, cmap='binary', origin='upper', alpha=1.0, vmin=-20, vmax=20, interpolation_stage='rgba')
    img2 = plt.imshow(DoD_plot, cmap='binary', origin='upper', alpha=1.0, vmin=-20, vmax=20, interpolation_stage='rgba')
    img2 = plt.imshow(envBAA_plot, alpha=0.5, cmap='cool', origin='upper', vmin=-10, vmax=+10, interpolation_stage='rgba')
    plt.title(run_name)
    plt.axis('off')
    plt.text(0.5, -0.4, f"Generated by: {os.path.basename(__file__)}", transform=plt.gca().transAxes, ha='center', va='center', fontsize=12)
    if 'full_BAA_stack' in run_mode:
        plt.savefig(os.path.join(output_path, run_name + '_DEM_DoD_envBAA_overlapping.tiff'), dpi=600)
    elif 'partial_envBAA_stack' in run_mode:
        plt.savefig(os.path.join(output_path, run_name + '_DEM_DoD_partial_envBAA_overlapping.tiff'), dpi=600)
    plt.show()

    # =========================================================================
    # STACK DEM, DOD, envBAA ACTIVE PERIODS, AND envBAA ACTIVE INTENSITY
    # =========================================================================
    if DoS_map_overlap == 1 :
        DEM_DoD_envBAA_stack = np.stack([DEM_rsz_rsc, DoD_rsz_rsc, envBAA_act_period_rsh, envBAA_intensity_rsh_filt, envBAA_act_period_4th, envDoS_intensity], axis=0)
    else:
        DEM_DoD_envBAA_stack = np.stack([DEM_rsz_rsc, DoD_rsz_rsc, envBAA_act_period_rsh, envBAA_intensity_rsh_filt, envBAA_act_period_4th], axis=0)
    
    # CREATE AND APPLY MASK ---------------------------------------------------
    mask = np.pad(mask, ((pad_rows, 0), (pad_cols, 0)), constant_values=np.nan)
    mask = mask[:DEM_DoD_envBAA_stack.shape[1], :DEM_DoD_envBAA_stack.shape[2]]
    
    for i in range(DEM_DoD_envBAA_stack.shape[0]):
        DEM_DoD_envBAA_stack[i, :, :] = DEM_DoD_envBAA_stack[i, :, :]#*mask
    
    # SAVE THE STACK ----------------------------------------------------------
    if not(os.path.exists(os.path.join(output_folder,set_name))):
        os.mkdir(os.path.join(output_folder,set_name))
        
    if 'full_BAA_stack' in run_mode:
        if DoS_map_overlap == 1 :
            path = os.path.join(output_folder, set_name, set_name + '_' + run_name + '_DEM_DoD_envBAA_DoS_stack.npy')
        else:
            path = os.path.join(output_folder, set_name, set_name + '_' + run_name + '_DEM_DoD_envBAA_stack.npy')
    elif 'partial_envBAA_stack' in run_mode:
        path = os.path.join(output_folder, set_name, set_name + '_' + run_name + '_DEM_DoD_partial_envBAA_stack.npy')
    np.save(path, DEM_DoD_envBAA_stack)


# STACK envBAA MAPS. ONE FOR EACH RUN -----------------------------------------

# This is the mean intensity for each run
envBAA_set_name_stack = np.stack(envBAA_intensity_list)
np.save(os.path.join(output_folder, set_name, set_name + '_envBAA_intensity_maps_single_runs.npy'), envBAA_set_name_stack)

# Number of active times
envBAA_act_period_stack = np.stack(envBAA_act_period_list)
np.save(os.path.join(output_folder, set_name, set_name + '_envBAA_act_period_maps_single_runs.npy'), envBAA_act_period_stack)

# Number of active times as 1/4 of the total run duration
envBAA_act_period_stack_4th = np.stack(envBAA_act_period_4th_list)
np.save(os.path.join(output_folder, set_name, set_name + '_envBAA_act_4th_period_maps_single_runs.npy'), envBAA_act_period_stack_4th)

#%%============================================================================
# COMPENSATION AND UNDER THRS envBAA OVERLAPPING
# =============================================================================
'''
The aim of this section is:
    Extract bedload activity information where compensation and filtering
    process affects the MAA.
    
    DoD information are used at timespan 1 the be coherent with the way I
    computed the compensation and the filtered effects maps (I need two DoD to
    comupte one comp or thrs map so I use a DoD that is coherent with the timescale)
    
    q15r2 and q15r7 fails but are still present. Maybe we shouls remove q15 and
    add q05.
'''
def calculate_envBAA_map_sums(envBAA_maps):
    # Ensure envBAA_maps is a numpy array
    envBAA_maps = np.array(envBAA_maps)
    
    # Create the result array with the same shape as the input minus the first time dimension
    result_shape = (envBAA_maps.shape[0] - 1, envBAA_maps.shape[1], envBAA_maps.shape[2])
    envBAA_map_sums = np.zeros(result_shape)

    # Iterate over the time dimension
    for t in range(result_shape[0]):
        # Extract consecutive time slices
        slice1 = envBAA_maps[t]
        slice2 = envBAA_maps[t+1]

        # Calculate the mean where both are zeros or both are greater than zero
        mean_mask = (slice1 == 0) & (slice2 == 0) | (slice1 > 0) & (slice2 > 0)
        envBAA_map_sums[t][mean_mask] = (slice1[mean_mask] + slice2[mean_mask]) / 2

        # Where only one entry is non-zero, take the non-zero entry
        non_zero_mask1 = (slice1 != 0) & (slice2 == 0)
        non_zero_mask2 = (slice1 == 0) & (slice2 != 0)
        envBAA_map_sums[t][non_zero_mask1] = slice1[non_zero_mask1]
        envBAA_map_sums[t][non_zero_mask2] = slice2[non_zero_mask2]

    return envBAA_map_sums


import os
import numpy as np

home_dir = os.getcwd()
report_dir = report_dir#os.path.join(home_dir, 'outputs')
stack_dir = stack_dir#os.path.join(home_dir, 'outputs', 'DoDs', 'DoDs_stack')






# =============================================================================
# 
# =============================================================================


# for set_name in ['q07_1', 'q10_2', 'q15_2', 'q20_2']:
# for set_name in ['q07_1']:
for set_name in ['q05_1', 'q07_1', 'q10_2', 'q20_2']:
    print(set_name + ' is running...')
    
    # IMPORT DATA -------------------------------------------------------------
    envBAA_intensity_maps = np.load(os.path.join(output_folder, set_name, set_name + '_envBAA_intensity_maps_single_runs.npy')) # Bedload mean intensity
    envBAA_active_period_maps = np.load(os.path.join(output_folder, set_name, set_name + '_envBAA_act_period_maps_single_runs.npy')) # Number of active 1/4 of the run duration. (4 means 4/4=1 always active)
    DoD_maps    = np.load(os.path.join(stack_dir,"DoD_stack_"+set_name+".npy"))
    comp_maps   = np.load(os.path.join(stack_dir, set_name + '_compensation_map_stack.npy'))
    thrs_maps   = np.load(os.path.join(stack_dir, set_name + '_threshold_pure_map_stack.npy'))
    
    # RESHAPING BEDLOAD ACTIVITY MAPS
    
    
    # REMOVE MAPS WHERE q15_2 FAILS -------------------------------------------
    if set_name == 'q15_2':
        DoD_maps = np.delete(DoD_maps, [2, 7], axis=0)
        comp_maps = np.delete(comp_maps, [2, 7], axis=0)
        thrs_maps = np.delete(thrs_maps, [2, 7], axis=0)
    
    # CREATE MAPS -------------------------------------------------------------
    
    # Number of active times (1 Exner time discretization for comp and thrs comparison)
    envBAA_map_activity_times = np.where(envBAA_intensity_maps>0,1,0)
    
    
    envBAA_map_activity_times_1Txnr = envBAA_map_activity_times[:-1,:,:] + envBAA_map_activity_times[1:,:,:]
    
    envBAA_active_period_maps_1Txnr = envBAA_active_period_maps[:-1,:,:] + envBAA_active_period_maps[1:,:,:]
    
    envBAA_map_activity_times_1Txnr = np.where(envBAA_map_activity_times_1Txnr>0,1,0)
    envBAA_map_activity_times_1Txnr = np.nansum(envBAA_map_activity_times_1Txnr, axis=0)
    
    # Mean activity time between two consecutive bedload map
    envBAA_intensity_map_sums = calculate_envBAA_map_sums(envBAA_intensity_maps)
    
    # Morphological active map
    DoD_maps  = DoD_maps[:-1,:,:,1] # DoD information are take at timespan 1. See the script notes for details
    
    # Map of compensation area
    comp_maps = comp_maps[:,:,:,0]
    
    # Map of filtering process area
    thrs_maps = thrs_maps[:,:,:,0]
    
    # RESIZE AND RESCALE TOPOGRAPHIC INFORMATION ------------------------------
    DoD_maps_rsz = np.copy(DoD_maps)[:,4:144,157:] # Resize DEM to fit the envBAA map domain 157
    DoD_maps_rsz_rsc = np.repeat(DoD_maps_rsz, 10, axis=2) # Rescale the resized DEM to fit the envBAA map scale
    
    comp_maps_rsz = np.copy(comp_maps)[:,4:144,157:] # Resize DEM to fit the envBAA map domain 157
    comp_maps_rsz_rsc = np.repeat(comp_maps_rsz, 10, axis=2) # Rescale the resized DEM to fit the envBAA map scale
    comp_maps_rsz_rsc = np.where(comp_maps_rsz_rsc==0, np.nan, comp_maps_rsz_rsc) # Set np.nan as transparent
    
    thrs_maps_rsz = np.copy(thrs_maps)[:,4:144,157:] # Resize DEM to fit the envBAA map domain 157
    thrs_maps_rsz_rsc = np.repeat(thrs_maps_rsz, 10, axis=2) # Rescale the resized DEM to fit the envBAA map scale
    thrs_maps_rsz_rsc = np.where(thrs_maps_rsz_rsc==0, np.nan, thrs_maps_rsz_rsc) # Set np.nan as transparent
    
    # APPLY MASK WHERE COMPENSATION OR FILTERING EFFECTS HAPPEN ---------------
    # Bedload intensity
    envBAA_intensity_map_sums = envBAA_intensity_map_sums[:,:,:DoD_maps_rsz_rsc.shape[2]]
    envBAA_comp = envBAA_intensity_map_sums*comp_maps_rsz_rsc
    envBAA_thrs = envBAA_intensity_map_sums*thrs_maps_rsz_rsc
    
    # Number of bedload active periods
    envBAA_active_period_maps_1Txnr = envBAA_active_period_maps_1Txnr[:,:,:DoD_maps_rsz_rsc.shape[2]]
    envBAA_periods_comp = envBAA_active_period_maps_1Txnr*comp_maps_rsz_rsc
    envBAA_periods_thrs = envBAA_active_period_maps_1Txnr*thrs_maps_rsz_rsc
    
    # SAVE STACKS -------------------------------------------------------------
    np.save(os.path.join(report_dir, set_name + '_DoD_1Txnr_stack.npy'), DoD_maps_rsz_rsc) # DoD maps
    np.save(os.path.join(report_dir, set_name + '_envBAA_intensity_1Txnr_stack.npy'), envBAA_intensity_map_sums) # Bedload intensity
    np.save(os.path.join(report_dir, set_name + '_envBAA_active_periods_1Txnr_stack.npy'), envBAA_active_period_maps_1Txnr) # Bedload active periods as the number of 1/4 of hte single run duration. 4 means 0.5 Txnr, 8 means 1 Txnr
    np.save(os.path.join(report_dir, set_name + '_envBAA_intensity_comp_stack.npy'), envBAA_comp) # Bedload intensity where comp occurs
    np.save(os.path.join(report_dir, set_name + '_envBAA_intensity_thrs_stack.npy'), envBAA_thrs) # Bedload intensity where filt process occurs
    np.save(os.path.join(report_dir, set_name + '_envBAA_periods_comp_stack.npy'), envBAA_periods_comp) # Bedload activity periods where comp occurs
    np.save(os.path.join(report_dir, set_name + '_envBAA_periods_thrs_stack.npy'), envBAA_periods_thrs) # Bedload activity periods filtering process occurs

    



#%%============================================================================
# MAPS PLOTS
# =============================================================================
'''
The aim of this section is produce plots of overlapped information about bdlaod
activity, morphological changes, compensation and filtering process effects.


q15_3 envBAA is not available
'''

from matplotlib.colors import ListedColormap

set_names = ['q05_1', 'q07_1', 'q10_2', 'q20_2']
# set_names = ['q07_1']

for set_name in set_names:
    

    # IMPORT STACK ------------------------------------------------------------
    # DoD maps stack
    DoD_maps_stack    = np.load(os.path.join(home_dir, 'output_data','DoDs', 'DoDs_stack',"DoD_stack_"+set_name+".npy"))
    # Bedlaod mean intensity maps stack
    envBAA_stack = np.load(os.path.join(report_dir, set_name + '_envBAA_intensity_1Txnr_stack.npy')) # envBAA maps envelope at 1 Exner time
    # Morphological compensation maps stack
    compensation_stack   = np.load(os.path.join(stack_dir, set_name + '_compensation_map_stack.npy'))
    # Filtering process maps stack
    threshold_pure_stack   = np.load(os.path.join(stack_dir, set_name + '_threshold_pure_map_stack.npy'))
    
    # RESIZE AND RESCALE TOPOGRAPHIC INFORMATION ------------------------------
    # Resize DoD to fit the envBAA map domain 157
    DoD_maps = np.copy(DoD_maps_stack)[:,4:144,157:]
    
    # Resize DoD to fit the envBAA map domain 157
    compensation_stack = np.copy(compensation_stack)[:,4:144,157:] 
    
    # Resize DoD to fit the envBAA map domain 157
    threshold_pure_stack = np.copy(threshold_pure_stack)[:,4:144,157:]
        
    # PLOT MAPS FOR EACH INSTANT ----------------------------------------------
    
    for t in range(0,compensation_stack.shape[0]):
        # Extract maps from the stack
        DoD_map = DoD_maps[t,:,:,1]
        compensation_map = compensation_stack[t,:,:,0]
        threshold_map = threshold_pure_stack[t,:,:,0]
        envBAA_map = envBAA_stack[t,:,:]
        
    
        # Convert 0 as np.nan to make it transparent
        DoD_map = np.where(DoD_map==0,np.nan,DoD_map)
        compensation_map = np.where(compensation_map==0,np.nan,compensation_map)
        threshold_map = np.where(threshold_map==0,np.nan,threshold_map)
        envBAA_map = np.where(envBAA_map==0,np.nan,envBAA_map)
        
        # Rescale the topographic informations
        DoD_map = np.repeat(DoD_map, 10, axis=1)
        compensation_map = np.repeat(compensation_map, 10, axis=1)
        threshold_map = np.repeat(threshold_map, 10, axis=1)
        
        
        
        # CREATE THE PLOT
        fig, axs = plt.subplots(3, 1, figsize=(15, 6))

        # Add a title for the entire chart
        fig.suptitle(set_name + ' - Run ' + str(t) +
                     '-' + str(t+1), fontsize=16)

        # Plot MAA, BAA and compensation map
        im1 = axs[0].imshow(envBAA_map, cmap=ListedColormap(['green']), alpha=0.2)
        # im1 = axs[0].imshow(DoD_map, cmap=ListedColormap(['darkgray']))
        im1 = axs[0].imshow(DoD_map, cmap='Greys', alpha = 0.5)
        im1 = axs[0].imshow(compensation_map, cmap=ListedColormap(['cyan']))
        
        axs[0].set_title('Morphological, Bedlaod active and Compensation map')

        # Plot the compensation map
        im2 = axs[1].imshow(compensation_map, cmap=ListedColormap(['cyan']))
        im2 = axs[1].imshow(DoD_map, cmap='RdBu', vmin=-25, vmax=25)
        axs[1].set_title('Morphological active and Compensation map')

        # Plot the DoD map and the threshold map
        im3 = axs[2].imshow(envBAA_map, cmap=ListedColormap(['green']), alpha=0.2)
        im3 = axs[2].imshow(DoD_map, cmap='RdBu', vmin=-25, vmax=25)
        # cbar3 = fig.colorbar(im3, ax=axs[2], orientation='horizontal', pad=0.1, fraction=0.05)
        # Add a colorbar for the third subplot anchored on the right, vertically
        # cbar = fig.colorbar(im3, ax=axs[2], location='right', pad=0.05)
        # cbar.set_label('Threshold Values')
        im3 = axs[2].imshow(threshold_map, cmap=ListedColormap(['fuchsia']))
        axs[2].set_title('Morphological activity and Threshold map')

        # # Add horizontal colorbars below each subplot
        # cbar1 = fig.colorbar(im1, ax=axs[0], orientation='horizontal', pad=0.1, fraction=0.05)
        # cbar2 = fig.colorbar(im2, ax=axs[1], orientation='horizontal', pad=0.1, fraction=0.05)
        # cbar3 = fig.colorbar(im3, ax=axs[2], orientation='horizontal', pad=0.1, fraction=0.05)

        # Add legends for each subplot
        legend_elements_0 = [
            plt.Line2D([0], [0], color='darkgray', lw=4, label='DoD Map'),
            plt.Line2D([0], [0], color='cyan', lw=4,
                       label='Compensation Map'),
            plt.Line2D([0], [0], color='green', alpha=0.2, lw=4, label='Env BAA Map')
        ]
        legend_elements_1 = [
            plt.Line2D([0], [0], color='cyan', lw=4, label='Compensation Map')
        ]
        legend_elements_2 = [
            plt.Line2D([0], [0], color='fuchsia', lw=4, label='Threshold Map'),
            plt.Line2D([0], [0], color='green', alpha=0.2, lw=4, label='Env BAA Map')
        ]

        # Position legends outside the plot area
        axs[0].legend(handles=legend_elements_0,
                      loc='center left', bbox_to_anchor=(1, 0.5))
        axs[1].legend(handles=legend_elements_1,
                      loc='center left', bbox_to_anchor=(1, 0.5))
        axs[2].legend(handles=legend_elements_2,
                      loc='center left', bbox_to_anchor=(1, 0.5))

        # Adjust layout
        # Adjust rect to make space for suptitle
        plt.tight_layout(pad=1.0, h_pad=1.0, rect=[0, 0, 0.85, 0.96])

        # Further adjust the spacing
        plt.subplots_adjust(hspace=0.4)

        plt.text(0.5, -0.3, f"Generated by: {os.path.basename(__file__)}",
                 transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)

        # Save the plot as a PDF file
        plt.savefig(os.path.join(home_dir, 'output_data', set_name + '_' + str(t) +
                    '_comp_thrs_DoD_envDoD_overlapping.pdf'), dpi=800, bbox_inches='tight')

        plt.show()

    
    # # =============================================================================
    # # TEST MAPS
    # # =============================================================================
    
    # dod10 = stack_bool[0,:,:,0]
    # dod21 = stack_bool[1,:,:,0]
    # dod20 = stack_bool[0,:,:,1]
    # pure_comp = compensation_pure_stack[0,:,:,0]
    # comp = compensation_stack[0,:,:,0]
    # thrs = threshold_pure_stack[0,:,:,0]
    
    # dod20 = np.repeat(dod20, 10, axis=1)
    # comp = np.repeat(comp, 10, axis=1)
    # thrs = np.repeat(thrs, 10, axis=1)
    

