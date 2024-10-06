#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 15:44:32 2023

@author: erri
"""

# Import packages
import numpy as np
# import matplotlib.pyplot as plt
# import time
import os
# import cv2
# from skimage import morphology
# from scipy import ndimage
# from scipy.ndimage import convolve
from PIL import Image
from PiQs_BAW_func_v1 import *

# DEFINE FOLDERS -------------------------------------------------------------#
folder_home = os.getcwd() # Setup home folder 

# SELECT RUN -----------------------------------------------------------------#

# run_names = ['q07r1','q07r2','q07r3','q07r4','q07r5','q07r6','q07r7','q07r8','q07r9'
#         ,'q10r1','q10r2','q10r3','q10r4','q10r5','q10r6','q10r7','q10r8','q10r9'
#         ,'q15r1','q15r2','q15r3','q15r4','q15r5','q15r6','q15r7','q15r8','q15r9'
#         ,'q20r1','q20r2','q20r3','q20r4','q20r5','q20r6','q20r7','q20r8','q20r9']

# run_names = ['q07r1', 'q10r1', 'q15r1', 'q20r1']

# run_names = ['q07r1']

#run_names = ['q05r1','q05r2','q05r3','q05r4','q05r5','q05r6','q05r7','q05r8','q05r9']
# run_names = ['q07rgm4']

# run_names = ['q07rgm', 'q10rgm2', 'q15rgm2', 'q20rgm2']

#run_names = ['q07r1','q07r2','q07r3','q07r4','q07r5','q07r6','q07r7','q07r8','q07r9']
#run_names  = ['q10r1','q10r2','q10r3','q10r4','q10r5','q10r6','q10r7','q10r8','q10r9']
#run_names  = ['q15r1','q15r2','q15r3','q15r4','q15r5','q15r6','q15r7','q15r8','q15r9']
#run_names = ['q20r1','q20r2','q20r3','q20r4','q20r5','q20r6','q20r7','q20r8','q20r9']

#run_names = ['q05_1','q07_1','q10_2','q20_2']

run_names  = ['q07rgm','q10rgm2','q15rgm2','q20rgm2']
# runs = ['q15rgm2']
# runs = ['q07rgm']


# ANALYSIS PARAMETERS --------------------------------------------------------#
analysis_list = ['envelope_timescale']

# SCRIPT MODE ----------------------------------------------------------------#
plt_show = 0

# SCRIPT PARAMETERS ----------------------------------------------------------#
downsampling_dim = 5

# EXNER TIME ------------------------------------------------------------------
Txnr_array = [20, 30, 48, 74, 109]

# RUN TIME --------------------------------------------------------------------
run_time_array = [16, 21, 28, 47, 81]  # Run duration in minutes

# RUN DURATION ----------------------------------------------------------------
run_frames = [10, 15, 20, 40, 88]  # Run duration in number of frames

# EXNER TIME -----------------------------------------------------------------#
# env_tscale_array = [5,5,5,5]
# env_tscale_array = [1,1,1,1]
#env_tscale_array = [2,3,4,9,20]      #1/4 of the frame number (~1/8*Txnr)
#env_tscale_array = [3,4,6,9,13]        # 1/8 Txnr
# env_tscale_array = [5,7,12,18]      # 1/4 Txnr
#env_tscale_array = [10,15,26,37,54]    # 1/2 Txnr
# env_tscale_array = [15,22,36,55]    # 3/4 Txnr
#env_tscale_array = [20,30,48,74, 109]    # 1 Txnr
# env_tscale_array = [25,37,60,92]    #5/4 Txnr
#env_tscale_array = [30,45,72,110, 163]   #6/4 Txnr
# env_tscale_array = [35,52,84,129]   #7/4 Txnr
env_tscale_array = [40,59,96,147, 218] # 2 Txnr
# env_tscale_array = [46, 67, 108, 166]
# env_tscale_array = [51,74,120,184]
# env_tscale_array = [56,82,132,202]
# env_tscale_array = [61,89,144,221]
# env_tscale_array = [66,97,156,239]
#env_tscale_array = [80,120,179,296,436] #4 Txnr
#env_tscale_array = [89,135,180,333, 490] #4.5 Txnr

# =============================================================================
# LOOP OVER THE RUNS
# =============================================================================
for run_name in run_names:
    
    # DEFINE ENVELOPE TIMESCALE
    if '05' in run_name:
        skip = 1
        env_tscale = env_tscale_array[4]
        Txnr = Txnr_array[4]
        run_time = run_time_array[4]
        run_frm = run_frames[4]
        set_name = 'q05_1'
        rgm_code = 'q05rgm5'
        rgm_skip_frame = 0
        MAW_set_name = 'q05_1'

    if '07' in run_name:
        skip = 1
        env_tscale = env_tscale_array[3]
        Txnr = Txnr_array[3]
        run_time = run_time_array[3]
        run_frm = run_frames[3]
        set_name = 'q07_1'
        rgm_code = 'q07rgm'
        rgm_skip_frame = 0
        MAW_set_name = 'q07_1'

    if '10' in run_name:
        skip = 1
        env_tscale = env_tscale_array[2]
        Txnr = Txnr_array[2]
        run_time = run_time_array[2]
        run_frm = run_frames[2]
        set_name = 'q10_2'
        rgm_code = 'q10rgm2'
        rgm_skip_frame = 0
        MAW_set_name = 'q10_2'

    if '15' in run_name:
        skip = 1
        env_tscale = env_tscale_array[1]
        Txnr = Txnr_array[1]
        run_time = run_time_array[1]
        run_frm = run_frames[1]
        set_name = 'q15_2'
        rgm_code = 'q15rgm2'
        rgm_skip_frame = 5
        MAW_set_name = 'q15_3'

    if '20' in run_name:
        skip = 1
        env_tscale = env_tscale_array[0]
        Txnr = Txnr_array[0]
        run_time = run_time_array[0]
        run_frm = run_frames[0]
        set_name = 'q20_2'
        rgm_code = 'q20rgm2'
        rgm_skip_frame = 0
        MAW_set_name = 'q20_2'
        
    print('****************')
    print(run_name, '  Timescale: ', env_tscale)
    print('****************')
    
    
    # SET COUNTERS
    k=0
    m=0
    tscale_counter = 1
     
    # SETUP DATA FOLDER
    path_input_data = os.path.join(folder_home, 'output_data')

    path_folder_stacks = os.path.join(path_input_data, '3_PiQs_BAW_2Dt_filter',set_name)
    path_folder_envelopes = os.path.join(path_folder_stacks,'envelopes_cleaned','space_discr_' + str(downsampling_dim), run_name)
    path_partial_envelopes = os.path.join(path_folder_envelopes, run_name + '_envTscale' + str(env_tscale)) # Path where to stock the envelopes taken at a given timescale
    if not(os.path.exists(path_folder_envelopes)):
        os.makedirs(path_folder_envelopes)
    if not(os.path.exists(path_partial_envelopes)):
        os.mkdir(path_partial_envelopes)
    if not(os.path.exists(path_partial_envelopes)):
        os.mkdir(path_partial_envelopes)
    

    meanBAW_array=[]

    
    # LOAD THE TOTAL LOW RESOLUTION STACK
    stack_path = os.path.join(path_folder_stacks, run_name + '_BAA_stack_LR' + str(downsampling_dim) + '_cld.npy')
    
    # LOAD THE DATA STACK
    stack = np.load(stack_path)
    # CONVERT THE STACK TO BOOLEAN
    stack_bool = np.where(stack>0, 1, stack)
        
    
    
    for i in range(stack_bool.shape[0]): # For each images in the folder...
        
        BAA_map_LR = stack_bool[i,:,:]
        

        '''
        THIS SECTION WILL FILL STACKS OF BAA MAPS AT A GIVEN TIMESCALE
        NB: BAA_map_LR can be used instead of BAA_map to get low resolution stack or 
        In this section an envelope stack will be created every timescale.
        This procedure became necessary to compute the bedload migration rate from images.
        First because I would like to scale the migration rate calculation for different runs,
        second to have slots that are easy to handle an where to perform a few statistics.
        '''
        if 'envelope_timescale' in analysis_list:
            # print('check')
            # print('timescale counter: ', tscale_counter)

            if tscale_counter == 1:
                # DEFINE THE FIRST ITERATION STACK
                partial_stack = np.expand_dims(BAA_map_LR, axis=0)

                tscale_counter+=1 # Counter update

            elif tscale_counter%env_tscale!=0: # If the residual of the division is not equal to zero tscale_counter 
                # STACK ALL THE DATA
                partial_stack = np.vstack((partial_stack, np.expand_dims(BAA_map_LR, axis=0)))

                tscale_counter+=1 # Counter update

            elif tscale_counter%env_tscale==0:
                # FILL THE PARTIAL STACK
                partial_stack = np.vstack((partial_stack, np.expand_dims(BAA_map_LR, axis=0)))
                
                partial_stack_bool = np.where(partial_stack>0,1,0)
 
                # ----------------------------- MAKE THE ENVELOPE ---------------------------- #
                partial_envelope_sum = np.nansum(partial_stack_bool, axis=0) # As the sum of active times
                partial_envelope = np.nansum(partial_stack, axis=0)
                partial_envelope_bool = np.where(partial_envelope>0, 1, 0) # As a boolean map
                # ---------------------------------------------------------------------------- #
                
                
                # COMPUTE THE BEDLOAD ACTIVE WIDTH AND SAVE IT INTO A NUMPY ARRAY
                channel_width = 120
                meanBAW = np.nansum(partial_envelope_bool)/partial_envelope_bool.shape[1]/channel_width
                meanBAW_array = np.append(meanBAW_array, meanBAW)
                

                # # # TRIM THE CELLS THAT ARE ACTIVE LESS THAN 1/5 OF THE PARTAL STACK LENGTH
                # t_lenght_thrs = int(round(env_tscale*0.2))
                # partial_envelope_thrs_mask = np.where(partial_envelope_sum>t_lenght_thrs, 1, 0)
                partial_envelope_thrs = partial_envelope #*partial_envelope_thrs_mask
                partial_envelope_thrs_bool = np.where(partial_envelope_thrs>0, 1, 0) # As a boolean map

                # # -------------------- SAVE THE PARTIAL ENVELOPE AS IMAGE -------------------- #
                # partial_envelop_thrs_bool_img = Image.fromarray(np.array(partial_envelop_thrs_bool*1).astype(np.uint16))
                # partial_envelop_thrs_bool_img.save(os.path.join(path_partial_envelopes, run_name + '_'+ str(m)+'_envBAA_partial.tiff'))

                # ---------------------- SAVE THE BOOL ENVELOPE AS .npy ---------------------- #
                np.save(os.path.join(path_partial_envelopes, str(m) + '_' + run_name + '_partialBAA_cld_env.npy'), np.around(partial_envelope_thrs_bool, decimals=2))
                

                # ---------------------------------------------------------------------------- #

                # Update counter
                m +=1
                tscale_counter = 1
                
                
        k += 1
    
    
    # GIVEN ALL THE _partialBAA_cld_env.npy FILES, BUILD THE STACK FOR EACH RUN
    # --- INPUT DATA: A SERIES OF PARTIAL ENVELOPES TAKEN AT DIFFERENT TIMESCALE -- #
    if 'envelope_timescale' in analysis_list:
        
        '''
        The input file of this section are the _partialBAA_cld_env.npy files.
        These maps are the result of the envelope process of a series of
        partial stack take at a give timescale.
        The maps are already provided as boolean map and the envelope is still
        epurated from the pixels that are active less than 20% of the partial
        envelope length.
        '''
        # Get a list of all files in the folder
        file_list = os.listdir(path_partial_envelopes)
        # List the envelope in the partial envelopes folder:
        file_names = [file for file in file_list if file.endswith(run_name + '_partialBAA_cld_env.npy')]

        # Sort the list using the custom sorting key function
        file_names = sorted(file_names, key=custom_sort_key)

        for i, file_name in enumerate(file_names):
            path = os.path.join(path_partial_envelopes, file_name)
            partial_env_bool = np.load(path)

            
            # STACK ALL THE ENVELOPES IN A STACK
            if i == 0:
                stack_bool = np.expand_dims(partial_env_bool, axis=0)
            else:
                stack_bool = np.vstack((stack_bool,np.expand_dims(partial_env_bool, axis=0)))

        # SAVE THE PARTIAL STACK AS A BINARY FILE

        np.save(os.path.join(path_partial_envelopes, run_name+ '_envT' + str(env_tscale) + '_partial_envBAA_stack.npy'), stack_bool) # Save the stack
        np.save(os.path.join(path_partial_envelopes, run_name+ '_envT' + str(env_tscale) + '_partial_envBAA_sum_stack.npy'), np.nansum(stack_bool, axis=0)) # Save the stack
        
        test = np.nansum(stack_bool, axis=0)
        
        print('max: ', np.max(test))
    
    
    # Save the meanBAW_array:
    np.savetxt(os.path.join(path_partial_envelopes, run_name + '_' + str(env_tscale) + '_meanBAW.txt'), np.around(meanBAW_array, decimals=2))
    print(np.round(np.nanmean(meanBAW_array), decimals=3))
    print(np.round(np.nanstd(meanBAW_array), decimals=3))

    # Save the mean and std of all the partialBAA
    mean = np.round(np.mean(meanBAW_array), decimals=3)
    std = np.round(np.std(meanBAW_array), decimals=3)
    mean_std_stack = [mean, std]
    np.savetxt(os.path.join(path_partial_envelopes, run_name + '_envT' + str(env_tscale) + '_mean_std_BAW_allenv.txt'), mean_std_stack)
    


    