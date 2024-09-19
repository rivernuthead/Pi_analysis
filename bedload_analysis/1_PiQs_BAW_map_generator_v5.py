#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 17:54:30 2023

@author: erri
"""
# Import packages
import numpy as np
import matplotlib.pyplot as plt
import time
import os
import cv2
from skimage import morphology
from scipy import ndimage
from scipy.ndimage import convolve
from PIL import Image
from PiQs_BAW_func_v1 import *

# SELECT RUN ------------------------------------------------------------------
# runs = ['q07r1','q07r2','q07r3','q07r4','q07r5','q07r6','q07r7','q07r8','q07r9'
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

runs = ['q05r1','q05r2','q05r3','q05r4','q05r5','q05r6','q05r7','q05r8','q05r9']
# runs = ['q07r1','q07r2','q07r3','q07r4','q07r5','q07r6','q07r7','q07r8','q07r9']
# runs = ['q10r1','q10r2','q10r3','q10r4','q10r5','q10r6','q10r7','q10r8','q10r9']
# runs = ['q15r1','q15r2','q15r3','q15r4','q15r5','q15r6','q15r7','q15r8','q15r9']
# runs = ['q20r1','q20r2','q20r3','q20r4','q20r5','q20r6','q20r7','q20r8','q20r9']
# runs = ['q07rgm', 'q10rgm2', 'q15rgm2', 'q20rgm2']
# runs = ['q07rgm', 'q10rgm2', 'q15rgm2']

# runs = ['q20rgm2']
# runs = ['q07rgm']

# SCRIPT MODE -----------------------------------------------------------------
plt_show = 0

# SCRIPT PARAMETERS -----------------------------------------------------------
survey_length = 6.140 # [meters]
dt = 1 # Time between shoots in minutes


folder_home = os.getcwd() # Setup home folder 

# ----------------------------- INITIALIZE ARRAYS ---------------------------- #
BAW_segm_array = []

# for run, env_timescale in zip(runs, Txnr):

for run in runs:

    if run[1:3] == '05':
        skip=1
    if run[1:3] == '07':
        skip=1
    if run[1:3] == '10':
        skip=1
    if run[1:3] == '15':
        skip=1
    if run[1:3] == '20':
        skip=1

    # -------------------- LONG TO REGIME SPECIFIC PARAMETERS -------------------- #
    if 'rgm' in run:
        if run[1:3] == '05':
            mov_avg = 5
            eps = 0.4
            skip=0
        if run[1:3] == '07':
            mov_avg = 35
            eps = 0.2
            skip=0
        if run[1:3] == '10':
            mov_avg = 20
            eps = 0.2
            skip=0
        if run[1:3] == '15':
            mov_avg = 20
            eps = 0.2
            skip=0
        if run[1:3] == '20':
            mov_avg = 20
            eps = 0.2
            skip=0
        
    else:
         mov_avg = 5
         eps = 0.4
    
    # LOAD MASK ---------------------------------------------------------------
    mask_path = os.path.join(folder_home, 'mask.tif') # Define image path
    if '05' in run:
        mask_path = os.path.join(os.getcwd(), 'Border_mask_' + run + '.tif')
    mask = Image.open(mask_path) # Open image as image
    mask_arr = np.array(mask)
    # -------------------------------------------------------------------------

    # SET COUNTERS
    k=0
    j=1
    r=0
    m=0
    tscale_counter = 1
     
    # SETUP DATA FOLDER
    diff_path_out = os.path.join(folder_home, 'Maps')
    path_diff = folder_home + '/Differences/' + run # Set the directory path where to pick up images
    path_img = folder_home + '/Photos/' + run
    path_report = os.path.join(folder_home, 'output_report', run)

    path_blurry_area = folder_home + '/Blurry_areas/'
    
    # Check if the folders already exist and create them
    if not(os.path.exists(path_diff)):
        os.mkdir(path_diff)
    if not(os.path.exists(path_img)):
        os.mkdir(path_img)
    if not(os.path.exists(path_report)):
        os.mkdir(path_report)
    if not(os.path.exists(path_blurry_area)):
        os.mkdir(path_blurry_area)
    if not(os.path.exists(os.path.join(path_blurry_area, run))):
        os.mkdir(os.path.join(path_blurry_area, run))
    if not(os.path.exists(os.path.join(diff_path_out, run))):
        os.mkdir(os.path.join(diff_path_out, run))

    
    # Create a file list with all the diff name
    diff_name_files = sorted(os.listdir(path_diff))
    diff_names = []
    for name in diff_name_files:
        if name.endswith('.png') and not(name.endswith('rsz.png')):
            diff_names = np.append(diff_names, name)
        else:
            pass
        
    # Create a file list with all the photo name
    photo_name_files = sorted(os.listdir(path_img))
    photo_names = []
    for name in photo_name_files:
        if name.endswith('.jpg') and not(name.endswith('rsz.jpg')):
            photo_names = np.append(photo_names, name)
        else:
            pass
        
        
    diff_path = os.path.join(path_diff, diff_names[1]) # Define image path
    diff = Image.open(diff_path) # Open image as image

    diff_arr = np.array(diff) # Convert image as numpy array
    dim_y, dim_x = diff_arr.shape
    
    # # MASK DIFFERENCE
    mask_arr = mask_arr[:dim_y,:dim_x]
    diff_arr = diff_arr*mask_arr
    
    # Load UPSTREAM-DOVNSTREAM MASK
    path_img = os.path.join(folder_home, 'Photos', run)
    maskUD_path = os.path.join(path_img, 'Mask.tif') # Define image path
    maskUD = Image.open(maskUD_path) # Open image as image
    maskUD_arr = np.array(maskUD)
    maskUD_arr = np.where(maskUD_arr==255, 1, 0)
    maskUD_arr = maskUD_arr[:dim_y,:dim_x]
    
    
    BAW_report = []  # Actiwe width report
    envBAA_report = [] # Envelope active width report
    diff_names = diff_names[skip:]
    photo_names = photo_names[skip:]
    BAA_stack_bool = []
    
    active_map_stack=[]
    
    env_mgrt_rate = []
    
    act_mean_array=[]
    act_stdev_array=[]
    diff_mean_array=[]
    diff_stdev_array=[]
    env_disc_BAW_array=[]
    ist_mgrt_rate_array=[]
    ist_mgrt_rate_positive_array=[]
    ist_mgrt_rate_negative_array=[]
    ist_mgrt_rate_abs_array=[]
    ist_mgrt_rate_negative_rel_array=[]
    ist_mgrt_rate_positive_rel_array=[]
    ist_mgrt_rate_rel_array=[]
    ist_mgrt_rate_negative_abs_array=[]
    ist_mgrt_rate_positive_abs_array=[]
    
    matrix_BAW = np.zeros((len(diff_names),dim_x))

    
    print('****************')
    print(run)
    print('****************')
        
    #  CHECK OUTLIERS
    # THIS SECTION CHECK ALL THE IMAGES AND FIND THE IMAGES WITH FLASH AND OTHER DISTURBANCY. THIS SECTION PROVIDES A NUMPY ARRAY WHERE THE INDEX OF THE IMAGES OUT
    # OF A CERTAIN BOUNDARY IS STORED. 

    for name in diff_names: # For each images in the folder...

        diff_path = os.path.join(path_diff, name) # Define image path
        diff = Image.open(diff_path) # Open image as image
        shaded_area_path = os.path.join(os.getcwd(), 'shaded_area_images', run)
        if not os.path.exists(shaded_area_path):
            os.makedirs(shaded_area_path)
        diff_arr = np.array(diff) # Convert image in array
        diff_mean_array=np.append(diff_mean_array, np.mean(diff))
        diff_stdev_array=np.append(diff_stdev_array, np.std(diff))
        
        # # MASK DIFFERENCE
        diff_arr = diff_arr*mask_arr


    sig = moving_average_filter(diff_mean_array, mov_avg)


    r_samp = resample_array(diff_mean_array, sig)
    elements = np.where(abs(r_samp-diff_mean_array)>eps)
    print('Over saturated images: ', int(np.array(elements).shape[1]), 'over ',np.array(diff_names).shape[0], ' images.' )
    
    # Clean diff_mena_array: outliers set as np.nan
    diff_mean_arr_cleaned = np.where(abs(r_samp-diff_mean_array)>eps, np.nan, diff_mean_array)

    
    diff_mean_arr_cleaned_c = np.copy(diff_mean_arr_cleaned) # Create a new array where perform linear interpolation of np.nan
    diff_mean_arr_lin_interp = interpolate_nans(diff_mean_arr_cleaned_c) # Linear interpolated array

    
    
    if plt_show==1:
        plt.plot(np.linspace(0,len(diff_mean_array),len(diff_mean_array)), diff_mean_array, label='raw')
        # plt.plot(np.linspace(0,len(diff_mean_array),len(sig)), sig)
        plt.plot(np.linspace(0,len(diff_mean_array),len(sig)), sig+eps, '--', c='black')
        plt.plot(np.linspace(0,len(diff_mean_array),len(sig)), sig-eps, '--', c='black')
        plt.plot(np.linspace(0,len(r_samp),len(r_samp)), r_samp, label='resampling')
        plt.plot(np.linspace(0,len(diff_mean_arr_lin_interp),len(diff_mean_arr_lin_interp)), diff_mean_arr_lin_interp, '.', c='red', label='interp')
        plt.plot(np.linspace(0,len(diff_mean_arr_cleaned_c),len(diff_mean_arr_cleaned_c)), diff_mean_arr_cleaned_c, '--', label='cleaned')
        plt.title(run + ' - moving average: ' + str(mov_avg))
        plt.legend()

        plt.show(block=False)
        plt.close()
        # time.sleep(1)
        
    


# 
    for name in diff_names: # For each images in the folder...
        # print('****************')
        print(name)
        # print('****************')
        
        # 1. OPEN DIFFERENCES AND IMAGES, AND THEN CONVERT TO NP.ARRAY
        diff_path = os.path.join(path_diff, name) # Define image path
        diff = Image.open(diff_path) # Open image as image
        diff_arr = np.array(diff) # Convert image as numpy array
        
        # # MASK ARRAY DOMAIN
        diff_arr = diff_arr*mask_arr
        
        index = int(np.where(diff_names==name)[0]) # Define a index that mach the image name
        
        # PERFORM THE IMAGE CORRECTION
        if index in np.array(elements): # Check if image is outlier (elements array contains all the outliers name)
            diff_arr = diff_arr-(diff_mean_array[index]-diff_mean_arr_lin_interp[index]) # Correct the differece

        Image.fromarray(np.array(diff_arr)).save(os.path.join(diff_path_out, run, run + '_' + str(name)[:-4] + 'diff.tiff')) # Save image
        
        # DEFINE BLURRY AREA FOR EACH DIFFERENCES AS THE SUM OF THE BLURRY AREAS OF THE THREE IMAGES THAT MAKE THE DIFFERENCE
        if run[1:3] == '15' or run[1:3] == '20':   
            blurry_stack = np.zeros((3,diff_arr.shape[0], diff_arr.shape[1])) # Initialize stack where allocate the three images that define the difference
            
            # # For each difference, find the name of the related photos
            name1 = photo_names[index]
            name2 = photo_names[index+1]
            name3 = photo_names[index+2]
            
            # 2. LOOP TO FILL THE STACK with the three blurry areas detected images
            i=0
            for img_name in (name1, name2, name3):
                # print(img_name)
                blurry_img_path = os.path.join(folder_home, 'blurry_area_images', run, img_name[:-4] + '_blurry_thrs_fsh2.png')
                blurry = Image.open(blurry_img_path)
                # img_path = os.path.join(path_img , img_name)
                # img = Image.open(img_path)
                # new_resolution = (img.width//5, img.height//5) # reduce resolution
                # img = img.resize(new_resolution, Image.ANTIALIAS) # Resize the image to the new resolution
                blurry_arr = np.array(blurry)[:diff_arr.shape[0],:diff_arr.shape[1]]
                

                
                # Fill stack
                blurry_stack[i,:,:] = blurry_arr
                i +=1
               
            blurry_areas = np.sum(blurry_stack, axis=0)
            blurry_areas = 1*(blurry_areas>0) # Set blurry areas map as a boolean map
            blurry = blurry_arr*mask_arr # Mask considering the channel domain
            
            blurry_areas = np.where(blurry_areas==0, np.nan, blurry_areas) # Trim zero values
            cv2.imwrite(os.path.join(path_blurry_area, run, name+'_blurry_map.png'), blurry)
            cv2.imwrite(os.path.join(path_blurry_area, run, name+'_blurry_map.tiff'), blurry)

        # DEFINE THE THRESHOLS
        # Set the threshold:
        if '05' in run:
            # print('SKIPPING BLURRY ANALYSIS...')
            thrs_actDS = 6.5 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active # 6.5
            thrs_actUS = 6.5 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active # 6.4
            thrs_act = thrs_actDS
            shad_coeff = 2.0
            thrs_rso1 = 10000
            thrs_rso2 = 10000
        if run[1:3] == '07':
            # print('SKIPPING BLURRY ANALYSIS...')
            thrs_actDS = 6.5 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active # old 6.5    new 5.7
            thrs_actUS = 6.4 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active # old 6.4    new 5.7
            thrs_act = thrs_actDS
            shad_coeff = 2.0
            thrs_rso1 = 5000
            thrs_rso2 = 10000
        if run[1:3] == '10':
            # print('SKIPPING BLURRY ANALYSIS...')
            thrs_actDS = 5.3 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
            thrs_actUS = 5.8 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
            thrs_act = thrs_actDS
            shad_coeff = 2.0
            thrs_rso1 = 5000
            thrs_rso2 = 10000
        if run[1:3] == '15':
            thrs_actDS = 7 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
            thrs_actUS = 7 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
            thrs_act = thrs_actDS
            shad_coeff = 2.5
            thrs_rso1 = 5000
            thrs_rso2 = 20000
        if run[1:3] == '20':
            thrs_actDS = 7.0 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
            thrs_actUS = 5.7 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
            thrs_act = thrs_actDS
            shad_coeff = 1.0
            thrs_rso1 = 5000
            thrs_rso2 = 10000
        
        # THRESHOLDS FOR REGIME RUNS
        if 'rgm' in run:
            if '05' in run:
                # print('SKIPPING BLURRY ANALYSIS...')
                thrs_actDS = 6.5 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active # 6.5
                thrs_actUS = 6.5 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active # 6.4
                thrs_act = thrs_actDS
                shad_coeff = 2.0
                thrs_rso1 = 10000
                thrs_rso2 = 10000
            if run[1:3] == '07':
                # print('SKIPPING BLURRY ANALYSIS...')
                thrs_actDS = 6.5 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
                thrs_actUS = 6.4 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
                thrs_act = thrs_actDS
                shad_coeff = 2
                thrs_rso1 = 5000
                thrs_rso2 = 10000
            if run[1:3] == '10':
                # print('SKIPPING BLURRY ANALYSIS...')
                thrs_actDS = 6.4 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
                thrs_actUS = 6.0 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
                thrs_act = thrs_actDS
                shad_coeff = 3.0
                thrs_rso1 = 5000
                thrs_rso2 = 15000
            if run[1:3] == '15':
                thrs_actDS = 6.0 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
                thrs_actUS = 6.0 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
                thrs_act = thrs_actDS
                shad_coeff = 2.5
                thrs_rso1 = 5000
                thrs_rso2 = 20000
            if run[1:3] == '20':
                thrs_actDS = 7.0 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
                thrs_actUS = 5.7 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
                thrs_act = thrs_actDS
                shad_coeff = 1.0
                thrs_rso1 = 5000
                thrs_rso2 = 10000
        
        # Threshold sensitivity analysis
        thrs_actDS = thrs_actDS + 0.2
        thrs_actUS = thrs_actUS + 0.2

        # THIS SECTION PROVIDES CORRECTED THRESHOLDS FOR OUTLIERS IMAGES CASE BY CASE
        if run == 'q07r1':
            if name[3:7]=='0036':
                thrs_actDS = 6.5 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
                thrs_actUS = 7.5 # Active threshold - For values equal or grater than thrs_act bedeload is considered as active
            

        # 1. NOISE REMOVING
        shaded_area_imege_path = os.path.join(folder_home, 'shaded_area_images', run, run + '_shaded_area.tiff')
        
        banks_noise_bool = Image.open(shaded_area_imege_path)
        banks_noise_bool = np.array(banks_noise_bool)[:dim_y,:dim_x]
        banks_noise_bool = (banks_noise_bool>0)*shad_coeff
        
        # REMOVE BANKS NOSIE
        diff_arr_msk = diff_arr-banks_noise_bool*1  # Remove noise array  (Difference between integer)
        
        diff_arr_msk = diff_arr_msk*(diff_arr_msk>0) # Trim negative values
        diff_arr_msk = diff_arr_msk.astype(np.uint16)  # Convert image as uint16
        
        # Convert and save
        # diff_msk = Image.fromarray(diff_arr_msk.astype(np.uint16)) # Convert array to image
        # diff_msk.save(os.path.join(diff_path_out, run, run + '_' + str(name)[:-4] + '_noise_rm.tiff'))
        

        # 1. DIFF AVERAGING
        diff_arr_avg0 = cv2.GaussianBlur(diff_arr_msk, (11,11), 0)
        
        n_row0, n_col0 = 11, 11  # 11x11 is fine
        kernel = np.ones((n_row0, n_col0), np.float32)/(n_row0*n_col0)
        diff_arr_avg0 = convolve(diff_arr_avg0, kernel, mode='same')
        
        # n_row0, n_col0 =11,11  # 7,7 is fine
        # kernel0=np.ones((n_row0,n_col0), np.float32)/(n_row0*n_col0)
        # diff_arr_avg0=cv2.filter2D(src=diff_arr_msk,ddepth=-1, kernel=kernel0) 
        
        # Convert and save
        # diff_avg0=Image.fromarray(np.array(diff_arr_avg0))
        # diff_avg0.save(os.path.join(diff_path_out, run, run + '_' + str(name)[:-4] + '_avg0.tiff'))
        
        # TODO
        # 2. THRESHOLDING
        # diff_arr_thrs = np.where(diff_arr_avg0>=thrs_act, diff_arr_avg0,0)
        # th, diff_arr_thrs = cv2.threshold(diff_arr_avg0, 0, 255, cv2.THRESH_TOZERO + cv2.THRESH_OTSU)
        
        # APPLY THRESHOLD UPSTREAM
        th, diff_arr_thrsUS = cv2.threshold(diff_arr_avg0, thrs_actUS, 255, cv2.THRESH_TOZERO)
        diff_arr_thrsUS = diff_arr_thrsUS*maskUD_arr
        
        # APPLY THRESHOLD DOWNSTREAM
        th, diff_arr_thrsDS = cv2.threshold(diff_arr_avg0, thrs_actDS, 255, cv2.THRESH_TOZERO)
        diff_arr_thrsDS = diff_arr_thrsDS*abs(maskUD_arr-1)
        
        #MERGE THE TWO IMAGES
        diff_arr_thrs = diff_arr_thrsUS + diff_arr_thrsDS
        
        # Convert and save
        # diff_thrs = Image.fromarray(diff_arr_thrs.astype(np.uint16))
        # diff_thrs.save(os.path.join(diff_path_out, run, run + '_' + str(name)[:-4] + '_thrs.tiff'))
        
        # BLURRY AREAS ARE COMPUTED ONLY FOR q20_2 AND q15_2 RUN
        if run[1:3] == '15' or run[1:3] == '20':
            # 3. INCLUDE BLURRY ZONES
            # SET AS ACTIVE EQUAL TO THE ACTIVITY THRESHOLD BLURRY AREAS
            diff_arr_blur = np.where(np.logical_and(diff_arr_thrs==0, blurry_areas==1), thrs_act, diff_arr_thrs)
            diff_arr_blur = diff_arr_blur*mask_arr
            # diff_arr_blur = diff_arr_blur.astype(np.uint16)
            
            # Convert and save
            # diff_blur_img = Image.fromarray(np.array(diff_arr_blur).astype(np.uint16))
            # diff_blur_img.save(os.path.join(diff_path_out, run, run + '_' + str(name)[:-4] + '_raw_blurry.tiff'))
        else:
            
            diff_arr_blur = diff_arr_thrs
        
        # 4. MORPHOLOGICAL ANALYSIS
        # 4.a FILL SMALL HOLES
        diff_arr_rm_hls, diff_arr_rm_hls_target = fill_small_holes(
            matrix=diff_arr_blur, avg_target_kernel=51, area_threshold=1000,
            connectivity=1, value=10)
        
        # ERODE AREAS TO FIND PENINSULA
        diff_arr_rm_hls = diff_arr_rm_hls*(ndimage.binary_erosion(diff_arr_rm_hls, iterations=2))
        
        # 4.d REMOVE SMALL OBJECTS
        diff_arr_rsm_mask = morphology.remove_small_objects(diff_arr_rm_hls>0, min_size=thrs_rso1, connectivity=1)
        diff_arr_rsm = diff_arr_rm_hls*diff_arr_rsm_mask
        
        # Convert and save
        # diff_rsm = Image.fromarray(np.array(diff_arr_rsm).astype(np.uint16))
        # diff_rsm.save(os.path.join(diff_path_out, run, run + '_' + str(name)[:-4] + '_rm_obj.tiff'))
        
        
        # AVERAGE
        n_row0, n_col0 =21,51  # 7,7 is fine
        kernel0=np.ones((n_row0,n_col0), np.float32)/(n_row0*n_col0)
        diff_arr_avg1 = cv2.filter2D(src=diff_arr_rsm,ddepth=-1, kernel=kernel0)
        
        
        
        # REAPPLY THRESHOLD
        # APPLY THRESHOLD UPSTREAM
        th, diff_arr_thrs2US = cv2.threshold(diff_arr_avg1, thrs_actUS, 255, cv2.THRESH_TOZERO)
        diff_arr_thrs2US = diff_arr_thrs2US*maskUD_arr
        
        # APPLY THRESHOLD DOWNSTREAM
        th, diff_arr_thrs2DS = cv2.threshold(diff_arr_avg1, thrs_actDS, 255, cv2.THRESH_TOZERO)
        diff_arr_thrs2DS = diff_arr_thrs2DS*abs(maskUD_arr-1)
        
        #MERGE THE TWO IMAGES
        diff_arr_avg1 = diff_arr_thrs2US + diff_arr_thrs2DS
        
        diff_arr_avg1_mask = morphology.remove_small_objects(diff_arr_avg1>0, min_size=thrs_rso2, connectivity=1)
        diff_arr_avg1 = diff_arr_avg1*diff_arr_avg1_mask
        
        # Save the ultimate image as a numpy array
        np.save(os.path.join(diff_path_out, run, run + '_' + str(name)[:-4]+ '_ultimate_map.npy'), diff_arr_avg1)
        
        # PERFORM LINEAR DOWNSAMPLING TO OBTAIN THE LOW RESOLUTION VERSION
        BAA_map_LR5 = non_overlapping_average(diff_arr_avg1, kernel_size=5)
        BAA_map_LR10 = non_overlapping_average(diff_arr_avg1, kernel_size=10)
        
        maskUD_arr_LR5 = non_overlapping_average(maskUD_arr, kernel_size=5)
        maskUD_arr_LR5 = np.where(maskUD_arr_LR5>0.5,1,0)
        
        maskUD_arr_LR10 = non_overlapping_average(maskUD_arr, kernel_size=10)
        maskUD_arr_LR10 = np.where(maskUD_arr_LR10>0.5,1,0)
        
        # TODO
        # REAPPLY THRESHOLD 
        # APPLY THRESHOLD UPSTREAM
        th, BAA_map_LR5_US = cv2.threshold(BAA_map_LR5, thrs_actUS, 255, cv2.THRESH_TOZERO)
        BAA_map_LR5_US = BAA_map_LR5_US*maskUD_arr_LR5
        
        # APPLY THRESHOLD DOWNSTREAM
        th, BAA_map_LR5_DS = cv2.threshold(BAA_map_LR5, thrs_actDS, 255, cv2.THRESH_TOZERO)
        BAA_map_LR5_DS = BAA_map_LR5_DS*abs(maskUD_arr_LR5-1)
        
        #MERGE THE TWO IMAGES
        BAA_map_LR5 = BAA_map_LR5_US + BAA_map_LR5_DS
        
        

        if plt_show==1:
            plt.hist(BAA_map_LR5.flatten(), bins=100, density=True, alpha=0.5, color='blue', edgecolor='black')
    
            # Add labels and title
            plt.xlabel('Value')
            plt.ylabel('Density')
            plt.title('Histogram with Density')
            
            # Show plot
            plt.show()
                
        
        
        # REAPPLY THRESHOLD
        # APPLY THRESHOLD UPSTREAM
        th, BAA_map_LR10_US = cv2.threshold(BAA_map_LR10, thrs_actUS, 255, cv2.THRESH_TOZERO)
        BAA_map_LR10_US = BAA_map_LR10_US*maskUD_arr_LR10
        
        # APPLY THRESHOLD DOWNSTREAM
        th, BAA_map_LR10_DS = cv2.threshold(BAA_map_LR10, thrs_actDS, 255, cv2.THRESH_TOZERO)
        BAA_map_LR10_DS = BAA_map_LR10_DS*abs(maskUD_arr_LR10-1)
        
        #MERGE THE TWO IMAGES
        BAA_map_LR10 = BAA_map_LR10_US + BAA_map_LR10_DS
        
        # SAVE LOW RESOLUTION MAPS
        np.save(os.path.join(diff_path_out, run, run + '_' + str(name)
                [:-4] + '_ultimate_map_LR5.npy'), BAA_map_LR5)
        np.save(os.path.join(diff_path_out, run, run + '_' + str(name)
                [:-4] + '_ultimate_map_LR10.npy'), BAA_map_LR10)



        
        # Convert and save
        diff_avg1 = Image.fromarray(np.array(diff_arr_avg1).astype(np.uint16))
        diff_avg1.save(os.path.join(diff_path_out, run, run + '_' + str(name)[:-4]
                                    + '_cld_avg.tiff'))


        # STATISTICS ON DELTA SATURATION AND ACTIVITY SIGNAL
        act_mean = []
        act_stdev = []
        diff_mean = []
        diff_stdev = []
        
        
        # BAA_map IS GIVEN AS A BEDLOAD INTENSITY MAP!!
        BAA_map = np.copy(diff_arr_avg1)
        
        act_mean=np.append(act_mean, np.mean(BAA_map>0))
        act_stdev=np.append(act_stdev, np.std(BAA_map>0))
        diff_mean=np.append(diff_mean, np.mean(diff_arr))
        diff_stdev=np.append(diff_stdev, np.std(diff_arr))
        
        
        # print('Bedload activity map:')
        # print('Mean: ', np.mean(BAA_map))
        # print('Stdev: ', np.std(BAA_map))
        # DATA ANALYSIS
        # BEDLOAD ACTIVE WIDTH (BAW):
        array_actW = np.sum(BAA_map>0,axis=0)/np.nansum(mask_arr, axis=0) # Calculate the crosswise number of active pixel and divide it by the total number of pixels in the section
        matrix_BAW[k,:] = array_actW
        BAW_report = np.append(BAW_report, np.mean(array_actW)) # Append the average active width to build up the report

        if plt_show ==1:
            # PLOTS
            fig, ax = plt.subplots(dpi=300, tight_layout=True)
            yData = BAW_report
            xData = np.linspace(0,len(yData),len(yData))*dt
            actW = ax.scatter(xData, yData, s=0.8)
            ax.set_title(r'BAW* over time - ' + run)
            ax.set_xlabel('Time [min]')
            ax.set_ylabel('BAW* [-]')
            plt.savefig(os.path.join(path_report,run +'_' + str(skip) + '_BAW_plot_report.pdf'), dpi=1000) # raster (png, jpg, rgb, tif), vector (pdf, eps), latex (pgf)
            plt.show(block=False)
            plt.close()
            # time.sleep(1)
    
            fig, ax = plt.subplots(dpi=300, tight_layout=True)
            yData = envBAA_report
            xData = np.linspace(0,len(yData),len(yData))
            actW = ax.scatter(xData, yData, s=0.8)
            ax.set_title(r'BAW* over time - ' + run)
            ax.set_xlabel('Time [min]')
            ax.set_ylabel('envBAW* [-]')
            plt.savefig(os.path.join(path_report,run + '_' + str(skip) + '_envBAW_plot_report.pdf'), dpi=1000)
            plt.show(block=False)
            plt.close()
