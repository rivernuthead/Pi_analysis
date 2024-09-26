#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 17:54:30 2023

@author: erri
"""
"""
Description:
-------------
This script processes images from various runs to detect shaded and blurry areas using image processing techniques in Python.

Input data:
------------
    1. Fused images
    2. Difference of saturation maps
    3. Flume borders masks
    4. Upstream and dowstream image mask

Output data:
------------
    1. Blurry areas maps
    2. Shadow areas maps
    
"""
# Import packages
import numpy as np
import matplotlib.pyplot as plt
import os
import cv2
from skimage import morphology
from scipy import ndimage
from PIL import Image
from PiQs_BAW_func_v1 import *

#%%
# runs = ['q07r1'] # Run name  q07rgm, q10rgm2, q15rgm2, q20rgm2
# runs = ['q07rgm', 'q10rgm2cut', 'q15rgm2', 'q20rgm2']
# runs = ['q07rgm_test2', 'q10rgm2_test', 'q15rgm2_test', 'q20rgm2_test'] # TEST2 database
# runs = ['q07r1'] # Define the run names for batch process
runs = ['q05r1','q05r2','q05r3','q05r4','q05r5','q05r6','q05r7','q05r8','q05r9']

# runs = ['q07r1','q07r2','q07r3','q07r4','q07r5','q07r6','q07r7','q07r8','q07r9'
#         ,'q10r1','q10r2','q10r3','q10r4','q10r5','q10r6','q10r7','q10r8','q10r9'
#         ,'q15r1','q15r2','q15r3','q15r4','q15r5','q15r6','q15r7','q15r8','q15r9'
#         ,'q20r1','q20r2','q20r3','q20r4','q20r5','q20r6','q20r7','q20r8','q20r9']
# runs = ['q07r1', 'q10r1', 'q15r1', 'q20r9']
# runs = ['q07rgm', 'q10rgm2', 'q15rgm2', 'q20rgm2']
# runs = ['q20r9']
# runs = ['q10rgm2']


# SET MAIN FOLDERS
folder_home = os.getcwd()

for run in runs:
    print(run)
    # Setup data folder
    path_input_data = os.path.join(folder_home, 'input_data')
    path_output_data = os.path.join(folder_home, 'output_data')
    path_diff = os.path.join(path_input_data, '2_Differences', run[0:3], run) # Set the directory path where to pick up images
    path_img = os.path.join(path_input_data, '1_Fused_images', run[0:3], run)
    path_blurry_area = os.path.join(path_output_data, '0_blurry_shaded_areas_detection','blurry_area_images')
    path_shaded_imgs = os.path.join(path_output_data, '0_blurry_shaded_areas_detection','shaded_area_images')
    
    
    # LOAD borders image MASK
    mask_path = os.path.join(path_input_data,'masks', 'mask.tif') # Define image path
    if '05' in run:
        mask_path = os.path.join(path_input_data,'masks', 'Border_mask_' + run + '.tif')
    mask = Image.open(mask_path) # Open image as image
    mask_arr = np.array(mask)
    
    if run == 'q20rgm2':
        mask_arr = mask_arr[:,:6280]
    dim_y,dim_x = mask_arr.shape
    
    # Load UPSTREAM-DOVNSTREAM MASK
    maskUD_path = os.path.join(path_img, 'Mask.tif') # Define image path
    maskUD = Image.open(maskUD_path) # Open image as image
    maskUD_arr = np.array(maskUD)
    maskUD_arr = np.where(maskUD_arr==255, 1, 0)

    # Check if the folders already exist and create them
    if not(os.path.exists(path_img)):
        os.makedirs(path_img)
    if not(os.path.exists(path_shaded_imgs)):
        os.makedirs(path_shaded_imgs)
    if not(os.path.exists(os.path.join(path_shaded_imgs, run))):
        os.makedirs(os.path.join(path_shaded_imgs, run))
    if not(os.path.exists(path_blurry_area)):
        os.makedirs(path_blurry_area)
    if not(os.path.exists(os.path.join(path_blurry_area, run))):
        os.makedirs(os.path.join(path_blurry_area, run))
    
    # CREATE A FILE LIST WITH ALL THE DIFF NAME
    diff_name_files = sorted(os.listdir(path_diff))
    diff_names = []
    for name in diff_name_files:
        if name.endswith('.png') and not(name.endswith('rsz.png')):
            diff_names = np.append(diff_names, name)
        else:
            pass
        
    # CREATE A FILE LIST WITH ALL THE PHOTO NAME
    photo_name_files = sorted(os.listdir(path_img))
    photo_names = []
    for name in photo_name_files:
        if name.endswith('.jpg') and not(name.endswith('rsz.jpg')):
            photo_names = np.append(photo_names, name)
        else:
            pass
        


#%%
    print('COMPUTE SHADED AREAS')
    print()
    for i in range(0,len(photo_names[:-2])): # For each images in the folder...
        
        print('****************')
        print(i)
        print('****************')
        
        # For each difference, find the name of the related photos
        name1 = photo_names[i]
        name2 = photo_names[i+1]
        name3 = photo_names[i+2]
        
        noise_stack = np.zeros((3,dim_y, dim_x)) # Initialize stack where allocate the three images that define the difference
        
        # 2. LOOP TO FILL THE STACK with the three images that have generated the difference
        # These will be the basis for the shaded area detection and then removal
        m=0
        for img_name in (name1, name2, name3):
            print(img_name)
            img_path = os.path.join(path_img , img_name)
            img = Image.open(img_path)
            img_arr = np.array(img)[:dim_y,:dim_x,:]
            
            # BANKS NOISE REMOVAL
            img_arr_value = img_arr[:,:,0] # Keep only RED band
            img_arr_value = img_arr_value.astype(np.uint8)
            
            # Fill stack
            noise_stack[m,:,:] = img_arr_value
            m+=1


        # STACK ANALYSIS
        img_arr_value = np.mean(noise_stack, axis=0) # Compute the mean of the stack, now is a 2D matrix
        img_arr_value = img_arr_value.astype(np.uint8) # Conversion
        blur = cv2.GaussianBlur(img_arr_value, (11,11), 0) # Create the blur image
        th, banks_noise = cv2.threshold(blur, 0, 255, cv2.THRESH_TOZERO_INV + cv2.THRESH_OTSU) # Compute thresholding analysis
        # banks_noise = banks_noise*(banks_noise<=105)
        banks_noise, matrix_target = fill_small_holes(banks_noise, 10, 1000, 2)
        
        iterations = 1
        banks_noise = banks_noise*(ndimage.binary_erosion(banks_noise, iterations=iterations))
        
        
        banks_noise_bool = np.array(banks_noise, dtype='bool') # Convert image in bool
        banks_noise_bool = morphology.remove_small_objects(banks_noise_bool, min_size = 10000, connectivity=1) # Morphological analysis
        banks_noise = banks_noise*banks_noise_bool # Apply morphological analysis
        banks_noise_bool = banks_noise>0
        
        
        banks_noise = banks_noise*mask_arr # Apply domain mask
        banks_noise_bool = banks_noise_bool*mask_arr # Apply domain mask
        
        # AVERAGE
        n_row0, n_col0 =5,101  # 7,7 is fine
        kernel0=np.ones((n_row0,n_col0), np.float32)/(n_row0*n_col0)
        banks_noise_bool = cv2.filter2D(src=banks_noise_bool,ddepth=-1, kernel=kernel0)
        banks_noise = cv2.filter2D(src=banks_noise,ddepth=-1, kernel=kernel0)
        
        # REMOVE SMALL OBJECTS
        banks_noise_bool_mask = np.array(banks_noise, dtype='bool') # Convert image in bool
        banks_noise_bool_rso = morphology.remove_small_objects(
            banks_noise_bool_mask, min_size=100, connectivity=1)  # Morphological analysis
        banks_noise = banks_noise*banks_noise_bool_rso
        
        
        Image.fromarray(banks_noise).save(os.path.join(path_shaded_imgs, run,run + '_shaded_area.tiff')) # Save image
        Image.fromarray(banks_noise_bool).save(os.path.join(path_shaded_imgs, run,run + '_shaded_bool_area.tiff')) # Save image


#%%
    print() # Space between two different part
    print('COMPUTE BLURRY AREAS')
    # TODO
    # BLURRY AREA DETECTION
    
    if run[1:3] == '20' or run[1:3] == '15':
        if run[1:3] == '20':
            thrs_US = 25
            thrs_DS = 100 
            rso_min_size = 100000
        if run[1:3] == '15':
            thrs_US = 18
            thrs_DS = 18
            rso_min_size = 1000
            
        # photo_names = ['Img0013.jpg'] # Single image analysis
        
        for i in range(0,len(photo_names)): # For each images in the folder...
            name = photo_names[i] # Define the image name
    
            print(name)
            img_path = os.path.join(path_img , name) # Define the image path
    
            maskUD_arr = maskUD_arr[:dim_y,:dim_x]
            blurry = detect_blurry_regions(img_path)
            blurry = blurry[:dim_y,:dim_x]
            blurry = blurry*mask_arr # Domain mask
            cv2.imwrite(os.path.join(path_blurry_area, run, name[:-4] + '_blurry.tiff'), blurry)
    
            blurry_avg, kernel = neigh_filter(blurry, 13, 'variance') # Apply neighborhood analysis according to PiQs_BAW_func_v1
            # cv2.imwrite(os.path.join(path_blurry_area, run, name[:-4] + '_laplacian_var.tiff'), laplacian_var) # Save the blurry image in the original size
            blurry_avg = blurry_avg*mask_arr # Mask the blurry areas
            cv2.imwrite(os.path.join(path_blurry_area, run, name[:-4] + '_blurry_avg.tiff'), blurry_avg)
    
            # TREAT DIFFERENTLY UPSTREAM(=1) AND DOWNSTREAM(=0) CHANNEL PORTION
            mask_arrUS = maskUD_arr
            mask_arrDS = (abs(maskUD_arr-1))
            blurry_US = blurry_avg*maskUD_arr
            blurry_DS = blurry_avg*(abs(maskUD_arr-1))
            
            blurry_US = blurry_US*mask_arr
            blurry_DS = blurry_DS*mask_arr
            
            cv2.imwrite(os.path.join(path_blurry_area, run, name[:-4] + '_blurryUS.tiff'), blurry_US)
            cv2.imwrite(os.path.join(path_blurry_area, run, name[:-4] + '_blurryDS.tiff'), blurry_DS)
    
            
            # APPLY THE THRESHOLD
            blurry_US_bool = np.where(blurry_US<thrs_US, 1, 0)
            blurry_US_bool = blurry_US_bool*mask_arrUS
            blurry_US_bool = blurry_US_bool*mask_arr
            blurry_DS_bool = np.where(blurry_DS<thrs_DS, 1, 0)
            blurry_DS_bool = blurry_DS_bool*mask_arrDS
            blurry_DS_bool = blurry_DS_bool*mask_arr
            
            blurry_thrs = blurry_US_bool + blurry_DS_bool # Set threshold
            
            blurry_thrs = blurry_thrs*mask_arr
            
            cv2.imwrite(os.path.join(path_blurry_area, run, name[:-4] + '_blurry_thrs.tiff'), blurry_thrs)
            
            # EROSION PROCESS
            iterations = 3
            blurry_thrs_smooth_bool = ndimage.binary_erosion(
                blurry_thrs, iterations=iterations)
            blurry_thrs_smooth = blurry_thrs*blurry_thrs_smooth_bool
            cv2.imwrite(os.path.join(path_blurry_area, run,
                        name[:-4] + '_blurry_thrs_smooth.tiff'), blurry_thrs_smooth)
    
            # REMOVE SMALL HOLES
            blurry_thrs_fsh, blurry_thrs_target = fill_small_holes(
                blurry_thrs_smooth, 51, 500, 1)
            cv2.imwrite(os.path.join(path_blurry_area, run,
                        name[:-4] + '_blurry_thrs_fsh.tiff'), blurry_thrs_fsh)
            
            # AVERAGE
            n_row0, n_col0 =25,151  # 7,7 is fine
            kernel0=np.ones((n_row0,n_col0), np.float32)/(n_row0*n_col0)
            blurry_thrs_fsh_avg = cv2.filter2D(src=blurry_thrs_fsh,ddepth=-1, kernel=kernel0)
            cv2.imwrite(os.path.join(path_blurry_area, run,
                        name[:-4] + '_blurry_thrs_fsh_avg.tiff'), blurry_thrs_fsh_avg)
            
            # REMOVE SMALL OBJECTS
            blurry_thrs_bool_rso = np.array(blurry_thrs_fsh_avg, dtype='bool') # Convert image in bool
            blurry_thrs_bool_rso = morphology.remove_small_objects(
                blurry_thrs_bool_rso, min_size=rso_min_size, connectivity=1)  # Morphological analysis
            blurry_thrs_rso = blurry_thrs_fsh_avg*blurry_thrs_bool_rso
            cv2.imwrite(os.path.join(path_blurry_area, run,
                        name[:-4] + '_blurry_thrs_rso.tiff'), blurry_thrs_rso)
    
            # REMOVE SMALL HOLES
            blurry_thrs_fsh2, blurry_thrs_target = fill_small_holes(
                blurry_thrs_rso, 51, rso_min_size, 1, value=1)
            cv2.imwrite(os.path.join(path_blurry_area, run,
                        name[:-4] + '_blurry_thrs_fsh2.tiff'), blurry_thrs_fsh2)
            cv2.imwrite(os.path.join(path_blurry_area, run,
                        name[:-4] + '_blurry_thrs_fsh2.png'), blurry_thrs_fsh2)
            
            