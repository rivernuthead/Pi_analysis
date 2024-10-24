#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 10:17:19 2024

@author: erri
"""
# =============================================================================
# IMPORT LIBRARIES
# =============================================================================
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors

# =============================================================================
# SET FOLDERS
# =============================================================================
home_dir      = os.getcwd() # Home directory
report_dir    = os.path.join(home_dir, 'output')
output_folder = os.path.join(report_dir,'DEM_DoD_envBAA_overlapping')
plot_folder   = os.path.join(os.getcwd(), 'plots')
run_dir       = os.path.join(home_dir, 'surveys')

# =============================================================================
# SCRIPT PARAMETERS
# =============================================================================
run_mode = [
    'full_BAA_stack',
    # 'partial_envBAA_stack'
    ]
'''
run_mode.
        - full_BAA_stack: full BAA map stack
        - partial_envBAA_stack: stack computed as the stack of partial BAA envelopes (see /PiQs_analysis/4_PiQs_BAW_envelopes_generator.py for details)
'''

analysis_mode = ['activity_time_analysis']

# run_names = ['q07r1', 'q07r2', 'q07r3', 'q07r4', 'q07r5', 'q07r6', 'q07r7', 'q07r8', 'q07r9',
#             'q10r1', 'q10r2', 'q10r3', 'q10r4', 'q10r5', 'q10r6', 'q10r7', 'q10r8', 'q10r9',
#             'q15r1', 'q15r2', 'q15r3', 'q15r4', 'q15r5', 'q15r6', 'q15r7', 'q15r8', 'q15r9',
#             'q20r1', 'q20r2', 'q20r3', 'q20r4', 'q20r5', 'q20r6', 'q20r7', 'q20r8', 'q20r9']
# set_names = ['q07_1', 'q07_1', 'q07_1', 'q07_1', 'q07_1', 'q07_1', 'q07_1', 'q07_1', 'q07_1',
#             'q10_2', 'q10_2', 'q10_2', 'q10_2', 'q10_2', 'q10_2', 'q10_2', 'q10_2', 'q10_2',
#             'q15_2', 'q15_2', 'q15_2', 'q15_2', 'q15_2', 'q15_2', 'q15_2', 'q15_2', 'q15_2',
#             'q20_2', 'q20_2', 'q20_2', 'q20_2', 'q20_2', 'q20_2', 'q20_2', 'q20_2', 'q20_2']

run_names = ['q05r1', 'q05r2', 'q05r3', 'q05r4', 'q05r5', 'q05r6', 'q05r7', 'q05r8', 'q05r9']
set_names = ['q05_1', 'q05_1', 'q05_1', 'q05_1', 'q05_1', 'q05_1', 'q05_1', 'q05_1', 'q05_1']

run_names = ['q07r1', 'q07r2', 'q07r3', 'q07r4', 'q07r5', 'q07r6', 'q07r7', 'q07r8', 'q07r9']
set_names = ['q07_1', 'q07_1', 'q07_1', 'q07_1', 'q07_1', 'q07_1', 'q07_1', 'q07_1', 'q07_1']

run_names = ['q10r1', 'q10r2', 'q10r3', 'q10r4', 'q10r5', 'q10r6', 'q10r7', 'q10r8', 'q10r9']
set_names = ['q10_2', 'q10_2', 'q10_2', 'q10_2', 'q10_2', 'q10_2', 'q10_2', 'q10_2', 'q10_2']

run_names = ['q15r1', 'q15r3', 'q15r4', 'q15r5', 'q15r6', 'q15r8', 'q15r9']
set_names = ['q15_2', 'q15_2', 'q15_2', 'q15_2', 'q15_2', 'q15_2', 'q15_2']

run_names = ['q20r1', 'q20r2', 'q20r3', 'q20r4', 'q20r5', 'q20r6', 'q20r7', 'q20r8', 'q20r9']
set_names = ['q20_2', 'q20_2', 'q20_2', 'q20_2', 'q20_2', 'q20_2', 'q20_2', 'q20_2', 'q20_2']

# run_names = ['q05r1']
# set_names = ['q05_1']

# run_names = ['q07r1']
# set_names = ['q07_1']

# run_names = ['q10r1']
# set_names = ['q10_2']

# run_names = ['q15r1']
# set_names = ['q15_2']

# run_names = ['q20r1']
# set_names = ['q20_2']

# =============================================================================
# # INITIALIZE ARRAY FOR DISTRIBUTION ANALYSIS
# =============================================================================
morpho_active_array_list                  = []
bedload_active_array_list                 = []
travelling_bedload_array_list             = []
bedload_active_time_array_list            = []
bedload_morph_active_time_array_list      = []
bedload_travelling_active_time_array_list = []

DEM_morph_act_hist_list               = []
DEM_travelling_bedload_hist_list      = []
DEM_never_active_hist_list            = []
DEM_travelling_8th_bedload_hist_list  = []
DEM_travelling_full_bedload_hist_list = []

travelling_active_area_list = []
total_active_area_list      = []

hist_bedload_intensity_list      = []
hist_bedload_intensity_8th_list  = []
hist_bedload_intensity_full_list = []

values_bedload_intensity_list        = [] # List of bedload intensity values
values_bedload_intensity_8th_list    = [] # List of bedload intensity values for area active 1/8th of the exner time
values_bedload_intensity_full_list   = [] # List of bedload intensity values for area active entire run
values_bedload_intensity_morpho_list = [] # List of bedload intensity values for area where morphological changes occurs
values_bedload_intensity_comp_list   = [] # List of bedlaod intensity values for area where compensation effects occurs
values_bedload_intensity_thrs_list   = [] # List of bedload intensity values for area where under threshold effects occurs

# Pixel counting:
DEM              = []
BAA              = []
travelling_full  = []
travelling_8th   = []
travelling_total = []
MAA              = []


ratio1_list = []
ratio2_list = []
ratio3_list = []
ratio4_list = []
ratio5_list = []
ratio6_list = []
ratio7_list = []

# =============================================================================
# LOOP PVER RUNS
# =============================================================================
for run_name,set_name in zip(run_names, set_names):
    print(run_name + ' as part of ' + set_name + ' is running...')
    
    # =========================================================================
    # SETUP FOLDERS AND PATH
    # =========================================================================
    plot_dir_out = os.path.join(plot_folder, 'DEM_DoD_envBAA_overlapping', set_name)
    if not(os.path.exists(plot_dir_out)):
        os.mkdir(plot_dir_out)
    
    # IMPORT DEM DoD envBAA stack
    # Script that produces the input file:
    # /home/erri/Documents/PhD/Research/5_research_repos/DoD_analysis/6_DEM_DoD_envBAA_overlapping_v10.py
    if 'full_BAA_stack' in run_mode:
        path = os.path.join(output_folder, set_name + '_' + run_name + '_DEM_DoD_envBAA_stack.npy')
    elif 'partial_envBAA_stack' in run_mode:
        path = os.path.join(output_folder, set_name + '_' + run_name + '_DEM_DoD_partial_envBAA_stack.npy')
    
    DEM_DoD_envBAA_stack = np.load(path)
    
    
    #==========================================================================
    # COMPUTE METRICS
    # =========================================================================
    
    # CREATE THE MORPHOLOGICAL ACTIVITY MASK ----------------------------------
    MA_mask_nan = np.where(DEM_DoD_envBAA_stack[1,:,:]==0, np.nan, DEM_DoD_envBAA_stack[1,:,:]) # Morphological activity mask
    MA_mask_nan = np.where(abs(MA_mask_nan)>0,1,MA_mask_nan)
    MA_mask = np.where(np.isnan(MA_mask_nan), 0, MA_mask_nan) # Morphological activity mask
    # MA_mask_nan_inverse = np.where(abs(DEM_DoD_envBAA_stack[1,:,:])>0,np.nan,1)
    MA_mask_nan_inverse = np.where(np.isnan(DEM_DoD_envBAA_stack[1,:,:]), np.nan, np.where(DEM_DoD_envBAA_stack[1,:,:] == 0, 1, np.nan))
    
    # CREATE THE BEDLOAD ACTIVITY MASK ----------------------------------------
    BA_mask_nan = np.where(DEM_DoD_envBAA_stack[2,:,:]==0, np.nan, DEM_DoD_envBAA_stack[2,:,:]) # Bedload activity mask
    BA_mask_nan = np.where(abs(BA_mask_nan)>0,1,BA_mask_nan) # 
    BA_mask = np.where(np.isnan(BA_mask_nan), 0, BA_mask_nan) # Bedload activity mask
    # BA_mask_nan_inverse = np.where(abs(DEM_DoD_envBAA_stack[2,:,:])>0,np.nan,1)
    BA_mask_nan_inverse = np.where(np.isnan(DEM_DoD_envBAA_stack[2,:,:]), np.nan, np.where(DEM_DoD_envBAA_stack[2,:,:] == 0, 1, np.nan))
    
    # CREATE THE 1/8th RUN DURATION BEDLOAD ACTIVE MASK -----------------------
    mask_8th = np.where(DEM_DoD_envBAA_stack[4,:,:]==1, 1, np.nan)
    
    # CREATE THE FULL RUN ACTIVE BEDLOAD MASK ---------------------------------
    mask_full = np.where(DEM_DoD_envBAA_stack[4,:,:]==4, 1, np.nan)
    
    # ALL DEM VALUES ----------------------------------------------------------
    DEM_values = DEM_DoD_envBAA_stack[0,:,:].flatten()
    DEM_values = DEM_values[~np.isnan(DEM_values)]
    
    # DEM VALUES WHERE ONLY MORPHOLOGY IS ACTIVE ------------------------------
    DEM_morph_act_matrix = DEM_DoD_envBAA_stack[0,:,:]*MA_mask_nan
    DEM_morph_act = DEM_morph_act_matrix.flatten()
    DEM_morph_act = DEM_morph_act[~np.isnan(DEM_morph_act)]
    
    # DEM VALUES WHERE ONLY MORPHOLOGY IS INACTIVE ----------------------------
    DEM_morph_inact_matrix = DEM_DoD_envBAA_stack[0,:,:]*MA_mask_nan_inverse
    DEM_morph_inact = DEM_morph_inact_matrix.flatten()
    DEM_morph_inact = DEM_morph_inact[~np.isnan(DEM_morph_inact)]
    
    # DEM VALUES WHERE ONLY BEDLOAD IS ACTIVE ---------------------------------
    DEM_bedload_act_matrix = DEM_DoD_envBAA_stack[0,:,:]*BA_mask_nan
    DEM_bedload_act = DEM_bedload_act_matrix.flatten()
    DEM_bedload_act = DEM_bedload_act[~np.isnan(DEM_bedload_act)]
    total_active_area_list.append(np.nansum(BA_mask_nan))
    
    # DEM VALUES WHERE ONLY BEDLOAD IS INACTIVE -------------------------------
    DEM_bedload_inact_matrix = DEM_DoD_envBAA_stack[0,:,:]*BA_mask_nan_inverse
    DEM_bedload_inact = DEM_bedload_inact_matrix.flatten()
    DEM_bedload_inact = DEM_bedload_inact[~np.isnan(DEM_bedload_inact)]
    
    # DEM VALUES THAT ARE BEDLOAD ACTIVE BUT NO MORPHOLOGICALLY ACTIVE (NOT MORPHOLOGICAL BEDLOAD)
    DEM_travelling_bedload_matrix = DEM_DoD_envBAA_stack[0,:,:]*BA_mask_nan*MA_mask_nan_inverse
    DEM_travelling_bedload = DEM_travelling_bedload_matrix.flatten()
    DEM_travelling_bedload = DEM_travelling_bedload[~np.isnan(DEM_travelling_bedload)]
    travelling_active_area_list.append(np.nansum((DEM_DoD_envBAA_stack[2,:,:]*BA_mask_nan*MA_mask_nan_inverse)>0))
    
    # DEM VALUES MORPHOLOGICALLY ACTIVE BUT NOT BEDLOAD ACTIVE (PiQs fails)
    DEM_morpho_bedload_matrix = DEM_DoD_envBAA_stack[0,:,:]*BA_mask_nan_inverse*MA_mask_nan
    DEM_morpho_bedload = DEM_morpho_bedload_matrix.flatten()
    DEM_morpho_bedload = DEM_morpho_bedload[~np.isnan(DEM_morpho_bedload)]
    
    # DEM VALUE WHERE OR BEDLOAD OR THE MORPOLOGY ARE ACTIVE ------------------
    BA_MA_mask = abs(BA_mask) + abs(MA_mask)
    BA_MA_mask_nan = np.where(BA_MA_mask>0,1,np.nan) # Mask that is comprehensive of bedload and morphological activity
    BA_MA_mask = np.where(np.isnan(BA_MA_mask_nan),0,1)
    DEM_moprpho_plus_bedload = DEM_DoD_envBAA_stack[0,:,:]*BA_MA_mask
    DEM_moprpho_plus_bedload_arr = DEM_moprpho_plus_bedload.flatten()
    DEM_moprpho_plus_bedload_arr = np.where(DEM_moprpho_plus_bedload_arr==0, np.nan, DEM_moprpho_plus_bedload_arr)
    DEM_moprpho_plus_bedload_arr = DEM_moprpho_plus_bedload_arr[~np.isnan(DEM_moprpho_plus_bedload_arr)]
    
    # DEM VALUES THAT ARE NEVER ACTIVE ----------------------------------------
    never_active_mask = 1- np.where(BA_MA_mask>0,1,0)
    # Set all 0 values to np.nan
    never_active_mask = np.where(never_active_mask==0,np.nan,never_active_mask)
    DEM_never_active = DEM_DoD_envBAA_stack[0,:,:]*never_active_mask
    
    # DEM VALUES THAT ARE NON-MORPHOLOGICAL ACTIVE ONLY 1/8 OF THE RUN LENGTH
    DEM_travelling_8th_bedload_matrix = DEM_DoD_envBAA_stack[0,:,:]*BA_mask_nan*MA_mask_nan_inverse*mask_8th
    DEM_travelling_8th_bedload = DEM_travelling_8th_bedload_matrix.flatten()
    DEM_travelling_8th_bedload = DEM_travelling_8th_bedload[~np.isnan(DEM_travelling_8th_bedload)]
    
    # DEM VALUES THAT ARE NON-MORPHOLOGICAL ACTIVE FOR THE FULL RUN LENGTH
    DEM_travelling_full_bedload_matrix = DEM_DoD_envBAA_stack[0,:,:]*BA_mask_nan*MA_mask_nan_inverse*mask_full
    DEM_travelling_full_bedload = DEM_travelling_full_bedload_matrix.flatten()
    DEM_travelling_full_bedload = DEM_travelling_full_bedload[(~np.isnan(DEM_travelling_full_bedload))]
    
    # NON-MORPHOLOGICAL 1/8th RUN ACTIVE INTENSITY VALUES
    bedload_intensity_active_8th_matrix = DEM_DoD_envBAA_stack[3,:,:]*BA_mask_nan*MA_mask_nan_inverse*mask_8th
    bedload_intensity_active_8th = bedload_intensity_active_8th_matrix.flatten()
    bedload_intensity_active_8th = bedload_intensity_active_8th[(~np.isnan(bedload_intensity_active_8th)) & (bedload_intensity_active_8th!=0)]
    
    # NON-MORPHOLOGICAL FULL RUN ACTIVE INTENSITY VALUES
    bedload_intensity_active_full_matrix = DEM_DoD_envBAA_stack[3,:,:]*BA_mask_nan*MA_mask_nan_inverse*mask_full
    bedload_intensity_active_full = bedload_intensity_active_full_matrix.flatten()
    bedload_intensity_active_full = bedload_intensity_active_full[(~np.isnan(bedload_intensity_active_full)) & (bedload_intensity_active_full!=0)]
    
    # BEDLAOD INTENSITY VALUES THAT ARE MORPHOLOGICALLY ACTIVE ONLY -----------
    bedload_intensity_active_morpho_matrix = DEM_DoD_envBAA_stack[3,:,:]*MA_mask
    bedload_intensity_active_morpho = bedload_intensity_active_morpho_matrix[(~np.isnan(bedload_intensity_active_morpho_matrix)) & (bedload_intensity_active_morpho_matrix!=0)]
    
    
    # =========================================================================
    # COLLECT DATA FOR EACH SET_NAME
    # =========================================================================
    # TOTAL BEDLOAD INTENSITY
    bedload_intensity_values = DEM_DoD_envBAA_stack[3,:,:]
    bedload_intensity_values = bedload_intensity_values[(~np.isnan(bedload_intensity_values)) & (bedload_intensity_values!=0)]
    values_bedload_intensity_list.append(bedload_intensity_values)
    
    # BEDLOAD INTESITY FOR 1/8th RUN ACTIVE ---------------------------------------
    values_bedload_intensity_8th_list.append(bedload_intensity_active_8th)
    
    # BEDLOAD INTESITY FOR FULL RUN ACTIVE ---------------------------------------
    values_bedload_intensity_full_list.append(bedload_intensity_active_full)
    
    # BEDLOAD INTESITY FOR MORPHOLOGICAL ACTIVE ---------------------------------------
    values_bedload_intensity_morpho_list.append(bedload_intensity_active_morpho)
    
    
    # =========================================================================
    # COMPUTE HISTOGRAM
    # =========================================================================
    
    # HISTOGRAM PARAMETERS
    bins = 100
    
    # ALL DEM VALUES
    DEM_values_hist, bin_edges = np.histogram(DEM_values, bins=bins, range=(-50, 50), density=True)
    DEM_values_hist_norm = DEM_values_hist/np.nansum(DEM_values_hist)
    
    # DEM VALUES THAT ARE MORPHOLOGICAL ACTIVE
    DEM_morph_act_hist, bin_edges = np.histogram(DEM_morph_act, bins=bins, range=(-50, 50), density=True)
    DEM_morph_hist_norm = DEM_morph_act_hist/np.nansum(DEM_morph_act_hist)
    morpho_active_array_list.append(DEM_morph_hist_norm) # Store the distribution values in the distribution matrix. Append the current array to the list
    
    # DEM VALUES THAT ARE ACTIVE FOR BEDLOAD
    DEM_bedload_act_hist, bin_edges = np.histogram(DEM_bedload_act, bins=bins, range=(-50, 50), density=True)
    DEM_bedload_act_hist_norm = DEM_bedload_act_hist/np.nansum(DEM_bedload_act_hist)
    bedload_active_array_list.append(DEM_bedload_act_hist_norm) # Store the distribution values in the distribution matrix. Append the current array to the list
    
    # DEM VALUES FOR NON-MORPHOLOGICAL BEDLOAD
    DEM_travelling_bedload_hist, bin_edges = np.histogram(DEM_travelling_bedload, bins=bins, range=(-50, 50), density=True)
    DEM_travelling_bedload_hist_norm = DEM_travelling_bedload_hist/np.nansum(DEM_travelling_bedload_hist)
    travelling_bedload_array_list.append(DEM_travelling_bedload_hist_norm) # Store the distribution values in the distribution matrix. Append the current array to the list
    
    # DEM VALUES THAT ARE BOTH MORPHOLOGICAL AND BEDLOAD ACTIVE
    DEM_morpho_bedload_hist, bin_edges = np.histogram(DEM_morpho_bedload, bins=bins, range=(-50, 50), density=True)
    DEM_morpho_bedload_hist_norm = DEM_morpho_bedload_hist/np.nansum(DEM_morpho_bedload_hist)
    
    # DEM VALUES THAT ARE BEDLOAD OR MORPHOLOGICAL ACTIVE
    DEM_moprpho_plus_bedload_hist, bin_edges = np.histogram(DEM_moprpho_plus_bedload_arr, bins=bins, range=(-50, 50), density=True)
    DEM_moprpho_plus_bedload_hist_norm = DEM_moprpho_plus_bedload_hist/np.nansum(DEM_moprpho_plus_bedload_hist)
    
    # =========================================================================
    # PLOT THE HISTOGRAM
    # =========================================================================
    
    lower_limit_x, upper_limit_x = -60, 25  # Replace with desired x-axis limits
    lower_limit_y, upper_limit_y = 0, 0.1  # Replace with desired y-axis limits
    bin_midpoints = 0.5 * (bin_edges[1:] + bin_edges[:-1]) # Calculate midpoints of bins
    linewidth = 1
    
    # PLOT HISTOGRAM USING LINES
    plt.figure(figsize=(8, 4))
    # plt.plot(bin_midpoints, DEM_values_hist_norm, drawstyle='steps-post', color='#648FFF', linewidth=linewidth, label = 'DEM values')
    plt.plot(bin_midpoints, DEM_morph_hist_norm, drawstyle='steps-post', color='#427c31', linewidth=linewidth, label = 'Morphologically \n active')
    plt.plot(bin_midpoints, DEM_bedload_act_hist_norm, drawstyle='steps-post', color='#FE6100', linewidth=linewidth, label = 'Bedload \n active')
    plt.plot(bin_midpoints, DEM_travelling_bedload_hist_norm, drawstyle='steps-post', color='#124d93', linewidth=linewidth, label = 'Non-morphological \n bedload')
    # plt.plot(bin_midpoints, DEM_morpho_bedload_hist_norm, drawstyle='steps-post', color='#000000', linewidth=linewidth, label = 'Morph not \n bedload')
    # plt.plot(bin_midpoints, DEM_moprpho_plus_bedload_hist_norm, drawstyle='steps-post', color='#785EF0', linewidth=linewidth, label = 'Morph + \n Bedload areas')
    plt.ylim(bottom=0, top=0.25)
    plt.legend(fontsize='small')
    plt.axvline(x=0, color='black', linestyle='--',linewidth=1 , label='x=0')
    plt.xlabel('DEM values [mm]')
    plt.ylabel('Values')
    plt.title(run_name + ' - DEM values distribution ')
    plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}", transform=plt.gca().transAxes, ha='center', va='center', fontsize=10)
    plt.savefig(os.path.join(plot_dir_out, run_name + '_DEM_values_distribution_analysis.pdf'), bbox_inches='tight')
    plt.show()
    

    # =========================================================================
    # DEM VALUES DISTRIBUTION ANLYSIS WITH 7 DEM CLASSES
    # =========================================================================
    
    # DISTRIBUTION PARAMETERS -------------------------------------------------
    num_bins = 7
    DEM_bin_edges = [-60.0,-8.64,-3.86,-0.93,1.36,3.36,5.44,60.0]
    DEM_class_labels = [f"{DEM_bin_edges[i]},{DEM_bin_edges[i+1]}" for i in range(len(DEM_bin_edges) - 1)]
    hist_range = (-60, 60)  # Range of values to include in the histogram
    
    # DEM VALUES --------------------------------------------------------------
    DEM_hist_array = np.copy(DEM_DoD_envBAA_stack[0,:,:])
    DEM_hist_array = DEM_hist_array[~np.isnan(DEM_hist_array)] # Trim 0 and np.nan
    DEM_hist, DEM_bin_edges = np.histogram(DEM_hist_array, bins=DEM_bin_edges, range=hist_range)
    
    # MORPHOLOGICALLY ACTIVE AND INACTIVE CELLS -------------------------------
    DEM_morph_act_hist, DEM_morph_bin_edges = np.histogram(DEM_morph_act, bins=DEM_bin_edges, range=hist_range)
    DEM_morph_inact_hist, DEM_morph_bin_edges = np.histogram(DEM_morph_inact, bins=DEM_bin_edges, range=hist_range)

    # BEDLOAD ACTIVE AND INACTIVE CELLS ---------------------------------------
    DEM_bedload_act_hist, DEM_morph_bin_edges = np.histogram(DEM_bedload_act, bins=DEM_bin_edges, range=hist_range)
    DEM_bedload_inact_hist, DEM_morph_bin_edges = np.histogram(DEM_bedload_inact, bins=DEM_bin_edges, range=hist_range)
    
    # MORPHOLOGICALLY ACTIVE, NON-MORPHOLOGICAL BEDLOAD AND FULLY INACTIVE CELLS
    dist_bins = DEM_bin_edges
    DEM_morph_act_hist, DEM_morph_bin_edges = np.histogram(DEM_morph_act, bins=dist_bins, range=hist_range)
    DEM_travelling_bedload_hist , DEM_travelling_bedload_bin_edges = np.histogram(DEM_travelling_bedload, bins=dist_bins, range=hist_range)
    DEM_never_active_hist, DEM_never_active_bin_edges = np.histogram(DEM_never_active, bins=dist_bins, range=hist_range)
    data_total_sum = DEM_morph_act_hist + DEM_travelling_bedload_hist + DEM_never_active_hist # Histogram  values sum for normalization
    DEM_morph_bin_midpoints = 0.5 * (DEM_morph_bin_edges[1:] + DEM_morph_bin_edges[:-1]) # Calculate midpoints of bins
    
    # NON-MORPHOLOGICAL BEDLOAD ACTIVE 1/8 OF THE EXNER TIME ------------------
    DEM_travelling_8th_bedload_hist, DEM_travelling_8th_bedload_bin_edges = np.histogram(DEM_travelling_8th_bedload, bins=dist_bins, range=hist_range)
    
    # TRAVELLIN BEDLOAD ACTIVE FOR THE FULL RUN DURATION ----------------------
    DEM_travelling_full_bedload_hist, DEM_travelling_full_bedload_bin_edges = np.histogram(DEM_travelling_full_bedload, bins=dist_bins, range=hist_range)
    
    # =============================================================================
    # COMPUTE THE RATIO BETWEEN ACTIVE METRICS AND THE TOTAL DEM VALUES
    # =============================================================================
    DEM_bedload_act_hist_percent = DEM_bedload_act_hist/DEM_hist
    DEM_travelling_full_bedload_hist_percent = DEM_travelling_full_bedload_hist/DEM_hist
    DEM_travelling_8th_bedload_hist_percent = DEM_travelling_8th_bedload_hist/DEM_hist
    DEM_travelling_bedload_hist_percent = DEM_travelling_bedload_hist/DEM_hist
    DEM_morph_act_hist_percent = DEM_morph_act_hist/DEM_hist
    
    DEM.append(len(DEM_values))
    travelling_full.append(len(DEM_travelling_full_bedload))
    travelling_8th.append(len(DEM_travelling_8th_bedload))
    MAA.append(np.nansum(MA_mask_nan))
    
    ratio1 = DEM_travelling_full_bedload_hist_percent/DEM_travelling_bedload_hist_percent
    ratio2 = DEM_travelling_bedload_hist_percent/DEM_bedload_act_hist_percent
    ratio3 = DEM_travelling_full_bedload_hist_percent/DEM_bedload_act_hist_percent
    ratio4 = DEM_travelling_bedload_hist_percent/DEM_morph_act_hist_percent
    ratio5 = DEM_travelling_full_bedload_hist_percent/DEM_morph_act_hist_percent
    ratio6 = DEM_travelling_8th_bedload_hist_percent/DEM_bedload_act_hist_percent
    ratio7 = DEM_travelling_8th_bedload_hist_percent/DEM_travelling_bedload_hist_percent
    
    ratio1_list.append(ratio1)
    ratio2_list.append(ratio2)
    ratio3_list.append(ratio3)
    ratio4_list.append(ratio4)
    ratio5_list.append(ratio5)
    ratio6_list.append(ratio6)
    ratio7_list.append(ratio7)
    
    # Normalize data
    DEM_morph_act_hist               = DEM_morph_act_hist/data_total_sum
    DEM_travelling_bedload_hist      =  DEM_travelling_bedload_hist/data_total_sum
    DEM_never_active_hist            = DEM_never_active_hist/data_total_sum
    DEM_travelling_8th_bedload_hist  = DEM_travelling_8th_bedload_hist/data_total_sum
    DEM_travelling_full_bedload_hist = DEM_travelling_full_bedload_hist/data_total_sum
    
    # APPEND HIST ARRAY TO THE LIST
    DEM_morph_act_hist_list.append(DEM_morph_act_hist)
    DEM_travelling_bedload_hist_list.append(DEM_travelling_bedload_hist)
    DEM_never_active_hist_list.append(DEM_never_active_hist)

    DEM_travelling_8th_bedload_hist_list.append(DEM_travelling_8th_bedload_hist)
    DEM_travelling_full_bedload_hist_list.append(DEM_travelling_full_bedload_hist)

    # =========================================================================
    # PLOT DISTRIBUTION OF FREQUENCY
    # =========================================================================
    
    # PLOT NON-MORPHOLOGICAL FUL ACTIVE BEDLOAD VALUES ------------------------
    data = np.copy(DEM_travelling_full_bedload)
    hist_travelling, bins_travelling = np.histogram(data, bins=range(-60, 62))
    pdf_travelling = hist_travelling / np.sum(hist_travelling)
    travelling_total_area = np.sum(hist_travelling)
    plt.figure(figsize=(10, 6))
    plt.bar(bins_travelling[:-1], pdf_travelling, width=1, align='edge', color='orange', edgecolor='black', alpha=0.7)
    plt.title('"Non-morphological bedload" full run Probability Distribution Function - ' +  run_name)
    plt.xlabel('Value')
    plt.ylabel('Probability')
    plt.ylim(0, 0.25)
    plt.grid(True, alpha=0.5)
    plt.text(0.6, 0.95, f'Non-morphological bedload pixels: {travelling_total_area}', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes, fontsize=14)
    plt.text(0.5, -0.15, f"Generated by: {os.path.basename(__file__)}", transform=plt.gca().transAxes, ha='center', va='center', fontsize=10)
    plt.savefig(os.path.join(plot_dir_out, run_name + '_DEM_values_travelling_bedload_PDF.pdf'), bbox_inches='tight')
    plt.show()


    # PLOT TOTAL BEDLOAD INTENSITY --------------------------------------------
    data_bedload_intensity = np.copy(DEM_DoD_envBAA_stack[3,:,:])
    data_bedload_intensity = data_bedload_intensity[(~np.isnan(data_bedload_intensity)) & (data_bedload_intensity!=0)]
    hist_bedload_intensity, bins_bedload_intensity = np.histogram(data_bedload_intensity, bins=range(0, 30))
    pdf_bedload_intensity = hist_bedload_intensity / np.sum(hist_bedload_intensity)
    plt.figure(figsize=(10, 6))
    plt.bar(bins_bedload_intensity[:-1], pdf_bedload_intensity, width=1, align='edge', color='orange', edgecolor='black', alpha=0.7)
    plt.title('Bedload intensity Probability Distribution Function - ' +  run_name)
    plt.xlabel('Bedlaod intensity')
    plt.ylabel('Probability')
    plt.ylim(0, 0.25)
    plt.grid(True, alpha=0.5)
    plt.text(0.5, -0.15, f"Generated by: {os.path.basename(__file__)}", transform=plt.gca().transAxes, ha='center', va='center', fontsize=10)
    plt.savefig(os.path.join(plot_dir_out, run_name + '_bedload_intensity_PDF.pdf'), bbox_inches='tight')
    plt.show()


    # PLOT DISTRIBUTION OF FREQUENCY ------------------------------------------
    data_bedload_intensity_8th = np.copy(bedload_intensity_active_8th)
    data_bedload_intensity_full = np.copy(bedload_intensity_active_full)
    hist_bedload_intensity_8th, bins_bedload_intensity_8th = np.histogram(data_bedload_intensity_8th, bins=range(0, 30))
    pdf_bedload_intensity_8th = hist_bedload_intensity_8th / np.sum(hist_bedload_intensity)
    hist_bedload_intensity_full, bins_bedload_intensity_full = np.histogram(data_bedload_intensity_full, bins=range(0, 30))
    pdf_bedload_intensity_full = hist_bedload_intensity_full / np.sum(hist_bedload_intensity)
    plt.figure(figsize=(10, 6))
    plt.bar(bins_bedload_intensity[:-1], pdf_bedload_intensity, width=1, align='edge', color='orange', edgecolor='black', alpha=0.7, label='Bedload intensity')
    plt.bar(bins_bedload_intensity_8th[:-1], pdf_bedload_intensity_8th, width=1, align='edge', color='blue', edgecolor='black', alpha=0.5, label='1/8th active')
    plt.bar(bins_bedload_intensity_full[:-1], pdf_bedload_intensity_full, width=1, align='edge', color='green', edgecolor='black', alpha=0.4, label='Full run active')
    plt.title('Bedload intensity\nProbability Distribution Function - ' +  run_name)
    plt.xlabel('Bedload intensity')
    plt.ylabel('Probability')
    plt.ylim(0, 0.25)
    plt.grid(True, color='gray', alpha=0.5)
    plt.legend()
    max_value = np.round(pdf_bedload_intensity[0],decimals=2)
    # plt.text(0.32, 0.95, f'1st bar value: {max_value}', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes, fontsize=14)
    # plt.text(0.5, 0.85, f'Non-morphological bedload pixels: {travelling_total_area}', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes, fontsize=14)
    plt.text(0.5, -0.15, f"Generated by: {os.path.basename(__file__)}", transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.savefig(os.path.join(plot_dir_out, run_name + '_bedload_intensity_8th_and_full_PDF.pdf'), bbox_inches='tight')
    plt.show()
    
    # APPEND HISTOGRAM VALUES FOR REPORTS -------------------------------------
    hist_bedload_intensity_list.append(hist_bedload_intensity)
    hist_bedload_intensity_8th_list.append(hist_bedload_intensity_8th)
    hist_bedload_intensity_full_list.append(hist_bedload_intensity_full)


    
    
    if not(isinstance(dist_bins, int)):
        fig, ax = plt.subplots(figsize=(10, 6))
        # Plot the bars
        bars1 = ax.bar(range(len(DEM_hist)), DEM_morph_act_hist, label='Morphologically active')
        ax.bar(range(len(DEM_hist)), DEM_travelling_bedload_hist, label='Non-morphological bedload', bottom=DEM_morph_act_hist)
        ax.bar(range(len(DEM_hist)), DEM_never_active_hist, label='Fully inactive', bottom=DEM_travelling_bedload_hist+DEM_morph_act_hist)
        
        # Set x axes ticks
        plt.xticks(range(len(DEM_hist)), DEM_class_labels, rotation=45)  # Rotate x-axis labels for readability
        
        # # Add values above each bar
        # for x, y in zip(range(len(DEM_hist)), DEM_morph_act_hist/2):
        #     if np.logical_not(np.isnan(y)):
        #         plt.text(x, y, str(np.round(y, decimals=3)), ha='center', va='bottom')
        
        # # Add values above each bar
        # for x, y in zip(range(len(DEM_hist)), DEM_morph_act_hist + DEM_travelling_bedload_hist/2):
        #     if np.logical_not(np.isnan(y)):
        #         plt.text(x, y, str(np.round(y, decimals=3)), ha='center', va='bottom')
        
        # # Add values above each bar
        # for x, y in zip(range(len(DEM_hist)), DEM_morph_act_hist + DEM_travelling_bedload_hist + DEM_never_active_hist/2):
        #     if np.logical_not(np.isnan(y)):
        #         plt.text(x, y, str(np.round(y, decimals=3)), ha='center', va='bottom')
         
        # Set the y limit
        # plt.ylim(0,0.50)
        # plt.ylim(0,np.max(DEM_hist*1.1))
        plt.ylim(0,1.1)
        plt.xlabel('DEM class')
        plt.ylabel('Y Values')
        plt.title(run_name + ' - DEM values distribution \n morphological, non-morphological and inactive')
        plt.legend(fontsize='small')
        plt.text(0.5, -0.3, f"Generated by: {os.path.basename(__file__)}", transform=plt.gca().transAxes, ha='center', va='center', fontsize=10)
        plt.savefig(os.path.join(plot_dir_out, run_name + '_morpho_nonmorpho_inactive.pdf'), bbox_inches='tight')
        plt.show()
    
    

    # =========================================================================
    # # BEDLOAD ACTIVITY TIMES ANALYSIS
    # =========================================================================
    if 'activity_time_analysis' in analysis_mode:
        '''
        This section produces the arrays of the distribution of frequency of
        the number of times pixels are revealed as active
        '''
        
        # ALL THE ACTIVE PIXEL ------------------------------------------------
        bedload_activity_time_matrix = DEM_DoD_envBAA_stack[4,:,:]/np.nanmax(DEM_DoD_envBAA_stack[4,:,:])
        # bedload_activity_time_matrix = np.load(os.path.join(os.getcwd(), 'output', set_name + '_envBAA_act_4th_period_maps_single_runs.npy'))
        bedload_activity_time_array = bedload_activity_time_matrix[~np.isnan(bedload_activity_time_matrix) & (bedload_activity_time_matrix != 0)] # Trim zero values
        # Compute histogram
        bedload_activity_time_hist, bedload_activity_times_bin_edges = np.histogram(bedload_activity_time_array, bins=20, range=(0, np.nanmax(bedload_activity_time_matrix)), density=True)
        bedload_activity_time_hist_norm = bedload_activity_time_hist/np.nansum(bedload_activity_time_hist)
        bedload_act_times_bin_midpoints = bedload_activity_times_bin_edges[1:]
        bedload_active_time_array_list.append(bedload_activity_time_hist_norm) # Append values for report
        

        # BEDLOAD ACTIVITY WHERE MORPHOLOGICALLY ACTIVE ONLY ------------------
        bedload_morph_activity_time_matrix = DEM_DoD_envBAA_stack[4,:,:]*MA_mask_nan/np.nanmax(DEM_DoD_envBAA_stack[4,:,:]) # Apply morphological activity mask
        bedload_morph_activity_time_array = bedload_morph_activity_time_matrix[~np.isnan(bedload_morph_activity_time_matrix) & (bedload_morph_activity_time_matrix != 0)] # Trim np.nan values
        # Compute the histogram
        bedload_morph_activity_time_hist, bedload_morph_activity_times_bin_edges = np.histogram(bedload_morph_activity_time_array, bins=20, range=(0, np.nanmax(bedload_morph_activity_time_array)), density=True)
        bedload_morph_activity_time_hist_norm = bedload_morph_activity_time_hist/np.nansum(bedload_morph_activity_time_hist)
        bedload_morph_act_times_bin_midpoints = 0.5 * (bedload_morph_activity_times_bin_edges[1:] + bedload_morph_activity_times_bin_edges[:-1])
        bedload_morph_active_time_array_list.append(bedload_morph_activity_time_hist_norm) # Append values for report
        

        # NON-MORPHOLOGICAL BEDLOAD ACTIVE ONLY
        bedload_travelling_activity_time_matrix = DEM_DoD_envBAA_stack[4,:,:]*BA_mask_nan*MA_mask_nan_inverse/np.nanmax(DEM_DoD_envBAA_stack[4,:,:]) # Apply non-morphological bedload mask
        bedload_travelling_activity_time_array = bedload_travelling_activity_time_matrix[~np.isnan(bedload_travelling_activity_time_matrix) & (bedload_travelling_activity_time_matrix != 0)] # Trim zero values
        bedload_travelling_activity_time_array_norm = bedload_travelling_activity_time_array
        # Compute the histogram
        bedload_travelling_activity_time_hist, bedload_travelling_activity_times_bin_edges = np.histogram(bedload_travelling_activity_time_array, bins=20, range=(0, np.nanmax(bedload_travelling_activity_time_array)), density=True)
        bedload_travelling_activity_time_hist = bedload_travelling_activity_time_hist/np.nansum(bedload_travelling_activity_time_hist)
        bedload_travelling_act_times_bin_midpoints = 0.5 * (bedload_travelling_activity_times_bin_edges[1:] + bedload_travelling_activity_times_bin_edges[:-1])
        bedload_travelling_active_time_array_list.append(bedload_travelling_activity_time_hist) # Append values for report
        
        
    # =============================================================================
    # PLOT MAP OF BEDLOAD ACTIVITY, 1/8th AND FULL ACTIVE BEDLOAD ACTIVIRY
    # =============================================================================
    # Define a custom colormap with constant blue color
    cmap_blue = mcolors.LinearSegmentedColormap.from_list("", ["blue", "blue"])
    cmap_red = mcolors.LinearSegmentedColormap.from_list("", ["red", "red"])
    fig, axs = plt.subplots(figsize=(15, 5))
    plt.imshow(DEM_DoD_envBAA_stack[2,:,:], cmap='Greys', alpha=0.5, label = '1')
    plt.imshow(DEM_travelling_8th_bedload_matrix, cmap=cmap_red, label = '2')
    plt.imshow(DEM_travelling_full_bedload_matrix, cmap=cmap_blue, label = '3')
    plt.title(run_name + ' - Bedload activity time duration')
    legend_elements = [
        mpatches.Patch(color='gray', alpha=0.8, label='Bedlaod activity'),
        mpatches.Patch(color='red', alpha=0.8, label='Non-morphological 1/8 Txnr'),
        mpatches.Patch(color='blue', alpha=0.8, label='Non-morphological full run')
    ]
    plt.legend(handles=legend_elements, bbox_to_anchor=(1.001, 1), loc='upper left')
    plt.text(0.5, -0.25, f"Generated by: {os.path.basename(__file__)}", transform=plt.gca().transAxes, ha='center', va='center', fontsize=10)
    plt.savefig(os.path.join(plot_dir_out, run_name + '_travelling_bedload_map.pdf'))
    plt.show()

# =============================================================================
# SAVE REPORT MATRIX
# =============================================================================
# STACK THE INFORMATION FOR EACH set_name
DEM_morph_act_hist_matrix = np.vstack(DEM_morph_act_hist_list)
DEM_travelling_bedload_hist_matrix = np.vstack(
    DEM_travelling_bedload_hist_list)
DEM_never_active_hist_matrix = np.vstack(DEM_never_active_hist_list)
DEM_travelling_8th_bedload_hist_matrix = np.vstack(
    DEM_travelling_8th_bedload_hist_list)
DEM_travelling_full_bedload_hist_matrix = np.vstack(
    DEM_travelling_full_bedload_hist_list)
hist_bedload_intensity_matrix = np.vstack(hist_bedload_intensity_list)
hist_bedload_intensity_8th_matrix = np.vstack(hist_bedload_intensity_8th_list)
hist_bedload_intensity_full_matrix = np.vstack(
    hist_bedload_intensity_full_list)
# COMPUTE MEAN AND STANDARD DEVIATION
DEM_morph_act_hist_mean = np.vstack((np.nanmean(
    DEM_morph_act_hist_matrix, axis=0), np.nanstd(DEM_morph_act_hist_matrix, axis=0)))
DEM_travelling_bedload_hist_mean = np.vstack((np.nanmean(
    DEM_travelling_bedload_hist_matrix, axis=0), np.nanstd(DEM_travelling_bedload_hist_matrix, axis=0)))
DEM_never_active_hist_mean = np.vstack((np.nanmean(
    DEM_never_active_hist_matrix, axis=0), np.nanstd(DEM_never_active_hist_matrix, axis=0)))
DEM_travelling_8th_bedload_hist_mean = np.vstack((np.nanmean(
    DEM_travelling_8th_bedload_hist_matrix, axis=0), np.nanstd(DEM_travelling_8th_bedload_hist_matrix, axis=0)))
DEM_travelling_full_bedload_hist_mean = np.vstack((np.nanmean(
    DEM_travelling_full_bedload_hist_matrix, axis=0), np.nanstd(DEM_travelling_full_bedload_hist_matrix, axis=0)))
hist_bedload_intensity_hist_mean = np.vstack((np.nanmean(
    hist_bedload_intensity_matrix, axis=0), np.nanstd(hist_bedload_intensity_matrix, axis=0)))
hist_bedload_intensity_8th_hist_mean = np.vstack((np.nanmean(
    hist_bedload_intensity_8th_matrix, axis=0), np.nanstd(hist_bedload_intensity_8th_matrix, axis=0)))
hist_bedload_intensity_full_hist_mean = np.vstack((np.nanmean(
    hist_bedload_intensity_full_matrix, axis=0), np.nanstd(hist_bedload_intensity_full_matrix, axis=0)))
# SAVE INFORMATION ABOUT THE RATIO
ratio1_list_mean = np.vstack(
    (np.nanmean(ratio1_list, axis=0), np.nanstd(ratio1_list, axis=0)))
ratio2_list_mean = np.vstack(
    (np.nanmean(ratio2_list, axis=0), np.nanstd(ratio2_list, axis=0)))
ratio3_list_mean = np.vstack(
    (np.nanmean(ratio3_list, axis=0), np.nanstd(ratio3_list, axis=0)))
ratio4_list_mean = np.vstack(
    (np.nanmean(ratio4_list, axis=0), np.nanstd(ratio4_list, axis=0)))
ratio5_list_mean = np.vstack(
    (np.nanmean(ratio5_list, axis=0), np.nanstd(ratio5_list, axis=0)))
ratio6_list_mean = np.vstack(
    (np.nanmean(ratio6_list, axis=0), np.nanstd(ratio6_list, axis=0)))
ratio7_list_mean = np.vstack(
    (np.nanmean(ratio7_list, axis=0), np.nanstd(ratio7_list, axis=0)))
ratio1_list_comprehensive = np.array(
    travelling_full)/np.array(travelling_active_area_list)
ratio2_list_comprehensive = np.array(
    travelling_active_area_list)/np.array(total_active_area_list)
ratio3_list_comprehensive = np.array(
    travelling_full)/np.array(total_active_area_list)
ratio4_list_comprehensive = np.array(travelling_active_area_list)/np.array(MAA)
ratio5_list_comprehensive = np.array(travelling_full)/np.array(MAA)
ratio6_list_comprehensive = np.array(
    travelling_8th)/np.array(total_active_area_list)
ratio7_list_comprehensive = np.array(
    travelling_8th)/np.array(travelling_active_area_list)
ratio1_list_comprehensive_mean = np.vstack((np.nanmean(
    ratio1_list_comprehensive, axis=0), np.nanstd(ratio1_list_comprehensive, axis=0)))
ratio2_list_comprehensive_mean = np.vstack((np.nanmean(
    ratio2_list_comprehensive, axis=0), np.nanstd(ratio2_list_comprehensive, axis=0)))
ratio3_list_comprehensive_mean = np.vstack((np.nanmean(
    ratio3_list_comprehensive, axis=0), np.nanstd(ratio3_list_comprehensive, axis=0)))
ratio4_list_comprehensive_mean = np.vstack((np.nanmean(
    ratio4_list_comprehensive, axis=0), np.nanstd(ratio4_list_comprehensive, axis=0)))
ratio5_list_comprehensive_mean = np.vstack((np.nanmean(
    ratio5_list_comprehensive, axis=0), np.nanstd(ratio5_list_comprehensive, axis=0)))
ratio6_list_comprehensive_mean = np.vstack((np.nanmean(
    ratio6_list_comprehensive, axis=0), np.nanstd(ratio6_list_comprehensive, axis=0)))
ratio7_list_comprehensive_mean = np.vstack((np.nanmean(
    ratio7_list_comprehensive, axis=0), np.nanstd(ratio7_list_comprehensive, axis=0)))

# SAVE RATIO FILES IN .txt FORMAT
np.savetxt(os.path.join(plot_dir_out,  set_name + '_ratio1_list_mean.txt'),
           ratio1_list_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_ratio2_list_mean.txt'),
           ratio2_list_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_ratio3_list_mean.txt'),
           ratio3_list_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_ratio4_list_mean.txt'),
           ratio4_list_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_ratio5_list_mean.txt'),
           ratio5_list_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_ratio6_list_mean.txt'),
           ratio6_list_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_ratio7_list_mean.txt'),
           ratio7_list_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
# This is the coprehensive ratio, as the mean and the stdev for each run.
# It is a value, and it is not divided in DEM values classes
np.savetxt(os.path.join(plot_dir_out,  set_name + '_ratio1_list_comprehensive_mean.txt'),
           ratio1_list_comprehensive_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_ratio2_list_comprehensive_mean.txt'),
           ratio2_list_comprehensive_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_ratio3_list_comprehensive_mean.txt'),
           ratio3_list_comprehensive_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_ratio4_list_comprehensive_mean.txt'),
           ratio4_list_comprehensive_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_ratio5_list_comprehensive_mean.txt'),
           ratio5_list_comprehensive_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_ratio6_list_comprehensive_mean.txt'),
           ratio6_list_comprehensive_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_ratio7_list_comprehensive_mean.txt'),
           ratio7_list_comprehensive_mean, delimiter='\t', header=set_name + ' - mean values and stdev')


# =============================================================================
# SAVE MEAN AND STDEV VALUES AS .txt FILES
# =============================================================================
np.savetxt(os.path.join(plot_dir_out,  set_name + '_DEM_morph_act_hist_mean.txt'),
           DEM_morph_act_hist_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_DEM_travelling_bedload_hist_mean.txt'),
           DEM_travelling_bedload_hist_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_DEM_never_active_hist_mean.txt'),
           DEM_never_active_hist_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_DEM_travelling_8th_bedload_hist_mean.txt'),
           DEM_travelling_8th_bedload_hist_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_DEM_travelling_full_bedload_hist_mean.txt'),
           DEM_travelling_full_bedload_hist_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_travelling_active_area.txt'),
           travelling_active_area_list, delimiter='\t', header=set_name + ' - number of non-morphological active pixel')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_bedload_actve_area.txt'),
           total_active_area_list, delimiter='\t', header=set_name + ' - number of bedload active pixel')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_hist_bedload_intensity_hist_mean.txt'),
           hist_bedload_intensity_hist_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_hist_bedload_intensity_8th_hist_mean.txt'),
           hist_bedload_intensity_8th_hist_mean, delimiter='\t', header=set_name + ' - mean values and stdev')
np.savetxt(os.path.join(plot_dir_out,  set_name + '_hist_bedload_intensity_full_hist_mean.txt'),
           hist_bedload_intensity_full_hist_mean, delimiter='\t', header=set_name + ' - mean values and stdev')

# =============================================================================
# STACK AND SAVE BEDLOAD, MORPH-BEDLAOD, 1/8th AND FULL DURATION INTENSITY VALUES
# =============================================================================
values_bedload_intensity = np.hstack(values_bedload_intensity_list)
values_bedload_intensity_8th = np.hstack(values_bedload_intensity_8th_list)
values_bedload_intensity_full = np.hstack(values_bedload_intensity_full_list)
values_bedload_intensity_morpho = np.hstack(
    values_bedload_intensity_morpho_list)
np.save(os.path.join(plot_dir_out,  set_name +
        '_values_bedload_intensity_list.npy'), np.asanyarray(values_bedload_intensity))
np.save(os.path.join(plot_dir_out,  set_name + '_values_bedload_intensity_8th_list.npy'),
        np.asanyarray(values_bedload_intensity_8th))
np.save(os.path.join(plot_dir_out,  set_name + '_values_bedload_intensity_full_list.npy'),
        np.asanyarray(values_bedload_intensity_full))
np.save(os.path.join(plot_dir_out,  set_name + '_values_bedload_intensity_morpho_list.npy'),
        np.asanyarray(values_bedload_intensity_morpho))
# =============================================================================



# =============================================================================
# PLOT OF THE DISTRIBUTION OF FREQUENCY OF DEM VALUES
# =============================================================================
DEM_values_arrays_list = []  # Initialize an empty list to store arrays
lower_limit_x, upper_limit_x = -60, 25  # Replace with desired x-axis limits
lower_limit_y, upper_limit_y = 0, 0.1  # Replace with desired y-axis limits
# Calculate midpoints of bins
bin_midpoints = 0.5 * (bin_edges[1:] + bin_edges[:-1])
for run_name, set_name in zip(run_names, set_names):
    plot_dir_out = os.path.join(plot_folder, 'DEM_DoD_envBAA_overlapping', set_name)
    # IMPORT DEM DoD envBAA STACK FROM DEM_DoD_envBAA_overlapping_v2.py
    path = os.path.join(output_folder, set_name + '_' +
                        run_name + '_DEM_DoD_envBAA_stack.npy')
    DEM_DoD_envBAA_stack_history = np.load(path)
    DEM_values = DEM_DoD_envBAA_stack_history[0, :, :].flatten()
    DEM_values = DEM_values[~np.isnan(DEM_values)]
    bins = 100
    DEM_values_hist, bin_edges = np.histogram(
        DEM_values, bins=bins, range=(-50, 50), density=True)
    DEM_values_hist_norm = DEM_values_hist/np.nansum(DEM_values_hist)
    plt.plot(bin_midpoints, DEM_values_hist_norm,
             drawstyle='steps-post', linewidth=1, label=run_name)
    # Append the current array to the list
    DEM_values_arrays_list.append(DEM_values_hist_norm)
plt.ylim(bottom=0, top=0.15)
plt.legend(fontsize='small')
plt.axvline(x=0, color='black', linestyle='--', linewidth=1, label='x=0')
plt.xlabel('DEM values [mm]')
plt.ylabel('Values')
plt.title(set_name + ' - DEM values')
plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
         transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
plt.savefig(os.path.join(plot_dir_out , set_name +
            '_DEM_values_distribution.pdf'), bbox_inches='tight')
plt.show()

# COMPUTE MEAN AND STDEV, THEN STACK AND SAVE REPORT --------------------------
DEM_values_matrix_result = np.vstack(DEM_values_arrays_list)
DEM_values_hist_average = np.nanmean(DEM_values_matrix_result, axis=0)
DEM_values_hist_stdev = np.nanstd(DEM_values_matrix_result, axis=0)
output_data = np.vstack(
    (bin_midpoints, DEM_values_hist_average, DEM_values_hist_stdev))
header = 'Midpoints, average, stdev \n Generated by ' + \
    str(os.path.basename(__file__))
np.savetxt(os.path.join(plot_dir_out,  set_name + '_DEM_values_distribution_average.txt'),
           output_data, delimiter='\t', header=header)
# =============================================================================



# =============================================================================
# COMPUTE MEAN AND STEDV, THEN STACK AND SAVE REPORT FOR ALL THE HISTOGRAM
# =============================================================================

# MORPHOLOGICALLY ACTIVE DATA -------------------------------------------------
# Convert the list of arrays to a NumPy matrix
morph_active_values_matrix_result = np.vstack(morpho_active_array_list)
morph_active_hist_average = np.nanmean(
    morph_active_values_matrix_result, axis=0)  # Compute the histogram average
# Compute the histogram standard deviation
morph_active_hist_stdev = np.nanstd(morph_active_values_matrix_result, axis=0)
output_data = np.vstack(
    (bin_midpoints, morph_active_hist_average, morph_active_hist_stdev))
header = 'Midpoints, average, stdev \n Generated by ' + \
    str(os.path.basename(__file__))
np.savetxt(os.path.join(plot_dir_out,  set_name + 'DEM_morph_active_values_distribution_average.txt'),
           output_data, delimiter='\t', header=header)


# BEDLOAD ACTIVE DATA ---------------------------------------------------------
# Convert the list of arrays to a NumPy matrix
bedload_active_values_matrix_result = np.vstack(bedload_active_array_list)
bedload_active_hist_average = np.nanmean(
    bedload_active_values_matrix_result, axis=0)  # Compute the histogram average
# Compute the histogram standard deviation
bedload_active_hist_stdev = np.nanstd(
    bedload_active_values_matrix_result, axis=0)
output_data = np.vstack(
    (bin_midpoints, bedload_active_hist_average, bedload_active_hist_stdev))
header = 'Midpoints, average, stdev \n Generated by ' + \
    str(os.path.basename(__file__))
np.savetxt(os.path.join(plot_dir_out,  set_name + 'DEM_bedload_active_values_distribution_average.txt'),
           output_data, delimiter='\t', header=header)


# NON-MORPHOLOGICAL BEDLOAD DATA ----------------------------------------------
travelling_bedload_active_values_matrix_result = np.vstack(
    travelling_bedload_array_list)  # Convert the list of arrays to a NumPy matrix
travelling_bedload_active_hist_average = np.nanmean(
    travelling_bedload_active_values_matrix_result, axis=0)  # Compute the histogram average
travelling_bedload_active_hist_stdev = np.nanstd(
    travelling_bedload_active_values_matrix_result, axis=0)  # Compute the histogram standard deviation
output_data = np.vstack(
    (bin_midpoints, travelling_bedload_active_hist_average, travelling_bedload_active_hist_stdev))
header = 'Midpoints, average, stdev \n Generated by ' + \
    str(os.path.basename(__file__))
np.savetxt(os.path.join(plot_dir_out,  set_name + 'DEM_travelling_bedload_active_values_distribution_average.txt'),
           output_data, delimiter='\t', header=header)


# BEDLOAD ACTIVITY TIMES ------------------------------------------------------
if 'activity_time_analysis' in analysis_mode:

    # ALL DEM VALUES
    bedload_activity_times_matrix_result = np.vstack(
        bedload_active_time_array_list)  # Convert the list of arrays to a NumPy matrix
    bedload_activity_times_hist_average = np.nanmean(
        bedload_activity_times_matrix_result, axis=0)  # Compute the histogram average
    bedload_activity_times_hist_stdev = np.nanstd(
        bedload_activity_times_matrix_result, axis=0)  # Compute the histogram standard deviation
    output_data_bed_act = np.vstack(
        (bedload_act_times_bin_midpoints, bedload_activity_times_hist_average, bedload_activity_times_hist_stdev))
    header = 'Midpoints, average, stdev \n Generated by ' + \
        str(os.path.basename(__file__))
    np.savetxt(os.path.join(plot_dir_out,  set_name + '_bedload_active_times_values_distribution_average.txt'),
               output_data_bed_act, delimiter='\t', header=header, fmt='%.4f')

    # BEDLOAD AND MORPHOLOGICAL ACTIVITY
    bedload_morph_activity_times_matrix_result = np.vstack(
        bedload_morph_active_time_array_list)  # Convert the list of arrays to a NumPy matrix
    bedload_morph_activity_times_hist_average = np.nanmean(
        bedload_morph_activity_times_matrix_result, axis=0)  # Compute the histogram average
    bedload_morph_activity_times_hist_stdev = np.nanstd(
        bedload_morph_activity_times_matrix_result, axis=0)  # Compute the histogram standard deviatioN
    output_data_bedload_morph_act = np.vstack(
        (bedload_morph_act_times_bin_midpoints, bedload_morph_activity_times_hist_average, bedload_morph_activity_times_hist_stdev))
    header = 'Midpoints, average, stdev \n Generated by ' + \
        str(os.path.basename(__file__))
    np.savetxt(os.path.join(plot_dir_out,  set_name + '_bedload_morph_active_times_values_distribution_average.txt'),
               output_data_bedload_morph_act, delimiter='\t', header=header, fmt='%.4f')

    # NON MORPHOLOGICAL BEDLOAD ACTIVITY
    bedload_travelling_activity_times_matrix_result = np.vstack(
        bedload_travelling_active_time_array_list)  # Convert the list of arrays to a NumPy matrix
    bedload_travelling_activity_times_hist_average = np.nanmean(
        bedload_travelling_activity_times_matrix_result, axis=0)  # Compute the histogram average
    bedload_travelling_activity_times_hist_stdev = np.nanstd(
        bedload_travelling_activity_times_matrix_result, axis=0)  # Compute the histogram standard deviation
    output_data_bedload_travelling = np.vstack((bedload_travelling_act_times_bin_midpoints,
                                               bedload_travelling_activity_times_hist_average, bedload_travelling_activity_times_hist_stdev))
    header = 'Midpoints, average, stdev \n Generated by ' + \
        str(os.path.basename(__file__))
    np.savetxt(os.path.join(plot_dir_out,  set_name + '_bedload_travelling_active_times_values_distribution_average.txt'),
               output_data_bedload_travelling, delimiter='\t', header=header, fmt='%.4f')




#%%============================================================================
# BEDLOAD INTENSITY AND ACTIVE PERIODS FOR COMPENSATION AND UNDER THRESHOLD EFFECTS
# =============================================================================
home_dir = os.getcwd()  # Home directory
report_dir = os.path.join(home_dir, 'output')

set_names = ['q05_1', 'q07_1', 'q10_2', 'q15_2', 'q20_2']

def clean_data(arr):
    flat_arr = arr.flatten()
    return flat_arr[(~np.isnan(flat_arr)) & (flat_arr != 0)]

for set_name in set_names:
    print(set_name + ' is running...')

    # IMPORT DATA
    DoD_maps = np.load(os.path.join(report_dir, set_name + '_DoD_1Txnr_stack.npy'))
    envBAA_intensity_1Txnr_stack = np.load(os.path.join(report_dir, set_name + '_envBAA_intensity_1Txnr_stack.npy'))  # Bedload intensity
    envBAA_act_periods_1Txnr_stack = np.load(os.path.join(report_dir, set_name + '_envBAA_active_periods_1Txnr_stack.npy'))  # Bedload active periods
    envBAA_intensity_comp_stack = np.load(os.path.join(report_dir, set_name + '_envBAA_intensity_comp_stack.npy'))  # Bedload intensity where comp occurs
    envBAA_intensity_thrs_stack = np.load(os.path.join(report_dir, set_name + '_envBAA_intensity_thrs_stack.npy'))  # Bedload intensity where filt process occurs
    envBAA_act_periods_comp_stack = np.load(os.path.join(report_dir, set_name + '_envBAA_periods_comp_stack.npy'))  # Bedload activity periods where comp occurs
    envBAA_act_periods_thrs_stack = np.load(os.path.join(report_dir, set_name + '_envBAA_periods_thrs_stack.npy'))  # Bedload activity periods filtering process occurs

    # CREATE AND APPLY MORPHOLOGICAL AND NON-MORPHOLOGICAL MASK AS STACK ------
    morph_active_mask     = np.where(np.isnan(DoD_maps), np.nan, np.where(DoD_maps != 0, 1, DoD_maps))
    not_morph_active_mask = np.where(np.isnan(DoD_maps), np.nan, np.where(DoD_maps != 0, 0, 1))
    envBAA_morph_intensity_1Txnr_stack = envBAA_intensity_1Txnr_stack * morph_active_mask
    envBAA_not_morph_intensity_1Txnr_stack = envBAA_intensity_1Txnr_stack * not_morph_active_mask

    # CLEAN DATA BY FLATTENING AND REMOVING NAN AND ZERO VALUES ---------------
    cleaned_intensity_data = {
        'envBAA_intensity_1Txnr_stack': clean_data(envBAA_intensity_1Txnr_stack),
        'envBAA_morph_intensity_1Txnr_stack': clean_data(envBAA_morph_intensity_1Txnr_stack),
        'envBAA_not_morph_intensity_1Txnr_stack': clean_data(envBAA_not_morph_intensity_1Txnr_stack),
        'envBAA_intensity_comp_stack': clean_data(envBAA_intensity_comp_stack),
        'envBAA_intensity_thrs_stack': clean_data(envBAA_intensity_thrs_stack)
    }
    cleaned_periods_data = {
        'envBAA_act_periods_1Txnr_stack': clean_data(envBAA_act_periods_1Txnr_stack),
        'envBAA_act_periods_comp_stack': clean_data(envBAA_act_periods_comp_stack),
        'envBAA_act_periods_thrs_stack': clean_data(envBAA_act_periods_thrs_stack)
    }

    # CREATE BOXPLOTS FOR INTENSITY DATA --------------------------------------
    plt.figure(figsize=(12, 6))
    plt.boxplot(
        [cleaned_intensity_data['envBAA_intensity_1Txnr_stack'],
         cleaned_intensity_data['envBAA_morph_intensity_1Txnr_stack'],
         cleaned_intensity_data['envBAA_not_morph_intensity_1Txnr_stack'],
         cleaned_intensity_data['envBAA_intensity_comp_stack'],
         cleaned_intensity_data['envBAA_intensity_thrs_stack']],
        labels=['Total Intensity', 'Morph Intensity', 'Non-Morph Intensity', 'Comp Intensity', 'Thrs Intensity']
    )
    plt.title(f'Bedload Intensity Comparison for {set_name} - 1Txnr discretization')
    plt.ylabel('Bedload Intensity')
    plt.xlabel('Categories')
    plt.grid(True)
    plt.text(0.5, -0.12, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.savefig(os.path.join(plot_dir_out, set_name +
                '_comp_thrs_bedload_activity.pdf'), dpi=800, bbox_inches='tight')
    plt.show()

    # CREATE BOXPLOTS FOR ACTIVE PERIODS DATA ---------------------------------
    plt.figure(figsize=(12, 6))
    plt.boxplot(
        [cleaned_periods_data['envBAA_act_periods_1Txnr_stack'],
         cleaned_periods_data['envBAA_act_periods_comp_stack'],
         cleaned_periods_data['envBAA_act_periods_thrs_stack']],
        labels=['Total Active Periods', 'Comp Active Periods', 'Thrs Active Periods']
    )
    plt.title(f'Bedload Active Periods Comparison for {set_name} - 1Txnr discretization')
    plt.ylabel('Active Periods')
    plt.xlabel('Categories')
    plt.grid(True)
    plt.text(0.5, -0.12, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.savefig(os.path.join(plot_dir_out, set_name +
                '_comp_thrs_bedload_active_period.pdf'), dpi=800, bbox_inches='tight')
    plt.show()


#%%============================================================================
# BEDLOAD INTENSITY AND ACTIVE PERIODS FOR COMPENSATION AND UNDER THRESHOLD EFFECTS
# CONFINEMENT EFFECT
# =============================================================================
import numpy as np
import os
import matplotlib.pyplot as plt

home_dir = os.getcwd()  # Home directory
report_dir = os.path.join(home_dir, 'output')

set_names = ['q05_1', 'q07_1', 'q10_2', 'q15_2', 'q20_2']

# Convert set_names to corresponding numeric values for x-axis
set_name_values = [float(set_name[1:3]) / 10 for set_name in set_names]

def clean_data(arr):
    flat_arr = arr.flatten()
    return flat_arr[(~np.isnan(flat_arr)) & (flat_arr != 0)]

# Initialize dictionaries to store cleaned data for all sets
all_cleaned_intensity_data = {
    'envBAA_intensity_1Txnr_stack': [],
    'envBAA_morph_intensity_1Txnr_stack': [],
    'envBAA_not_morph_intensity_1Txnr_stack': [],
    'envBAA_intensity_comp_stack': [],
    'envBAA_intensity_thrs_stack': []
}

all_cleaned_periods_data = {
    'envBAA_act_periods_1Txnr_stack': [],
    'envBAA_act_periods_comp_stack': [],
    'envBAA_act_periods_thrs_stack': []
}

for set_name in set_names:
    print(set_name + ' is running...')

    # IMPORT DATA -------------------------------------------------------------
    DoD_maps = np.load(os.path.join(report_dir, set_name + '_DoD_1Txnr_stack.npy'))
    envBAA_intensity_1Txnr_stack = np.load(os.path.join(report_dir, set_name + '_envBAA_intensity_1Txnr_stack.npy'))  # Bedload intensity
    envBAA_act_periods_1Txnr_stack = np.load(os.path.join(report_dir, set_name + '_envBAA_active_periods_1Txnr_stack.npy'))  # Bedload active periods
    envBAA_intensity_comp_stack = np.load(os.path.join(report_dir, set_name + '_envBAA_intensity_comp_stack.npy'))  # Bedload intensity where comp occurs
    envBAA_intensity_thrs_stack = np.load(os.path.join(report_dir, set_name + '_envBAA_intensity_thrs_stack.npy'))  # Bedload intensity where filt process occurs
    envBAA_act_periods_comp_stack = np.load(os.path.join(report_dir, set_name + '_envBAA_periods_comp_stack.npy'))  # Bedload activity periods where comp occurs
    envBAA_act_periods_thrs_stack = np.load(os.path.join(report_dir, set_name + '_envBAA_periods_thrs_stack.npy'))  # Bedload activity periods filtering process occurs

    # CREATE  AND APPLY MORPHOLOGICAL AND NON-MORPHOLOGICAL MASK AS STACK -----
    morph_active_mask     = np.where(np.isnan(DoD_maps), np.nan, np.where(DoD_maps != 0, 1, DoD_maps))
    not_morph_active_mask = np.where(np.isnan(DoD_maps), np.nan, np.where(DoD_maps != 0, 0, 1))
    envBAA_morph_intensity_1Txnr_stack = envBAA_intensity_1Txnr_stack * morph_active_mask
    envBAA_not_morph_intensity_1Txnr_stack = envBAA_intensity_1Txnr_stack * not_morph_active_mask

    # CLEAN DATA BY FLATTENING AND REMOVING NAN AND ZERO VALUES ---------------
    all_cleaned_intensity_data['envBAA_intensity_1Txnr_stack'].append(clean_data(envBAA_intensity_1Txnr_stack))
    all_cleaned_intensity_data['envBAA_morph_intensity_1Txnr_stack'].append(clean_data(envBAA_morph_intensity_1Txnr_stack))
    all_cleaned_intensity_data['envBAA_not_morph_intensity_1Txnr_stack'].append(clean_data(envBAA_not_morph_intensity_1Txnr_stack))
    all_cleaned_intensity_data['envBAA_intensity_comp_stack'].append(clean_data(envBAA_intensity_comp_stack))
    all_cleaned_intensity_data['envBAA_intensity_thrs_stack'].append(clean_data(envBAA_intensity_thrs_stack))
    all_cleaned_periods_data['envBAA_act_periods_1Txnr_stack'].append(clean_data(envBAA_act_periods_1Txnr_stack))
    all_cleaned_periods_data['envBAA_act_periods_comp_stack'].append(clean_data(envBAA_act_periods_comp_stack))
    all_cleaned_periods_data['envBAA_act_periods_thrs_stack'].append(clean_data(envBAA_act_periods_thrs_stack))

# CREATE BOXPLOTS FOR INTENSITY DATA ------------------------------------------
for key in all_cleaned_intensity_data:
    plt.figure(figsize=(12, 6))
    data = [all_cleaned_intensity_data[key][i] for i in range(len(set_names))]
    plt.boxplot(data, positions=set_name_values, widths=0.05)
    plt.title(f'Comparison of {key} Data for Different Sets')
    plt.ylabel('Intensity')
    plt.xlabel('Discharge [l/s]')
    plt.xticks(set_name_values, set_name_values)
    plt.grid(True)
    plt.text(0.5, -0.12, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.savefig(os.path.join(plot_dir_out, f'{key}_intensity_comparison.pdf'), dpi=800, bbox_inches='tight')
    plt.show()

# CREATE BOXPLOTS FOR PERIODS DATA --------------------------------------------
for key in all_cleaned_periods_data:
    plt.figure(figsize=(12, 6))
    data = [all_cleaned_periods_data[key][i] for i in range(len(set_names))]
    plt.boxplot(data, positions=set_name_values, widths=0.05)
    plt.title(f'Comparison of {key} Data for Different Sets')
    plt.ylabel('Active Periods')
    plt.xlabel('Discharge [l/s]')
    plt.xticks(set_name_values, set_name_values)
    plt.grid(True)
    plt.text(0.5, -0.12, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.savefig(os.path.join(plot_dir_out, f'{key}_periods_comparison.pdf'), dpi=800, bbox_inches='tight')
    plt.show()


