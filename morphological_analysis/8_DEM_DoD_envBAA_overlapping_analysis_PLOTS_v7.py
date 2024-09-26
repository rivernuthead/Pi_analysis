#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 10:17:19 2024

@author: erri
"""
# =============================================================================
# IMPORT LIBRARIES
# =============================================================================
from scipy import stats
import numpy as np
import os
import matplotlib.pyplot as plt

# =============================================================================
# FUNCTIONS
# =============================================================================
def generate_colors(color_map_name, num_colors):
    '''
    Given a color ramp, this function extract a list of a given number of 
    colors equally separated inside the color ramp

    Parameters
    ----------
    color_map_name : TYPE
        DESCRIPTION.
    num_colors : TYPE
        DESCRIPTION.

    Returns
    -------
    colors : TYPE
        DESCRIPTION.

    '''
    # Get colormap from name
    color_map = plt.get_cmap(color_map_name)

    # Generate equally spaced values from 0 to 1
    values = np.linspace(0, 1, num_colors)

    # Generate colors from colormap
    colors = [color_map(value) for value in values]

    return colors

# =============================================================================
# SETUP PATHS AND FOLDERS
# =============================================================================
home_dir = os.getcwd()  # Home directory
report_dir = os.path.join(home_dir, 'output')
output_folder = os.path.join(report_dir, 'topographic_PiQs_overlapping')
run_dir = os.path.join(home_dir, 'surveys')
plot_dir = os.path.join(os.getcwd(), 'plots', 'DEM_DOD_envBAA_overlapping')
if not(os.path.exists(plot_dir)):
        os.mkdir(plot_dir)

# =============================================================================
# SCRIPT PARAMETERS
# =============================================================================
run_mode = [
    # 'full_BAA_stack',
    'partial_envBAA_stack'
]
'''
run_mode.
        - full_BAA_stack: full BAA map stack
        - partial_envBAA_stack: stack computed as the stack of partial BAA envelopes (see /PiQs_analysis/4_PiQs_BAW_envelopes_generator.py for details)
'''

analysis_mode = ['activity_time_analysis']

run_names = ['q05r1', 'q05r2', 'q05r3', 'q05r4',
             'q05r5', 'q05r6', 'q05r7', 'q05r8', 'q05r9']
set_names = ['q05_1', 'q05_1', 'q05_1', 'q05_1',
             'q05_1', 'q05_1', 'q05_1', 'q05_1', 'q05_1']

run_names = ['q07r1', 'q07r2', 'q07r3', 'q07r4',
              'q07r5', 'q07r6', 'q07r7', 'q07r8', 'q07r9']
set_names = ['q07_1', 'q07_1', 'q07_1', 'q07_1',
              'q07_1', 'q07_1', 'q07_1', 'q07_1', 'q07_1']

run_names = ['q10r1', 'q10r2', 'q10r3', 'q10r4',
              'q10r5', 'q10r6', 'q10r7', 'q10r8', 'q10r9']
set_names = ['q10_2', 'q10_2', 'q10_2', 'q10_2',
              'q10_2', 'q10_2', 'q10_2', 'q10_2', 'q10_2']

run_names = ['q15r1', 'q15r3', 'q15r4', 'q15r5', 'q15r6', 'q15r8', 'q15r9']
set_names = ['q15_2', 'q15_2', 'q15_2', 'q15_2', 'q15_2', 'q15_2', 'q15_2']

run_names = ['q20r1', 'q20r2', 'q20r3', 'q20r4',
              'q20r5', 'q20r6', 'q20r7', 'q20r8', 'q20r9']
set_names = ['q20_2', 'q20_2', 'q20_2', 'q20_2',
              'q20_2', 'q20_2', 'q20_2', 'q20_2', 'q20_2']

# run_names = ['q07r1']
# set_names = ['q07_1']

# run_names = ['q10r1']
# set_names = ['q10_2']

# run_names = ['q15r1']
# set_names = ['q15_2']

# run_names = ['q20r1']
# set_names = ['q20_2']

# =============================================================================
# INITIALIZE ARRAY FOR DISTRIBUTION ANALYSIS
# =============================================================================
morpho_active_array_list = []
bedload_active_array_list = []
travelling_bedload_array_list = []
bedload_active_time_array_list = []
bedload_morph_active_time_array_list = []
bedload_travelling_active_time_array_list = []

DEM_morph_act_hist_list = []
DEM_travelling_bedload_hist_list = []
DEM_never_active_hist_list = []
DEM_travelling_8th_bedload_hist_list = []
DEM_travelling_full_bedload_hist_list = []

travelling_active_area_list = []
total_active_area_list = []

'''
THIS VERY FIRST SECTION IS JUST IN PREPARATION FOR THE PLOT SECTION THAT IS THE
SCRIPT CORE
'''
# =============================================================================
# DATA PREPARATION
# =============================================================================
for run_name, set_name in zip(run_names, set_names):

    plot_dir_out = os.path.join(plot_dir, set_name)
    if not(os.path.exists(plot_dir_out)):
        os.mkdir(plot_dir_out)

    # IMPORT DEM DoD envBAA STACK FROM DEM_DoD_envBAA_overlapping_v3.py
    if 'full_BAA_stack' in run_mode:
        path = os.path.join(output_folder, set_name + '_' +
                            run_name + '_DEM_DoD_envBAA_stack.npy')
    elif 'partial_envBAA_stack' in run_mode:
        path = os.path.join(output_folder, set_name + '_' +
                            run_name + '_DEM_DoD_partial_envBAA_stack.npy')
    DEM_DoD_envBAA_stack = np.load(path)

    # CREATE THE MORPHOLOGICAL ACTIVITY MASK ----------------------------------
    # Morphological activity mask
    MA_mask_nan = np.where(
        DEM_DoD_envBAA_stack[1, :, :] == 0, np.nan, DEM_DoD_envBAA_stack[1, :, :])
    MA_mask_nan = np.where(abs(MA_mask_nan) > 0, 1, MA_mask_nan)
    MA_mask = np.where(np.isnan(MA_mask_nan), 0, MA_mask_nan)
    MA_mask_nan_inverse = np.where(np.isnan(DEM_DoD_envBAA_stack[1, :, :]), np.nan, np.where(
        DEM_DoD_envBAA_stack[1, :, :] == 0, 1, np.nan))

    # CREATE THE BEDLOAD ACTIVITY MASK ----------------------------------------
    # Bedload activity mask
    BA_mask_nan = np.where(
        DEM_DoD_envBAA_stack[3, :, :] == 0, np.nan, DEM_DoD_envBAA_stack[3, :, :])
    BA_mask_nan = np.where(abs(BA_mask_nan) > 0, 1, BA_mask_nan)
    BA_mask = np.where(np.isnan(BA_mask_nan), 0,
                       BA_mask_nan)  # Bedload activity mask
    BA_mask_nan_inverse = np.where(np.isnan(DEM_DoD_envBAA_stack[2, :, :]), np.nan, np.where(
        DEM_DoD_envBAA_stack[2, :, :] == 0, 1, np.nan))

    # ALL DEM VALUES ----------------------------------------------------------
    DEM_values = DEM_DoD_envBAA_stack[0, :, :].flatten()
    DEM_values = DEM_values[~np.isnan(DEM_values)]

    # MORPHOLOGICAL ACTIVE DEM VALUES -----------------------------------------
    DEM_morph_act_matrix = DEM_DoD_envBAA_stack[0, :, :]*MA_mask_nan
    DEM_morph_act = DEM_morph_act_matrix.flatten()
    DEM_morph_act = DEM_morph_act[~np.isnan(DEM_morph_act)]

    # =========================================================================
    # DISTRIBUTION ANLYSIS WITH 7 DEM CLASSES
    # =========================================================================
    # DISTRIBUTION PARAMETERS -------------------------------------------------
    num_bins = 7
    DEM_bin_edges = [-60.0, -8.64, -3.86, -0.93, 1.36, 3.36, 5.44, 60.0]
    DEM_class_labels = [
        f"{DEM_bin_edges[i]},{DEM_bin_edges[i+1]}" for i in range(len(DEM_bin_edges) - 1)]
    hist_range = (-60, 60)  # Range of values to include in the histogram
    dist_bins = DEM_bin_edges
    # COMPUTE HISTOGRAM -------------------------------------------------------
    DEM_morph_act_hist, DEM_morph_bin_edges = np.histogram(
        DEM_morph_act, bins=dist_bins, range=hist_range)
    # DEM_travelling_bedload_hist , DEM_travelling_bedload_bin_edges = np.histogram(DEM_travelling_bedload, bins=dist_bins, range=hist_range)
    # DEM_never_active_hist, DEM_never_active_bin_edges = np.histogram(DEM_never_active, bins=dist_bins, range=hist_range)
    DEM_morph_bin_midpoints = 0.5 * \
        (DEM_morph_bin_edges[1:] + DEM_morph_bin_edges[:-1]) # Calculate midpoints of bins


#%%============================================================================
# RUN PLOT SCRIPT
# =============================================================================
''' 
In this section I want to investigate if it is visible a discharge effect
evaluating the DEM values that:
    1. All the DEM values
    2. All the DEM values that are morphologically active
    3. All the DEM values that are bedload active
    4. All the DEM values that experienced non-morphological bedload
To get the mean and the stdev remember to run the code using the entire run set
not the single run.
'''

set_names = ['q05_1', 'q07_1', 'q10_2', 'q15_2', 'q20_2']
# set_names = ['q07_1']
# set_names = ['q20_2']

color_map_name = 'viridis'  # Specify the name of the colormap
num_colors = len(set_names)  # Specify the number of colors you want
colors = generate_colors(color_map_name, num_colors)


# =============================================================================
# PLOT DEM VALUES DISTRIBUTION OF FREQUENCY
# (AVERAGE AND STDEV) FOR DIFFERENT DISCHARGES
# =============================================================================
for set_name, color in zip(set_names, colors):
    plot_dir_out = os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', set_name)
    DEM_values_hist_data = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_DEM_values_distribution_average.txt'))
    x_values = DEM_values_hist_data[0, :]
    mean_values = DEM_values_hist_data[1, :]
    stdev_values = DEM_values_hist_data[2, :]
    plt.plot(x_values, mean_values, label=set_name, color=color)
    # Plot the shaded area for +/- one standard deviation
    plt.fill_between(x_values, mean_values - stdev_values, mean_values +
                     stdev_values, color=color, alpha=0.3)  # , label='± 1 Stdev')
plt.xlabel('DEM values [mm]')
# plt.ylabel('Y values')
plt.title('All DEM values distribution')
plt.legend()
plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
         transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
plt.savefig(os.path.join(plot_dir_out,
            'all_DEM_values_average_distribution.pdf'))
plt.show()


# =============================================================================
# # PLOT MORPHOLOGICAL ACTIVE DEM VALUES DISTRIBUTION OF FREQUENCY
# (AVERAGE AND STDEV) FOR DIFFERENT DISCHARGES
# =============================================================================
for set_name, color in zip(set_names, colors):
    plot_dir_out = os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', set_name)
    morph_active_values_hist_data = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + 'DEM_morph_active_values_distribution_average.txt'))

    x_values = morph_active_values_hist_data[0, :]
    mean_values = morph_active_values_hist_data[1, :]
    stdev_values = morph_active_values_hist_data[2, :]
    plt.plot(x_values, mean_values, label=set_name, color=color)
    # Plot the shaded area for +/- one standard deviation
    plt.fill_between(x_values, mean_values - stdev_values, mean_values +
                     stdev_values, color=color, alpha=0.3)  # , label='± 1 Stdev')
plt.xlabel('DEM values [mm]')
# plt.ylabel('Y values')
plt.title('Morphologically active DEM values distribution')
plt.legend()
plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
         transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
plt.savefig(os.path.join(plot_dir_out,
            'morph_active_DEM_values_average_distribution.pdf'))
plt.show()


# =============================================================================
# # PLOT BEDLOAD ACTIVE DEM VALUES DISTRIBUTION OF FREQUENCY
# (AVERAGE AND STDEV) FOR DIFFERENT DISCHARGES
# =============================================================================
for set_name, color in zip(set_names, colors):
    plot_dir_out = os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', set_name)
    bedload_active_values_hist_data = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + 'DEM_bedload_active_values_distribution_average.txt'))

    x_values = bedload_active_values_hist_data[0, :]
    mean_values = bedload_active_values_hist_data[1, :]
    stdev_values = bedload_active_values_hist_data[2, :]
    plt.plot(x_values, mean_values, label=set_name, color=color)
    plt.fill_between(x_values, mean_values - stdev_values, mean_values +
                     stdev_values, color=color, alpha=0.3)  # , label='± 1 Stdev')
plt.xlabel('DEM values [mm]')
# plt.ylabel('Y values')
plt.title('Bedload active DEM values distribution')
plt.legend()
plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
         transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
plt.savefig(os.path.join(plot_dir_out,
            'bedload_active_DEM_values_average_distribution.pdf'))
plt.show()


# =============================================================================
# PLOT NON-MORPHOLOCIAL ACTIVE DEM VALUES DISTRIBUTION OF FREQUENCY
# (AVERAGE AND STDEV) FOR DIFFERENT DISCHARGES
# =============================================================================
for set_name, color in zip(set_names, colors):
    plot_dir_out = os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', set_name)
    travelling_bedload_active_values_hist_data = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + 'DEM_travelling_bedload_active_values_distribution_average.txt'))
    x_values = travelling_bedload_active_values_hist_data[0, :]
    mean_values = travelling_bedload_active_values_hist_data[1, :]
    stdev_values = travelling_bedload_active_values_hist_data[2, :]
    plt.plot(x_values, mean_values, label=set_name, color=color)
    # Plot the shaded area for +/- one standard deviation
    plt.fill_between(x_values, mean_values - stdev_values, mean_values +
                     stdev_values, color=color, alpha=0.3)  # , label='± 1 Stdev')
plt.xlabel('DEM values [mm]')
# plt.ylabel('Y values')
plt.title('No morph bedload bedload active DEM values distribution')
plt.legend()
plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
         transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
plt.savefig(os.path.join(plot_dir_out,
            'no_morph_bedload_active_DEM_values_average_distribution.pdf'))
plt.show()


# =============================================================================
# PLOT CHARTS
# =============================================================================
for set_name, color in zip(set_names, colors):
    
    # =========================================================================
    # PLOT DEM VALUES DISTRIBUTION OF FREQUENCY - FOR EACH DISCHARGE
    # MORPHOLOGICAL, NON-MORPHOLOGICAL, NEVER ACTIVE
    # =========================================================================
    plot_dir_out = os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', set_name)
    DEM_morph_act_hist_mean_data = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_DEM_morph_act_hist_mean.txt'))
    DEM_travelling_bedload_hist_mean_data = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_DEM_travelling_bedload_hist_mean.txt'))
    DEM_never_active_hist_mean_data = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_DEM_never_active_hist_mean.txt'))

    x_values = DEM_morph_bin_midpoints
    # All the data have to be shown as stached line:
    DEM_morph_act_hist_mean_values_data = DEM_morph_act_hist_mean_data[0, :]
    DEM_morph_act_hist_stdev_values_data = DEM_morph_act_hist_mean_data[1, :]
    DEM_travelling_bedload_hist_mean_values_data = DEM_morph_act_hist_mean_values_data + \
        DEM_travelling_bedload_hist_mean_data[0, :]
    DEM_travelling_bedload_hist_stdev_values_data = DEM_travelling_bedload_hist_mean_data[
        1, :]
    DEM_never_active_hist_mean_values_data = DEM_travelling_bedload_hist_mean_values_data + \
        DEM_never_active_hist_mean_data[0, :]
    DEM_never_active_hist_stdev_values_data = DEM_never_active_hist_mean_data[1, :]

    # Make the plot dimensionless over the total number of pixel in each class
    # hist_total_sum = np.nansum((DEM_morph_act_hist_mean_values,DEM_travelling_bedload_hist_mean_values,DEM_never_active_hist_mean_values))


    # =========================================================================
    # PLOT RELATIVE MORPHOLOGICAL, TRAVELLING AND NEVER ACTIVE
    # =========================================================================
    # Plot the mean line
    plt.plot(x_values, DEM_morph_act_hist_mean_values_data, label='Morph active')
    plt.plot(x_values, DEM_travelling_bedload_hist_mean_values_data,
             label='"No morph bedload"')
    plt.plot(x_values, DEM_never_active_hist_mean_values_data)
    # Plot the shaded area for +/- one standard deviation
    plt.fill_between(x_values, DEM_morph_act_hist_mean_values_data - DEM_morph_act_hist_stdev_values_data,
                     DEM_morph_act_hist_mean_values_data + DEM_morph_act_hist_stdev_values_data, color=color, alpha=0.3)  # , label='± 1 Stdev')
    plt.fill_between(x_values, DEM_travelling_bedload_hist_mean_values_data - DEM_travelling_bedload_hist_stdev_values_data,
                     DEM_travelling_bedload_hist_mean_values_data + DEM_travelling_bedload_hist_stdev_values_data, color=color, alpha=0.3)  # , label='± 1 Stdev')
    # plt.fill_between(x_values, DEM_never_active_hist_mean_values_data - DEM_never_active_hist_stdev_values_data,
    #                  DEM_never_active_hist_mean_values_data + DEM_never_active_hist_stdev_values_data, color=color, alpha=0.3)  # , label='± 1 Stdev')

    plt.xticks(DEM_morph_bin_midpoints, DEM_class_labels, rotation=60)
    plt.xlabel('DEM values [mm]')
    # plt.ylabel('Y values')
    plt.title('Morphological and no morphological active  - lines - ' + set_name)
    plt.legend()
    plt.ylim(bottom=0, top=1)
    plt.text(0.8, 0.92, "Never active", transform=plt.gca().transAxes,
             ha='center', va='center', fontsize=10)
    plt.text(0.5, -0.4, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.savefig(os.path.join(plot_dir_out, set_name +
                '_morpho_no_morph_never_active_DEM_values_average_distribution.pdf'), bbox_inches='tight')
    plt.show()

    # =========================================================================
    # STACKED BARS CHART
    # =========================================================================
    categories = DEM_class_labels
    values1 = DEM_morph_act_hist_mean_data[0, :]
    values2 = DEM_travelling_bedload_hist_mean_data[0, :]
    values3 = DEM_never_active_hist_mean_data[0, :]

    plt.figure(figsize=(8, 6))
    bars1 = plt.bar(categories, values1, label='Morph active')
    bars2 = plt.bar(categories, values2, bottom=values1,
                    label='no morphological bedload')
    bars3 = plt.bar(categories, values3, bottom=[
                    i+j for i, j in zip(values1, values2)], label='Never active')
    plt.xlabel('DEM values [mm]')  # Rotate x-labels
    plt.ylabel('Percentage')
    plt.title('Morphological and no morphological active  - bars - ' + set_name)
    plt.legend()
    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()
            x_pos = bar.get_x() + bar.get_width() / 2
            y_pos = bar.get_y() + height / 2
            plt.text(x_pos, y_pos, '{:.2f}'.format(height),
                     ha='center', va='center', color='white', fontsize=14)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.text(0.5, -0.4, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.savefig(os.path.join(plot_dir_out, set_name +
                '_morpho_no_morph_never_active_DEM_values_average_distribution_stacked_bars.pdf'), bbox_inches='tight')
    plt.show()


    # =========================================================================
    # PLOT ABSOLUTE MORPHOLOGICAL ACTIVE NAD TRAVELLING BEDLOAD
    # =========================================================================
    plt.plot(x_values, DEM_morph_act_hist_mean_values_data, label='Morph active')
    plt.plot(
        x_values, DEM_travelling_bedload_hist_mean_data[0, :], label='No morph bedload')
    plt.plot(
        x_values, DEM_never_active_hist_mean_data[0, :], label='Never active')

    # Plot the shaded area for +/- one standard deviation
    plt.fill_between(x_values, DEM_morph_act_hist_mean_values_data - DEM_morph_act_hist_stdev_values_data,
                     DEM_morph_act_hist_mean_values_data + DEM_morph_act_hist_stdev_values_data, alpha=0.3)  # , label='± 1 Stdev')
    plt.fill_between(x_values, DEM_travelling_bedload_hist_mean_data[0, :] - DEM_travelling_bedload_hist_stdev_values_data,
                     DEM_travelling_bedload_hist_mean_data[0, :] + DEM_travelling_bedload_hist_stdev_values_data, alpha=0.3)  # , label='± 1 Stdev')
    plt.fill_between(x_values, DEM_never_active_hist_mean_data[0, :] - DEM_never_active_hist_stdev_values_data,
                     DEM_never_active_hist_mean_data[0, :] + DEM_never_active_hist_stdev_values_data, alpha=0.3)  # , label='± 1 Stdev')
    plt.xticks(DEM_morph_bin_midpoints, DEM_class_labels, rotation=60)
    plt.xlabel('DEM values [mm]')
    # plt.ylabel('Y values')
    plt.title(' Morphological and no morphological active absolute  - lines - ' + set_name)
    plt.legend()
    plt.ylim(bottom=0, top=1)
    plt.text(0.5, -0.4, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.savefig(os.path.join(plot_dir_out, set_name +
                '_morpho_no_morph_never_active_DEM_values_average_distribution_absolute.pdf'), bbox_inches='tight')
    plt.show()

    # =========================================================================
    # PLOT STACKED BARS CHART
    # =========================================================================
    # DATA IMPORT -------------------------------------------------------------
    categories = DEM_class_labels
    values1 = DEM_morph_act_hist_mean_values_data
    values2 = DEM_travelling_bedload_hist_mean_data[0, :]
    values3 = DEM_never_active_hist_mean_data[0, :]

    # Calculate the width of each bar
    bar_width = 0.25
    index = np.arange(len(categories))

    # Plot
    plt.figure(figsize=(8, 6))
    bars1 = plt.bar(index - bar_width, values1,
                    bar_width, label='Morph active')
    bars2 = plt.bar(index, values2, bar_width, label='No morphological\nbedload')
    bars3 = plt.bar(index + bar_width, values3,
                    bar_width, label='Never active')
    plt.xlabel('DEM values [mm]')  # Rotate x-labels
    plt.ylabel('Percentage')
    plt.title('Morphological and no morphological active absolute  - bars - ' + set_name)
    plt.legend()
    plt.ylim(0, 0.8)
    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()
            x_pos = bar.get_x() + bar.get_width() / 2
            y_pos = bar.get_y() + height / 2
            plt.text(x_pos, y_pos, '{:.2f}'.format(height),
                     ha='center', va='center', rotation=90, color='black', fontsize=14)
    plt.xticks(index, categories, rotation=45)
    plt.tight_layout()
    plt.text(0.5, -0.4, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.savefig(os.path.join(plot_dir_out, set_name +
                '_morpho_no_morph_never_active_DEM_values_average_distribution_absolute_bars.pdf'), bbox_inches='tight')
    plt.show()

    # =============================================================================
    # PLOT 1/8th AND 1/2nd BEDLOAD ACTIVE DEM VALUES
    # =============================================================================

    # DATA IMPORT -------------------------------------------------------------
    DEM_travelling_8th_bedload_hist_mean_data = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_DEM_travelling_8th_bedload_hist_mean.txt'))
    DEM_travelling_full_bedload_hist_mean_data = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_DEM_travelling_full_bedload_hist_mean.txt'))

    # All the data have to be shown as stached line:
    DEM_travelling_8th_bedload_hist_mean_values_data = DEM_morph_act_hist_mean_values_data + \
        DEM_travelling_8th_bedload_hist_mean_data[0, :]
    DEM_travelling_8th_bedload_hist_stdev_values_data = DEM_travelling_8th_bedload_hist_mean_data[
        1, :]

    DEM_travelling_full_bedload_hist_mean_values_data = DEM_travelling_8th_bedload_hist_mean_values_data + \
        DEM_travelling_full_bedload_hist_mean_data[0, :]
    DEM_travelling_full_bedload_hist_stdev_values_data = DEM_travelling_full_bedload_hist_mean_data[
        1, :]

    # Plot the mean line
    plt.plot(x_values, DEM_morph_act_hist_mean_values_data, label='morph active')
    plt.plot(x_values, DEM_travelling_8th_bedload_hist_mean_values_data,
             label='1/8 Txnr active\nno morph bedload')
    plt.plot(x_values, DEM_travelling_full_bedload_hist_mean_values_data,
             label='full run active\nno morph bedload')
    plt.plot(x_values, DEM_travelling_bedload_hist_mean_values_data,
             label='total no morph bedload')
    plt.plot(x_values, DEM_never_active_hist_mean_values_data)

    # Plot the shaded area for +/- one standard deviation
    plt.fill_between(x_values, DEM_morph_act_hist_mean_values_data - DEM_morph_act_hist_stdev_values_data,
                     DEM_morph_act_hist_mean_values_data + DEM_morph_act_hist_stdev_values_data, alpha=0.3)  # , label='± 1 Stdev')
    plt.fill_between(x_values, DEM_travelling_8th_bedload_hist_mean_values_data - DEM_travelling_8th_bedload_hist_stdev_values_data,
                     DEM_travelling_8th_bedload_hist_mean_values_data + DEM_travelling_8th_bedload_hist_stdev_values_data, alpha=0.3)  # , label='± 1 Stdev')
    plt.fill_between(x_values, DEM_travelling_full_bedload_hist_mean_values_data - DEM_travelling_full_bedload_hist_stdev_values_data,
                     DEM_travelling_full_bedload_hist_mean_values_data + DEM_travelling_full_bedload_hist_stdev_values_data, alpha=0.3)  # , label='± 1 Stdev')
    plt.fill_between(x_values, DEM_travelling_bedload_hist_mean_values_data - DEM_travelling_bedload_hist_stdev_values_data,
                     DEM_travelling_bedload_hist_mean_values_data + DEM_travelling_bedload_hist_stdev_values_data, alpha=0.3)  # , label='± 1 Stdev')
    plt.ylim(bottom=0, top=1)
    plt.xticks(DEM_morph_bin_midpoints, DEM_class_labels, rotation=60)
    plt.xlabel('DEM values [mm]')
    # plt.ylabel('Y values')
    plt.title('Morphological and divided no morph active - lines - ' + set_name)
    plt.legend()
    plt.text(0.8, 0.92, "Never active", transform=plt.gca().transAxes,
             ha='center', va='center', fontsize=10)
    plt.text(0.5, -0.4, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.savefig(os.path.join(plot_dir_out, set_name +
                '_morpho_no_morph_never_active_DEM_values_average_distribution_mod.pdf'), bbox_inches='tight')
    plt.show()


    # =========================================================================
    # STACKED BARS CHART 
    # =========================================================================
    categories = DEM_class_labels
    values1 = DEM_morph_act_hist_mean_data[0, :]
    values2 = DEM_travelling_8th_bedload_hist_mean_data[0, :]
    values3 = DEM_travelling_full_bedload_hist_mean_data[0, :]
    values4 = DEM_travelling_bedload_hist_mean_data[0, :] - values3 - values2
    values5 = DEM_never_active_hist_mean_data[0, :]

    plt.figure(figsize=(8, 6))
    bars1 = plt.bar(categories, values1, label='Morph active')
    bars2 = plt.bar(categories, values2, bottom=values1,
                    label='1/8 no morph bedload')
    bars3 = plt.bar(categories, values3, bottom=[
                    i+j for i, j in zip(values1, values2)], label='full run no morph bedload')
    bars4 = plt.bar(categories, values4, bottom=values1 +
                    values2 + values3, label='Remaining no morph bedload')
    bars5 = plt.bar(categories, values5, bottom=values1 +
                    values2 + values3 + values4, label='Never active')
    plt.xlabel('DEM values [mm]')  # Rotate x-labels
    plt.ylabel('Percentage')
    plt.title(
        'Morphological and divided no morph active - stacked bars - ' + set_name)
    plt.legend()

    # Annotate values on top of each bar
    for bars in [bars1, bars2, bars3, bars4, bars5]:
        for bar in bars:
            height = bar.get_height()
            x_pos = bar.get_x() + bar.get_width() / 2
            y_pos = bar.get_y() + height / 2
            plt.text(x_pos, y_pos, '{:.2f}'.format(height),
                     ha='center', va='center', color='white', fontsize=14)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.text(0.5, -0.4, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.savefig(os.path.join(plot_dir_out, set_name +
                '_morpho_no_morph_never_active_DEM_values_average_distribution_mod_stacked_bars.pdf'), bbox_inches='tight')
    plt.show()


# =============================================================================
# BEDLOAD INTENSITY STACKED BARS CHART
# =============================================================================
    # DATA IMPORT -------------------------------------------------------------
    bedload_intensity_hist_mean_data = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_hist_bedload_intensity_hist_mean.txt'))
    bedload_intensity_8th_hist_mean_data = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_hist_bedload_intensity_8th_hist_mean.txt'))
    bedload_intensity_full_hist_mean_data = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_hist_bedload_intensity_full_hist_mean.txt'))

    bedload_intensity_hist_mean = bedload_intensity_hist_mean_data[0, :]
    bedload_intensity_8th_hist_mean = bedload_intensity_8th_hist_mean_data[0, :]
    bedload_intensity_full_hist_mean = bedload_intensity_full_hist_mean_data[0, :]

    values1 = (bedload_intensity_hist_mean - bedload_intensity_8th_hist_mean -
               bedload_intensity_full_hist_mean)/bedload_intensity_hist_mean
    values2 = bedload_intensity_8th_hist_mean/bedload_intensity_hist_mean
    values3 = bedload_intensity_full_hist_mean/bedload_intensity_hist_mean

    plt.figure(figsize=(8, 6))
    x_data = np.linspace(0, len(values1), num=len(values1))
    bars1 = plt.bar(
        x_data, values3, label='full run bedload intensity', color='green', alpha=0.8)
    bars2 = plt.bar(x_data, values2, bottom=values3,
                    label='"1/8 run bedload intensity', color='blue', alpha=0.8)
    bars3 = plt.bar(x_data, values1, bottom=[
                    i+j for i, j in zip(values2, values3)], label='Total bedload', color='orange', alpha=0.8)
    plt.xlabel('Bedload intensity')  # Rotate x-labels
    plt.ylabel('Percentage')
    plt.title('Bedload intensity - stacked bars - ' + set_name)
    plt.legend()
    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = np.round(bar.get_height(), decimals=4)
            x_pos = np.round(bar.get_x() + bar.get_width() / 2, decimals=4)
            y_pos = np.round(bar.get_y() + height / 2, decimals=4)
            plt.text(x_pos, y_pos, '{:.2f}'.format(height),
                     ha='center', va='center', color='black', fontsize=10, rotation=90)
    # plt.xticks(rotation=45)
    plt.tight_layout()
    plt.xlim(0, 30)
    plt.text(0.5, -0.15, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.savefig(os.path.join(plot_dir_out, set_name +
                '_bedload_intensity_average_distribution_mod_stacked_bars.pdf'), bbox_inches='tight')
    plt.show()


# =============================================================================
# PLOT LINE CHART - 1/8th AND 1/2nd OF THE Txnr
# =============================================================================
    ''' Same plot as above but with:
        - difference between morphological active and 1/8 no morph bedload
        - difference between morphological active and full no morph bedload'''

    # DATA IMPORT -------------------------------------------------------------
    DEM_travelling_8th_bedload_hist_mean_data = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_DEM_travelling_8th_bedload_hist_mean.txt'))
    DEM_travelling_full_bedload_hist_mean_data = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_DEM_travelling_full_bedload_hist_mean.txt'))
    # All the data have to be shown as stacked line:
    DEM_travelling_8th_bedload_hist_mean_values_data = DEM_travelling_8th_bedload_hist_mean_data[
        0, :]
    DEM_travelling_8th_bedload_hist_stdev_values_data = DEM_travelling_8th_bedload_hist_mean_data[
        1, :]
    DEM_travelling_full_bedload_hist_mean_values_data = DEM_travelling_full_bedload_hist_mean_data[
        0, :]
    DEM_travelling_full_bedload_hist_stdev_values_data = DEM_travelling_full_bedload_hist_mean_data[
        1, :]

    plt.plot(x_values, DEM_travelling_8th_bedload_hist_mean_values_data,
             label='1/8 Txnr active\nno morph bedload')
    plt.plot(x_values, DEM_travelling_full_bedload_hist_mean_values_data,
             label='full run active\nno morph bedload')
    # Plot the shaded area for +/- one standard deviation
    plt.fill_between(x_values, DEM_travelling_8th_bedload_hist_mean_values_data - DEM_travelling_8th_bedload_hist_stdev_values_data,
                     DEM_travelling_8th_bedload_hist_mean_values_data + DEM_travelling_8th_bedload_hist_stdev_values_data, alpha=0.3)  # , label='± 1 Stdev')
    plt.fill_between(x_values, DEM_travelling_full_bedload_hist_mean_values_data - DEM_travelling_full_bedload_hist_stdev_values_data,
                     DEM_travelling_full_bedload_hist_mean_values_data + DEM_travelling_full_bedload_hist_stdev_values_data, alpha=0.3)  # , label='± 1 Stdev')
    plt.xticks(DEM_morph_bin_midpoints, DEM_class_labels, rotation=60)
    plt.ylim(0, 0.4)
    plt.xlabel('DEM values [mm]')
    # plt.ylabel('Y values')
    plt.title('Divided no morph bedload active, absolute - lines - ' + set_name)
    plt.legend()
    plt.text(0.5, -0.4, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.savefig(os.path.join(plot_dir_out, set_name +
                '_8th_and_full_activity_distribution.pdf'), bbox_inches='tight')
    plt.show()

    # =========================================================================
    # PLOT BARS CHART
    # =========================================================================
    categories = DEM_class_labels
    values1 = DEM_travelling_8th_bedload_hist_mean_values_data
    values2 = DEM_travelling_full_bedload_hist_mean_values_data
    bar_width = 0.25
    index = np.arange(len(categories))
    plt.figure(figsize=(8, 6))
    bars1 = plt.bar(index - 0.5*bar_width, values1,
                    bar_width, label='1/8 run no morph')
    bars2 = plt.bar(index + 0.5*bar_width, values2,
                    bar_width, label='Full run no morph')
    plt.ylim(0, 0.4)
    plt.xlabel('DEM values [mm]')  # Rotate x-labels
    plt.ylabel('Percentage')
    plt.title('Divided no morph bedload active, absolute - bars - ' + set_name)
    plt.legend()
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            x_pos = bar.get_x() + bar.get_width() / 2
            y_pos = bar.get_y() + height / 2
            plt.text(x_pos, y_pos, '{:.2f}'.format(height),
                     ha='center', va='center', rotation=90, color='black', fontsize=14)
    plt.xticks(index, categories, rotation=45)
    plt.tight_layout()
    plt.text(0.5, -0.4, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.savefig(os.path.join(plot_dir_out, set_name +
                '_8th_and_full_activity_distribution_stacked_bars.pdf'), bbox_inches='tight')
    plt.show()

    # =========================================================================
    # BEDLOAD INTENSITY BOXPLOT FOR EACH DISCHARGE
    # =========================================================================
    # LOAD DATA ---------------------------------------------------------------
    values_bedload_intensity_list = np.load(os.path.join(
        plot_dir_out, set_name + '_values_bedload_intensity_list.npy'))
    values_bedload_intensity_8th_list = np.load(os.path.join(
        plot_dir_out, set_name + '_values_bedload_intensity_8th_list.npy'))
    values_bedload_intensity_full_list = np.load(os.path.join(
        plot_dir_out, set_name + '_values_bedload_intensity_full_list.npy'))
    values_bedload_intensity_morpho = np.load(os.path.join(
        plot_dir_out, set_name + '_values_bedload_intensity_morpho_list.npy'))

    # Combine the arrays into a single list of lists
    data = [values_bedload_intensity_list,
            values_bedload_intensity_8th_list,
            values_bedload_intensity_full_list,
            values_bedload_intensity_morpho]
    # =========================================================================
    # PLOT BOXPLOT
    # =========================================================================
    box_colors = ['orange', 'blue', 'green', 'gray']
    box_alpha = 0.5
    boxplot = plt.boxplot(data, patch_artist=True)
    for box, color in zip(boxplot['boxes'], box_colors):
        box.set(color=color, alpha=box_alpha)
    for flier in boxplot['fliers']:
        flier.set(markersize=4)  # Adjust the markersize as needed
    plt.xticks([1, 2, 3, 4], ['bedload\nintensity', 'no morph\nintensity 1/8th',
               'no morph\nintensity full', 'bedload morph\nintensity'])
    plt.title('Boxplot of Bedload Intensity Values - ' + set_name)
    # plt.xlabel('Array')
    plt.ylabel('Bedload intensity')
    plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.savefig(os.path.join(plot_dir_out, set_name +
                '_bedload_travelling_boxplot.pdf'), dpi=800, bbox_inches='tight')
    plt.show()

# =============================================================================
# TEST ANOVA
# =============================================================================
    '''
    Given the large number of points in the dataset, the ANOVA test fails to
    correctly assess differences between distributions. It therefore considers
    the distributions all statistically similar.
    To overcome this problem, a new distribution, consisting of 100 percentiles
    of the original value distribution, was used instead of the original data.
    '''
    import numpy as np
    import pandas as pd
    from scipy import stats
    import statsmodels.api as sm
    from statsmodels.formula.api import ols
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    
    # LOAD DATA
    group1 = np.array([values_bedload_intensity_list]).flatten()
    group2 = np.array([values_bedload_intensity_8th_list]).flatten()
    group3 = np.array([values_bedload_intensity_full_list]).flatten()
    group4 = np.array([values_bedload_intensity_morpho]).flatten()

    percentiles = np.arange(1, 101)  # Percentiles from 1 to 100
    group1 = np.percentile(group1, percentiles)
    group2 = np.percentile(group2, percentiles)
    group3 = np.percentile(group3, percentiles)
    group4 = np.percentile(group4, percentiles)
    
    
    # CREATE DATAFRAME --------------------------------------------------------
    data = {
        'value': np.concatenate([group1, group2, group3, group4]),
        'group': (['total bedload'] * len(group1) + 
                  ['not morph 1/8th'] * len(group2) + 
                  ['not morph full'] * len(group3) + 
                  ['morph bedload'] * len(group4))
    }
    df = pd.DataFrame(data)
    
    # Perform ANOVA
    model = ols('value ~ C(group)', data=df).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    # print(anova_table)
    # Perform Tukey's HSD test
    tukey = pairwise_tukeyhsd(endog=df['value'], groups=df['group'], alpha=0.05)
    # print(tukey)

    # Save the results to a .txt file
    path_report_txt = os.path.join(output_folder, set_name + '_anova_tukey_results_bedlaod.txt')
    with open(path_report_txt, 'w') as file:
        file.write(set_name + " - Morphological bedload intensity values:\n")
        file.write(f"Generated by: {os.path.basename(__file__)}\n")
        file.write("ANOVA Table:\n")
        file.write(anova_table.to_string())
        file.write("\n\nTukey HSD Test:\n")
        file.write(tukey.summary().as_text())
    print("_anova_tukey_results_bedlaod.txt'")

# =============================================================================
# PLOT
# =============================================================================
# PLOT PERCENTAGE OF BEDLOAD ACTIVE TIMES
for set_name in set_names:
    plot_dir_out = os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', set_name)
    DEM_morph_act_hist_mean_data = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_DEM_morph_act_hist_mean.txt'))
    x_values = DEM_morph_bin_midpoints
    # All the data have to be shown as stacked line:
    DEM_morph_act_hist_mean_values_data = DEM_morph_act_hist_mean_data[0, :]
    DEM_morph_act_hist_stdev_values_data = DEM_morph_act_hist_mean_data[1, :]
    plt.plot(x_values, DEM_morph_act_hist_mean_values_data, label=set_name)
    plt.fill_between(x_values, DEM_morph_act_hist_mean_values_data - DEM_morph_act_hist_stdev_values_data,
                     DEM_morph_act_hist_mean_values_data + DEM_morph_act_hist_stdev_values_data, alpha=0.3)  # , label='± 1 Stdev')
plt.xticks(DEM_morph_bin_midpoints, DEM_class_labels, rotation=60)
plt.xlabel('DEM values [mm]')
# plt.ylabel('Y values')
plt.title('Morphological active DEM values')
plt.legend()
plt.text(0.5, -0.4, f"Generated by: {os.path.basename(__file__)}",
         transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
plt.savefig(os.path.join(plot_dir,
            'REPORT_morphological_active_DEM_values_average_distribution.pdf'), bbox_inches='tight')
plt.show()



#%%============================================================================
# BEDLOAD WITHOUT MORPHOLOGICLA CHANGES
# =============================================================================
# LOAD DATA -------------------------------------------------------------------
q05_data_bedload_raw = np.loadtxt(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q05_1',  'q05_1' + '_bedload_active_times_values_distribution_average.txt'))
q05_data_bedload_mean = q05_data_bedload_raw[1, :]
q05_data_bedload_mean = q05_data_bedload_mean[q05_data_bedload_mean != 0]
q05_data_bedload_stdev = q05_data_bedload_raw[2, :]
q05_data_bedload_stdev = q05_data_bedload_stdev[q05_data_bedload_stdev != 0]

q07_data_bedload_raw = np.loadtxt(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q07_1',  'q07_1' + '_bedload_active_times_values_distribution_average.txt'))
q07_data_bedload_mean = q07_data_bedload_raw[1, :]
q07_data_bedload_mean = q07_data_bedload_mean[q07_data_bedload_mean != 0]
q07_data_bedload_stdev = q07_data_bedload_raw[2, :]
q07_data_bedload_stdev = q07_data_bedload_stdev[q07_data_bedload_stdev != 0]

q10_data_bedload_raw = np.loadtxt(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q10_2',  'q10_2' + '_bedload_active_times_values_distribution_average.txt'))
q10_data_bedload_mean = q10_data_bedload_raw[1, :]
q10_data_bedload_mean = q10_data_bedload_mean[q10_data_bedload_mean != 0]
q10_data_bedload_stdev = q10_data_bedload_raw[2, :]
q10_data_bedload_stdev = q10_data_bedload_stdev[q10_data_bedload_stdev != 0]

q15_data_bedload_raw = np.loadtxt(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q15_2',  'q15_2' + '_bedload_active_times_values_distribution_average.txt'))
q15_data_bedload_mean = q15_data_bedload_raw[1, :]
q15_data_bedload_mean = q15_data_bedload_mean[q15_data_bedload_mean != 0]
q15_data_bedload_stdev = q15_data_bedload_raw[2, :]
q15_data_bedload_stdev = q15_data_bedload_stdev[q15_data_bedload_stdev != 0]

q20_data_bedload_raw = np.loadtxt(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q20_2',  'q20_2' + '_bedload_active_times_values_distribution_average.txt'))
q20_data_bedload_mean = q20_data_bedload_raw[1, :]
q20_data_bedload_mean = q20_data_bedload_mean[q20_data_bedload_mean != 0]
q20_data_bedload_stdev = q20_data_bedload_raw[2, :]
q20_data_bedload_stdev = q15_data_bedload_stdev[q15_data_bedload_stdev != 0]

# PLOT BARS CHART -------------------------------------------------------------
plt.figure(figsize=(10, 6))
categories = ['1/8', '1/4', '3/8', '1/2']
plt.ylim(bottom=0, top=1)
bar_width = 0.15
x = np.arange(len(categories))
# Plotting the bars for each dataset
plt.bar(x - 2*bar_width, q05_data_bedload_mean,
        yerr=q05_data_bedload_stdev, width=bar_width, label='q05_1')
plt.bar(x - 1*bar_width, q07_data_bedload_mean,
        yerr=q07_data_bedload_stdev, width=bar_width, label='q07_1')
plt.bar(x, q10_data_bedload_mean,
        yerr=q10_data_bedload_stdev, width=bar_width, label='q10_2')
plt.bar(x + 1*bar_width, q15_data_bedload_mean,
        yerr=q15_data_bedload_stdev, width=bar_width, label='q15_2')
plt.bar(x + 2*bar_width, q20_data_bedload_mean,
        yerr=q20_data_bedload_stdev, width=bar_width, label='q20_2')
plt.xlabel('Exner time')
plt.ylabel('Values')
plt.title('Bedload activity')
plt.xticks(x, categories)
plt.legend()
plt.text(0.5, -0.15, f"Generated by: {os.path.basename(__file__)}",
         transform=plt.gca().transAxes, ha='center', va='center', fontsize=10)
plt.savefig(os.path.join(plot_dir,
            'bedload_activity_percentage_values_average_distribution_bars.pdf'))
plt.show()

#%%============================================================================
# BEDLOAD WITHOUT MORPHOLOGICLA CHANGES
# =============================================================================
# LOAD DATA -------------------------------------------------------------------
q05_data_bedload_morph_raw = np.loadtxt(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q05_1',  'q05_1' + '_bedload_morph_active_times_values_distribution_average.txt'))
q05_data_bedload_morph_mean = q05_data_bedload_morph_raw[1, :]
q05_data_bedload_morph_mean = q05_data_bedload_morph_mean[q05_data_bedload_morph_mean != 0]
q05_data_bedload_morph_stdev = q05_data_bedload_morph_raw[2, :]
q05_data_bedload_morph_stdev = q05_data_bedload_morph_stdev[q05_data_bedload_morph_stdev != 0]

q07_data_bedload_morph_raw = np.loadtxt(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q07_1',  'q07_1' + '_bedload_morph_active_times_values_distribution_average.txt'))
q07_data_bedload_morph_mean = q07_data_bedload_morph_raw[1, :]
q07_data_bedload_morph_mean = q07_data_bedload_morph_mean[q07_data_bedload_morph_mean != 0]
q07_data_bedload_morph_stdev = q07_data_bedload_morph_raw[2, :]
q07_data_bedload_morph_stdev = q07_data_bedload_morph_stdev[q07_data_bedload_morph_stdev != 0]

q10_data_bedload_morph_raw = np.loadtxt(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q10_2',  'q10_2' + '_bedload_morph_active_times_values_distribution_average.txt'))
q10_data_bedload_morph_mean = q10_data_bedload_morph_raw[1, :]
q10_data_bedload_morph_mean = q10_data_bedload_morph_mean[q10_data_bedload_morph_mean != 0]
q10_data_bedload_morph_stdev = q10_data_bedload_morph_raw[2, :]
q10_data_bedload_morph_stdev = q10_data_bedload_morph_stdev[q10_data_bedload_morph_stdev != 0]

q15_data_bedload_morph_raw = np.loadtxt(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q15_2',  'q15_2' + '_bedload_morph_active_times_values_distribution_average.txt'))
q15_data_bedload_morph_mean = q15_data_bedload_morph_raw[1, :]
q15_data_bedload_morph_mean = q15_data_bedload_morph_mean[q15_data_bedload_morph_mean != 0]
q15_data_bedload_morph_stdev = q15_data_bedload_morph_raw[2, :]
q15_data_bedload_morph_stdev = q15_data_bedload_morph_stdev[q15_data_bedload_morph_stdev != 0]

q20_data_bedload_morph_raw = np.loadtxt(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q20_2',  'q20_2' + '_bedload_morph_active_times_values_distribution_average.txt'))
q20_data_bedload_morph_mean = q20_data_bedload_morph_raw[1, :]
q20_data_bedload_morph_mean = q20_data_bedload_morph_mean[q20_data_bedload_morph_mean != 0]
q20_data_bedload_morph_stdev = q20_data_bedload_morph_raw[2, :]
q20_data_bedload_morph_stdev = q20_data_bedload_morph_stdev[q20_data_bedload_morph_stdev != 0]

# PLOT BARS CHART -------------------------------------------------------------
plt.figure(figsize=(12, 6))
plt.ylim(bottom=0, top=1)
bar_width = 0.15
x = np.arange(len(categories))
plt.bar(x - 2*bar_width, q05_data_bedload_morph_mean,
        yerr=q05_data_bedload_morph_stdev, width=bar_width, label='q05_1')
plt.bar(x - 1*bar_width, q07_data_bedload_morph_mean,
        yerr=q07_data_bedload_morph_stdev, width=bar_width, label='q07_1')
plt.bar(x - 0*bar_width, q10_data_bedload_morph_mean,
        yerr=q10_data_bedload_morph_stdev, width=bar_width, label='q10_2')
plt.bar(x + 1*bar_width, q15_data_bedload_morph_mean,
        yerr=q15_data_bedload_morph_stdev, width=bar_width, label='q15_2')
plt.bar(x + 2*bar_width, q20_data_bedload_morph_mean,
        yerr=q20_data_bedload_morph_stdev, width=bar_width, label='q20_2')
plt.xlabel('Exner time')
plt.ylabel('Values')
plt.title('Bedload activity with morphological changes')
plt.xticks(x, categories)
plt.legend()
plt.text(0.5, -0.15, f"Generated by: {os.path.basename(__file__)}",
         transform=plt.gca().transAxes, ha='center', va='center', fontsize=10)
plt.savefig(os.path.join(plot_dir,
            'bedload_morph_activity_percentage_values_average_distribution_bars.pdf'))
plt.show()

#%%============================================================================
# BEDLOAD WITHOUT MORPHOLOGICLA CHANGES
# =============================================================================
# LOAD DATA -------------------------------------------------------------------
q05_data_travel_raw = np.loadtxt(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q05_1',  'q05_1' + '_bedload_travelling_active_times_values_distribution_average.txt'))
q05_data_travel_mean = q05_data_travel_raw[1, :]
q05_data_travel_mean = q05_data_travel_mean[q05_data_travel_mean != 0]
q05_data_travel_stdev = q05_data_travel_raw[2, :]
q05_data_travel_stdev = q05_data_travel_stdev[q05_data_travel_stdev != 0]

q07_data_travel_raw = np.loadtxt(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q07_1',  'q07_1' + '_bedload_travelling_active_times_values_distribution_average.txt'))
q07_data_travel_mean = q07_data_travel_raw[1, :]
q07_data_travel_mean = q07_data_travel_mean[q07_data_travel_mean != 0]
q07_data_travel_stdev = q07_data_travel_raw[2, :]
q07_data_travel_stdev = q07_data_travel_stdev[q07_data_travel_stdev != 0]

q10_data_travel_raw = np.loadtxt(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q10_2',  'q10_2' + '_bedload_travelling_active_times_values_distribution_average.txt'))
q10_data_travel_mean = q10_data_travel_raw[1, :]
q10_data_travel_mean = q10_data_travel_mean[q10_data_travel_mean != 0]
q10_data_travel_stdev = q10_data_travel_raw[2, :]
q10_data_travel_stdev = q10_data_travel_stdev[q10_data_travel_stdev != 0]

q15_data_travel_raw = np.loadtxt(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q15_2',  'q15_2' + '_bedload_travelling_active_times_values_distribution_average.txt'))
q15_data_travel_mean = q15_data_travel_raw[1, :]
q15_data_travel_mean = q15_data_travel_mean[q15_data_travel_mean != 0]
q15_data_travel_stdev = q15_data_travel_raw[2, :]
q15_data_travel_stdev = q15_data_travel_stdev[q15_data_travel_stdev != 0]

q20_data_travel_raw = np.loadtxt(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q20_2',  'q20_2' + '_bedload_travelling_active_times_values_distribution_average.txt'))
q20_data_travel_mean = q20_data_travel_raw[1, :]
q20_data_travel_mean = q20_data_travel_mean[q20_data_travel_mean != 0]
q20_data_travel_stdev = q20_data_travel_raw[2, :]
q20_data_travel_stdev = q20_data_travel_stdev[q20_data_travel_stdev != 0]

# PLOT BARS CHART -------------------------------------------------------------
plt.figure(figsize=(12, 6))
plt.ylim(bottom=0, top=1)
bar_width = 0.15
x = np.arange(len(categories))
plt.bar(x - 2*bar_width, q05_data_travel_mean,
        yerr=q05_data_travel_stdev, width=bar_width, label='q05_1')
plt.bar(x - 1*bar_width, q07_data_travel_mean,
        yerr=q07_data_travel_stdev, width=bar_width, label='q07_1')
plt.bar(x , q10_data_travel_mean,
        yerr=q10_data_travel_stdev, width=bar_width, label='q10_2')
plt.bar(x + 1*bar_width, q15_data_travel_mean,
        yerr=q15_data_travel_stdev, width=bar_width, label='q15_2')
plt.bar(x + 2*bar_width, q20_data_travel_mean,
        yerr=q20_data_travel_stdev, width=bar_width, label='q20_2')
plt.xlabel('Exner time')
plt.ylabel('Values')
plt.title('Bedload activity without morphological changes')
plt.xticks(x, categories)
plt.legend()
plt.text(0.5, -0.15, f"Generated by: {os.path.basename(__file__)}",
         transform=plt.gca().transAxes, ha='center', va='center', fontsize=10)
plt.savefig(os.path.join(plot_dir,
            'Bedload_travelling_activity_percentage_values_average_distribution_bars.pdf'), dpi=1600)
plt.show()


#%%============================================================================
# BOXPLOT TOTAL BEDLOAD
# =============================================================================

# LOAD DATA -------------------------------------------------------------------
q05_1_values_bedload_intensity_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q05_1',  'q05_1' + '_values_bedload_intensity_list.npy'))
q07_1_values_bedload_intensity_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q07_1',  'q07_1' + '_values_bedload_intensity_list.npy'))
q10_2_values_bedload_intensity_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q10_2',  'q10_2' + '_values_bedload_intensity_list.npy'))
q15_2_values_bedload_intensity_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q15_2',  'q15_2' + '_values_bedload_intensity_list.npy'))
q20_2_values_bedload_intensity_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q20_2',  'q20_2' + '_values_bedload_intensity_list.npy'))

data = [q05_1_values_bedload_intensity_list,
        q07_1_values_bedload_intensity_list,
        q10_2_values_bedload_intensity_list,
        q15_2_values_bedload_intensity_list,
        q20_2_values_bedload_intensity_list]

# PLOT BOXPLOT ----------------------------------------------------------------
box_colors = ['#648FFF', '#785EF0', '#DC267F', '#FE6100', '#FFB000']
box_alpha = 0.5
plt.figure(figsize=(12, 6))
boxplot = plt.boxplot(data, patch_artist=True)
# Set box colors and alpha individually
for box, color in zip(boxplot['boxes'], box_colors):
    box.set(color=color, alpha=box_alpha)
# Set outlier marker size
for flier in boxplot['fliers']:
    flier.set(markersize=4)  # Adjust the markersize as needed
plt.xticks([1, 2, 3, 4, 5], ['Q=0.5 l/s', 'Q=0.7 l/s', 'Q=1.0 l/s', 'Q=1.5 l/s', 'Q=2.0 l/s'])
plt.ylim(bottom=0, top=40)
plt.title('Total bedload intensity')
# plt.xlabel('Array')
plt.ylabel('Bedload intensity')
plt.text(0.5, -0.1, f"Generated by: {os.path.basename(__file__)}",
         transform=plt.gca().transAxes, ha='center', va='center', fontsize=8)
plt.savefig(os.path.join(plot_dir, 'REPORT_total_bedload_boxplot.pdf'),
            dpi=800, bbox_inches='tight')
plt.show()

#%%============================================================================
# BOXPLOT 1/8TH DURATION NO MORPHOLOGICAL BEDLOAD
# =============================================================================
# LOAD DATA -------------------------------------------------------------------
q05_1_values_bedload_intensity_8th_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q05_1',  'q05_1' + '_values_bedload_intensity_8th_list.npy'))
q07_1_values_bedload_intensity_8th_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q07_1',  'q07_1' + '_values_bedload_intensity_8th_list.npy'))
q10_2_values_bedload_intensity_8th_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q10_2',  'q10_2' + '_values_bedload_intensity_8th_list.npy'))
q15_2_values_bedload_intensity_8th_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q15_2',  'q15_2' + '_values_bedload_intensity_8th_list.npy'))
q20_2_values_bedload_intensity_8th_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q20_2',  'q20_2' + '_values_bedload_intensity_8th_list.npy'))
data = [q05_1_values_bedload_intensity_8th_list,
        q07_1_values_bedload_intensity_8th_list,
        q10_2_values_bedload_intensity_8th_list,
        q15_2_values_bedload_intensity_8th_list,
        q20_2_values_bedload_intensity_8th_list]

# PLOT BOXPLOT ----------------------------------------------------------------
box_colors = ['#648FFF', '#785EF0', '#DC267F', '#FE6100', '#FFB000']
box_alpha = 0.5
plt.figure(figsize=(12, 6))
boxplot = plt.boxplot(data, patch_artist=True)
# Set box colors and alpha individually
for box, color in zip(boxplot['boxes'], box_colors):
    box.set(color=color, alpha=box_alpha)
# Set outlier marker size
for flier in boxplot['fliers']:
    flier.set(markersize=4)  # Adjust the markersize as needed
plt.xticks([1, 2, 3, 4, 5], ['Q=0.5 l/s', 'Q=0.7 l/s', 'Q=1.0 l/s', 'Q=1.5 l/s', 'Q=2.0 l/s'])
plt.ylim(bottom=0, top=40)
plt.title('1/8th duration non-morphological bedload intensity')
# plt.xlabel('Array')
plt.ylabel('Bedload intensity')
plt.text(0.5, -0.1, f"Generated by: {os.path.basename(__file__)}",
         transform=plt.gca().transAxes, ha='center', va='center', fontsize=8)
plt.savefig(os.path.join(plot_dir,
            'REPORT_8th_no_morph_bedload_boxplot.pdf'), dpi=800, bbox_inches='tight')
plt.show()

#%%============================================================================
# BOXPLOT FULL DURATION NO MORPHOLOGICAL BEDLOAD
# =============================================================================
# LOAD DATA -------------------------------------------------------------------
q05_1_values_bedload_intensity_full_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q05_1',  'q05_1' + '_values_bedload_intensity_full_list.npy'))
q07_1_values_bedload_intensity_full_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q07_1',  'q07_1' + '_values_bedload_intensity_full_list.npy'))
q10_2_values_bedload_intensity_full_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q10_2',  'q10_2' + '_values_bedload_intensity_full_list.npy'))
q15_2_values_bedload_intensity_full_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q15_2',  'q15_2' + '_values_bedload_intensity_full_list.npy'))
q20_2_values_bedload_intensity_full_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q20_2',  'q20_2' + '_values_bedload_intensity_full_list.npy'))
data = [q05_1_values_bedload_intensity_full_list,
        q07_1_values_bedload_intensity_full_list,
        q10_2_values_bedload_intensity_full_list,
        q15_2_values_bedload_intensity_full_list,
        q20_2_values_bedload_intensity_full_list]

# PLOT BOXPLOT ----------------------------------------------------------------
box_colors = ['#648FFF', '#785EF0', '#DC267F', '#FE6100', '#FFB000']
box_alpha = 0.5
plt.figure(figsize=(12, 6))
boxplot = plt.boxplot(data, patch_artist=True)
# Set box colors and alpha individually
for box, color in zip(boxplot['boxes'], box_colors):
    box.set(color=color, alpha=box_alpha)
# Set outlier marker size
for flier in boxplot['fliers']:
    flier.set(markersize=4)  # Adjust the markersize as needed
plt.xticks([1, 2, 3, 4, 5], ['Q=0.5 l/s', 'Q=0.7 l/s', 'Q=1.0 l/s', 'Q=1.5 l/s', 'Q=2.0 l/s'])
plt.ylim(bottom=0, top=40)
plt.title('Full duration non-morphological bedload intensity')
# plt.xlabel('Array')
plt.ylabel('Bedload intensity')
plt.text(0.5, -0.1, f"Generated by: {os.path.basename(__file__)}",
         transform=plt.gca().transAxes, ha='center', va='center', fontsize=8)
plt.savefig(os.path.join(plot_dir,
            'REPORT_morphological_bedload_boxplot.pdf'), dpi=800, bbox_inches='tight')
plt.show()

#%%============================================================================
# BOXPLOT MORPHOLOGICAL ONLY BEDLOAD
# =============================================================================
# LOAD DATA -------------------------------------------------------------------
q05_1_values_bedload_intensity_morpho_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q05_1',  'q05_1' + '_values_bedload_intensity_morpho_list.npy'))
q07_1_values_bedload_intensity_morpho_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q07_1',  'q07_1' + '_values_bedload_intensity_morpho_list.npy'))
q10_2_values_bedload_intensity_morpho_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q10_2',  'q10_2' + '_values_bedload_intensity_morpho_list.npy'))
q15_2_values_bedload_intensity_morpho_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q15_2',  'q15_2' + '_values_bedload_intensity_morpho_list.npy'))
q20_2_values_bedload_intensity_morpho_list = np.load(os.path.join(
    os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q20_2',  'q20_2' + '_values_bedload_intensity_morpho_list.npy'))

data = [q05_1_values_bedload_intensity_morpho_list,
        q07_1_values_bedload_intensity_morpho_list,
        q10_2_values_bedload_intensity_morpho_list,
        q15_2_values_bedload_intensity_morpho_list,
        q20_2_values_bedload_intensity_morpho_list]

# PLOT BOXPLOT ----------------------------------------------------------------
box_colors = ['#648FFF', '#785EF0', '#DC267F', '#FE6100', '#FFB000']
box_alpha = 0.5
plt.figure(figsize=(12, 6))
boxplot = plt.boxplot(data, patch_artist=True)
for box, color in zip(boxplot['boxes'], box_colors):
    box.set(color=color, alpha=box_alpha)
for flier in boxplot['fliers']: # Set outlier marker size
    flier.set(markersize=4)  # Adjust the markersize as needed
plt.xticks([1, 2, 3, 4, 5], ['Q=0.5 l/s', 'Q=0.7 l/s', 'Q=1.0 l/s', 'Q=1.5 l/s', 'Q=2.0 l/s'])
plt.ylim(bottom=0, top=40)
plt.title('Morphological bedload intensity')
# plt.xlabel('Array')
plt.ylabel('Bedload intensity')
plt.text(0.5, -0.1, f"Generated by: {os.path.basename(__file__)}",
         transform=plt.gca().transAxes, ha='center', va='center', fontsize=8)
plt.savefig(os.path.join(plot_dir,
            'REPORT_morphological_bedload_boxplot.pdf'), dpi=800, bbox_inches='tight')
plt.show()

#%%============================================================================
# ANOVA STATISTICAL ANALYSIS - FULL RUN DURATION NO MORPHOLOGICAL BEDLOAD INTENSITY
# =============================================================================
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# LOAD DATA
q05_1_full_no_morph_bedload = np.load(os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q05_1',  'q05_1' + '_values_bedload_intensity_full_list.npy')).flatten()
q07_1_full_no_morph_bedload = np.load(os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q07_1',  'q07_1' + '_values_bedload_intensity_full_list.npy')).flatten()
q10_2_full_no_morph_bedload = np.load(os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q10_2',  'q10_2' + '_values_bedload_intensity_full_list.npy')).flatten()
q15_2_full_no_morph_bedload = np.load(os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q15_2',  'q15_2' + '_values_bedload_intensity_full_list.npy')).flatten()
q20_2_full_no_morph_bedload = np.load(os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q20_2',  'q20_2' + '_values_bedload_intensity_full_list.npy')).flatten()

# TODO TEST use percentiles
percentiles = np.arange(1, 101)  # Percentiles from 1 to 100

q05_1_full_no_morph_bedload = np.percentile(q05_1_full_no_morph_bedload, percentiles)
q07_1_full_no_morph_bedload = np.percentile(q07_1_full_no_morph_bedload, percentiles)
q10_2_full_no_morph_bedload = np.percentile(q10_2_full_no_morph_bedload, percentiles)
q15_2_full_no_morph_bedload = np.percentile(q15_2_full_no_morph_bedload, percentiles)
q20_2_full_no_morph_bedload = np.percentile(q20_2_full_no_morph_bedload, percentiles)

# Prepare the data
data_full_no_morph_bedload = {
    'value': np.concatenate([q05_1_full_no_morph_bedload, q07_1_full_no_morph_bedload, q10_2_full_no_morph_bedload, q15_2_full_no_morph_bedload, q20_2_full_no_morph_bedload]),
    'group': (['Q=0.5 l/s'] * len(q05_1_full_no_morph_bedload) +
              ['Q=0.7 l/s'] * len(q07_1_full_no_morph_bedload) + 
              ['Q=1.0 l/s'] * len(q10_2_full_no_morph_bedload) + 
              ['Q=1.5 l/s'] * len(q15_2_full_no_morph_bedload) + 
              ['Q=2.0 l/s'] * len(q20_2_full_no_morph_bedload))
}
df_full_no_morph_bedload = pd.DataFrame(data_full_no_morph_bedload)

# Perform ANOVA
model = ols('value ~ C(group)', data=df_full_no_morph_bedload).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print(anova_table)

# Perform Tukey's HSD test
tukey = pairwise_tukeyhsd(endog=df_full_no_morph_bedload['value'], groups=df_full_no_morph_bedload['group'], alpha=0.05)
print(tukey)

# Save the results to a .txt file
path_report_txt = os.path.join(report_dir, 'ANOVA_tukey_full_no_morph_bedload_discharge.txt')
with open(path_report_txt, 'w') as file:
    file.write("No morphological bedload intensity values - full run duration:\n")
    file.write(f"Generated by: {os.path.basename(__file__)}\n")
    file.write("ANOVA Table:\n")
    file.write(anova_table.to_string())
    file.write("\n\nTukey HSD Test:\n")
    file.write(tukey.summary().as_text())


#%%
# =============================================================================
# ANOVA STATISTICAL ANALYSIS - MORPHOLOGICAL BEDLOAD INTENSITY
# =============================================================================
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# LOAD DATA
q05_1_morphological_bedload = np.load(os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q05_1',  'q05_1' + '_values_bedload_intensity_morpho_list.npy')).flatten()
q07_1_morphological_bedload = np.load(os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q07_1',  'q07_1' + '_values_bedload_intensity_morpho_list.npy')).flatten()
q10_2_morphological_bedload = np.load(os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q10_2',  'q10_2' + '_values_bedload_intensity_morpho_list.npy')).flatten()
q15_2_morphological_bedload = np.load(os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q15_2',  'q15_2' + '_values_bedload_intensity_morpho_list.npy')).flatten()
q20_2_morphological_bedload = np.load(os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q20_2',  'q20_2' + '_values_bedload_intensity_morpho_list.npy')).flatten()


# TODO TEST use percentiles
percentiles = np.arange(1, 101)  # Percentiles from 1 to 100

q05_1_morphological_bedload = np.percentile(q05_1_morphological_bedload, percentiles)
q07_1_morphological_bedload = np.percentile(q07_1_morphological_bedload, percentiles)
q10_2_morphological_bedload = np.percentile(q10_2_morphological_bedload, percentiles)
q15_2_morphological_bedload = np.percentile(q15_2_morphological_bedload, percentiles)
q20_2_morphological_bedload = np.percentile(q20_2_morphological_bedload, percentiles)


# Prepare the data
data_morphological_bedload = {
    'value': np.concatenate([q05_1_morphological_bedload, q07_1_morphological_bedload, q10_2_morphological_bedload, q15_2_morphological_bedload, q20_2_morphological_bedload]),
    'group': (['Q=0.5 l/s'] * len(q05_1_morphological_bedload) +
              ['Q=0.7 l/s'] * len(q07_1_morphological_bedload) + 
              ['Q=1.0 l/s'] * len(q10_2_morphological_bedload) + 
              ['Q=1.5 l/s'] * len(q15_2_morphological_bedload) + 
              ['Q=2.0 l/s'] * len(q20_2_morphological_bedload))
}
df_morphological_bedload = pd.DataFrame(data_morphological_bedload)

# Perform ANOVA
model = ols('value ~ C(group)', data=df_morphological_bedload).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print(anova_table)

# Perform Tukey's HSD test
tukey = pairwise_tukeyhsd(endog=df_morphological_bedload['value'], groups=df_morphological_bedload['group'], alpha=0.05)
print(tukey)

# Save the results to a .txt file
path_report_txt = os.path.join(report_dir, 'ANOVA_tukey_morphological_bedload_discharge.txt')
with open(path_report_txt, 'w') as file:
    file.write("Morphological bedload intensity values:\n")
    file.write(f"Generated by: {os.path.basename(__file__)}\n")
    file.write("ANOVA Table:\n")
    file.write(anova_table.to_string())
    file.write("\n\nTukey HSD Test:\n")
    file.write(tukey.summary().as_text())

#%%
# =============================================================================
# ANOVA STATISTICAL ANALYSIS - TOTAL BEDLOAD INTENSITY
# =============================================================================
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# LOAD DATA
q05_1_total_bedload = np.load(os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q05_1',  'q05_1' + '_values_bedload_intensity_list.npy')).flatten()
q07_1_total_bedload = np.load(os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q07_1',  'q07_1' + '_values_bedload_intensity_list.npy')).flatten()
q10_2_total_bedload = np.load(os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q10_2',  'q10_2' + '_values_bedload_intensity_list.npy')).flatten()
q15_2_total_bedload = np.load(os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q15_2',  'q15_2' + '_values_bedload_intensity_list.npy')).flatten()
q20_2_total_bedload = np.load(os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', 'q20_2',  'q20_2' + '_values_bedload_intensity_list.npy')).flatten()

# TODO TEST use percentiles
percentiles = np.arange(1, 101)  # Percentiles from 1 to 100

q05_1_total_bedload = np.percentile(q05_1_total_bedload, percentiles)
q07_1_total_bedload = np.percentile(q07_1_total_bedload, percentiles)
q10_2_total_bedload = np.percentile(q10_2_total_bedload, percentiles)
q15_2_total_bedload = np.percentile(q15_2_total_bedload, percentiles)
q20_2_total_bedload = np.percentile(q20_2_total_bedload, percentiles)

# Prepare the data
data_total_bedload = {
    'value': np.concatenate([q05_1_total_bedload, q07_1_total_bedload, q10_2_total_bedload, q15_2_total_bedload, q20_2_total_bedload]),
    'group': (['Q=0.5 l/s'] * len(q05_1_total_bedload) +
              ['Q=0.7 l/s'] * len(q07_1_total_bedload) + 
              ['Q=1.0 l/s'] * len(q10_2_total_bedload) + 
              ['Q=1.5 l/s'] * len(q15_2_total_bedload) + 
              ['Q=2.0 l/s'] * len(q20_2_total_bedload))
}
df_total_bedload = pd.DataFrame(data_total_bedload)

# Perform ANOVA
model = ols('value ~ C(group)', data=df_total_bedload).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print(anova_table)

# Perform Tukey's HSD test
tukey = pairwise_tukeyhsd(endog=df_total_bedload['value'], groups=df_total_bedload['group'], alpha=0.05)
print(tukey)

# Save the results to a .txt file
path_report_txt = os.path.join(report_dir, 'ANOVA_tukey_total_bedload_discharge.txt')
with open(path_report_txt, 'w') as file:
    file.write("total bedload intensity values:\n")
    file.write(f"Generated by: {os.path.basename(__file__)}\n")
    file.write("ANOVA Table:\n")
    file.write(anova_table.to_string())
    file.write("\n\nTukey HSD Test:\n")
    file.write(tukey.summary().as_text())

# %%
# =============================================================================
#   METRICS RATIO PLOTS
# =============================================================================
report_dir = os.path.join(os.getcwd(), 'output')
output_folder = os.path.join(report_dir, 'topographic_PiQs_overlapping')
run_dir = os.path.join(os.getcwd(), 'surveys')

titles = ["full / total no morph bedload", "total no morph bedload / total bedload", "full no morph bedload / total bedload",
          "total no morph bedload / MAW", "full no morph bedload / MAW", "1/8th no morph bedload / total bedlaod", "1/8th no morph bedload / total no morph bedload"]

for set_name in set_names:
    # LOAD DATA ---------------------------------------------------------------
    plot_dir_out = os.path.join(os.getcwd(),'plots', 'DEM_DoD_envBAA_overlapping', set_name)
    ratio1_list_mean = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_ratio1_list_mean.txt'), delimiter='\t')
    ratio2_list_mean = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_ratio2_list_mean.txt'), delimiter='\t')
    ratio3_list_mean = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_ratio3_list_mean.txt'), delimiter='\t')
    ratio4_list_mean = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_ratio4_list_mean.txt'), delimiter='\t')
    ratio5_list_mean = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_ratio5_list_mean.txt'), delimiter='\t')
    ratio6_list_mean = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_ratio6_list_mean.txt'), delimiter='\t')
    ratio7_list_mean = np.loadtxt(os.path.join(
        plot_dir_out,  set_name + '_ratio7_list_mean.txt'), delimiter='\t')

    data = [ratio1_list_mean, ratio2_list_mean, ratio3_list_mean,
            ratio4_list_mean, ratio5_list_mean, ratio6_list_mean, ratio7_list_mean]
    
    # PLOT DATA
    # Define y-axis limits for each chart
    y_limits = [
        (0, 1),  # Limits for Chart 1
        (0, 1),  # Limits for Chart 2
        (0, 0.5),  # Limits for Chart 3
        (0, 3),  # Limits for Chart 4
        (0, 1.2),  # Limits for Chart 5
        (0, 0.5),  # Limits for Chart 6
        (0, 0.5)   # Limits for Chart 7
    ]

    # Create subplots
    fig, axes = plt.subplots(nrows=7, ncols=1, figsize=(10, 20))

    # Plot each array on a separate chart
    for i, ax in enumerate(axes):
        mean = data[i][0]
        stdev = data[i][1]
        ax.errorbar(DEM_class_labels, mean, yerr=stdev,
                    fmt='o', capsize=5, label='Mean with stdev')
        ax.set_title(titles[i])
        ax.set_ylim(y_limits[i])
        ax.set_xlabel('DEM Class Labels')
        ax.set_ylabel('Values')
        ax.legend()

        # Set the overall title
        fig.suptitle(set_name, fontsize=16, y=0.90)
    
    plt.text(0.5, -0.3, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    
    plt.tight_layout(rect=[0, 0, 1, 0.9])
    plt.savefig(os.path.join(plot_dir_out,
                set_name + '_ratios_DEM_classes.pdf'), dpi=800, bbox_inches='tight')
    plt.show()


# =============================================================================
# PLOT THE RATIOS METRICS WITH THE DISCHARGE
# =============================================================================
# MERGE THE SINGLE REPORT -----------------------------------------------------
files_ratio = ["ratio1", "ratio2", "ratio3",
               "ratio4", "ratio5", "ratio6", "ratio7"]
for ratio in files_ratio:
    # Define the filenames and the corresponding labels
    file_labels = {
        os.getcwd() + '/plots' + '/DEM_DoD_envBAA_overlapping' + "/q05_1/q05_1_" + ratio + "_list_comprehensive_mean.txt": "0.5",
        os.getcwd() + '/plots' + '/DEM_DoD_envBAA_overlapping' + "/q07_1/q07_1_" + ratio + "_list_comprehensive_mean.txt": "0.7",
        os.getcwd() + '/plots' + '/DEM_DoD_envBAA_overlapping' + "/q10_2/q10_2_" + ratio + "_list_comprehensive_mean.txt": "1.0",
        os.getcwd() + '/plots' + '/DEM_DoD_envBAA_overlapping' + "/q15_2/q15_2_" + ratio + "_list_comprehensive_mean.txt": "1.5",
        os.getcwd() + '/plots' + '/DEM_DoD_envBAA_overlapping' + "/q20_2/q20_2_" + ratio + "_list_comprehensive_mean.txt": "2.0"
    }

    files = list(file_labels.keys()) # List of files to process    
    data = {label: [] for label in file_labels.values()} # Initialize an empty dictionary to store the data

    # Read each file and extract all values
    for file in files:
        label = file_labels[file]
        with open(file, 'r') as f:
            lines = f.readlines()
            # Assume all values except the first line are mean values
            mean_values = [line.strip() for line in lines[1:]]
            data[label].extend(mean_values)

    # Write the merged data to a new file
    output_file = os.path.join(output_folder, ratio + "_merged_data.txt")
    with open(output_file, 'w') as f_out:
        # Write the header
        f_out.write(", ".join(data.keys()) + "\n")

        # Write the data
        max_length = max(len(data[label]) for label in data)
        for i in range(max_length):
            row_values = [data[label][i] if i < len(
                data[label]) else "" for label in data]
            f_out.write(", ".join(row_values) + "\n")

# LOAD THE MERGED REPORT
ratio1_comprehensive = np.loadtxt(
    '/home/erri/Documents/PhD/Research/5_research_repos/DoD_analysis/output/topographic_PiQs_overlapping/ratio1_merged_data.txt', delimiter=',')
ratio2_comprehensive = np.loadtxt(
    '/home/erri/Documents/PhD/Research/5_research_repos/DoD_analysis/output/topographic_PiQs_overlapping/ratio2_merged_data.txt', delimiter=',')
ratio3_comprehensive = np.loadtxt(
    '/home/erri/Documents/PhD/Research/5_research_repos/DoD_analysis/output/topographic_PiQs_overlapping/ratio3_merged_data.txt', delimiter=',')
ratio4_comprehensive = np.loadtxt(
    '/home/erri/Documents/PhD/Research/5_research_repos/DoD_analysis/output/topographic_PiQs_overlapping/ratio4_merged_data.txt', delimiter=',')
ratio5_comprehensive = np.loadtxt(
    '/home/erri/Documents/PhD/Research/5_research_repos/DoD_analysis/output/topographic_PiQs_overlapping/ratio5_merged_data.txt', delimiter=',')
ratio6_comprehensive = np.loadtxt(
    '/home/erri/Documents/PhD/Research/5_research_repos/DoD_analysis/output/topographic_PiQs_overlapping/ratio6_merged_data.txt', delimiter=',')
ratio7_comprehensive = np.loadtxt(
    '/home/erri/Documents/PhD/Research/5_research_repos/DoD_analysis/output/topographic_PiQs_overlapping/ratio7_merged_data.txt', delimiter=',')
data = {"full / total no morph bedload": np.array(ratio1_comprehensive),
        "total no morph bedload / total bedload": np.array(ratio2_comprehensive),
        "full no morph bedload / total bedload": np.array(ratio3_comprehensive),
        "total no morph bedload / MAW": np.array(ratio4_comprehensive),
        "full no morph bedload / MAW": np.array(ratio5_comprehensive),
        "1/8th no morph bedload / total bedlaod": np.array(ratio6_comprehensive),
        "1/8th no morph bedload / total no morph bedload": np.array(ratio7_comprehensive)
        }
titles = ["full / total no morph bedload", "total no morph bedload / total bedload", "full no morph bedload / total bedload",
          "total no morph bedload / MAA", "full no morph bedload / MAA", "1/8th no morph bedload / total bedlaod", "1/8th no morph bedload / total no morph bedload"]

# Plot each dataset
for key, dataset in data.items():
    x_values = dataset[0]
    y_values = dataset[1]
    stdev = dataset[2]
    plt.errorbar(x_values, y_values, yerr=stdev, fmt='o-', label=key)
plt.xlabel('Discharge [l/s]')
# plt.ylabel('Y values')
plt.xticks(data["full / total no morph bedload"][0])
# plt.title('Ratio')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.grid(True)
plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
         transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
plt.savefig(os.path.join(plot_dir,
            'REPORT_ratio_discharge.pdf'), dpi=800, bbox_inches='tight')
plt.show()


#%%============================================================================
# COMPENSATION AND FILTERING EFFECTS ANALYSIS
# =============================================================================
import numpy as np
import os
import matplotlib.pyplot as plt

home_dir = os.getcwd()  # Home directory
report_dir = os.path.join(home_dir, 'output')
output_folder = os.path.join(report_dir, 'topographic_PiQs_overlapping')

'''
The aim of this section is analyze the mean bedload activity where compensation
and filtering process affect our morphological changes detection capabilities.
Furthermore we repeat the same analysis considering number of bedlaod active
periods.
'''
for set_name in ['q07_1', 'q10_2', 'q15_2', 'q20_2']:
# for set_name in ['q07_1']:

    # IMPORT DATA -------------------------------------------------------------
    envBAA_maps                  = np.load(os.path.join(report_dir, set_name + '_envBAA_1Txnr_stack.npy')) # Bedload intensity
    envBAA_map_activity_periods = np.load(os.path.join(report_dir, set_name + '_envBAA_active_periods_1Txnr_stack.npy')) # Bedload active periods
    envBAA_comp                 = np.load(os.path.join(report_dir, set_name + '_envBAA_comp_stack.npy')) # Bedload intensity where comp occurs
    envBAA_thrs                 = np.load(os.path.join(report_dir, set_name + '_envBAA_thrs_stack.npy')) # Bedload intensity where filt process occurs
    envBAA_periods_comp         = np.load(os.path.join(report_dir, set_name + '_envBAA_periods_comp_stack.npy')) # Bedload activity periods where comp occurs
    envBAA_periods_thrs         = np.load(os.path.join(report_dir, set_name + '_envBAA_periods_thrs_stack.npy')) # Bedload activity periods filtering process occurs

    # COMBINE ALL LAYERS INTO A SINGLE ARRAY FOR EACH DATASET------------------
    # Remove np.nana and 0 values
    envBAA_maps_combined = envBAA_maps.flatten()
    envBAA_maps_combined = envBAA_maps_combined[~np.isnan(envBAA_maps_combined) & (envBAA_maps_combined != 0)]
    envBAA_comp_combined = envBAA_comp.flatten()
    envBAA_comp_combined = envBAA_comp_combined[~np.isnan(envBAA_comp_combined) & (envBAA_comp_combined != 0)]
    envBAA_thrs_combined = envBAA_thrs.flatten()
    envBAA_thrs_combined = envBAA_thrs_combined[~np.isnan(envBAA_thrs_combined) & (envBAA_thrs_combined != 0)]
    envBAA_periods_comp_combined = envBAA_periods_comp.flatten()
    envBAA_periods_comp_combined = envBAA_periods_comp_combined[~np.isnan(envBAA_periods_comp_combined) & (envBAA_periods_comp_combined != 0)]
    envBAA_periods_thrs_combined = envBAA_periods_thrs.flatten()
    envBAA_periods_thrs_combined = envBAA_periods_thrs_combined[~np.isnan(envBAA_periods_thrs_combined) & (envBAA_periods_thrs_combined != 0)]
    
# =============================================================================
# COMPUTE DISTRIBUTION OF FREQUENCY BARS PLOT - BEDLOAD INTENSITY
# =============================================================================

    num_bins = 30 # Define the number of bins for the histograms
    value_range = (0, 30) # Define the range of values for the histograms
    # Compute histograms using np.histogram
    envBAA_hist, envBAA_bins = np.histogram(envBAA_maps_combined, bins=num_bins, range=value_range)
    envBAA_comp_hist, envBAA_comp_hist_bins = np.histogram(envBAA_comp_combined, bins=envBAA_bins, range=value_range)
    envBAA_thrs_hist, envBAA_thrs_hist_bins = np.histogram(envBAA_thrs_combined, bins=envBAA_bins, range=value_range)
    # Normalize 
    envBAA_hist_norm_single         = envBAA_hist/envBAA_hist
    envBAA_comp_hist_norm_single    = envBAA_comp_hist/envBAA_hist
    envBAA_thrs_hist_norm_single    = envBAA_thrs_hist/envBAA_hist
    remaining_envBAA_hist_norm_single = envBAA_hist_norm_single - envBAA_comp_hist_norm_single - envBAA_thrs_hist_norm_single
    
    # PLOT HISTOGRAM WITH BARS
    plt.figure(figsize=(12, 8))
    # width = (envBAA_bins[1] - envBAA_bins[0])
    width = 0.6
    plt.bar(envBAA_bins[:-1], envBAA_thrs_hist_norm_single, width=width, alpha=0.6, label='envBAA Thrs Data', color='red')
    plt.bar(envBAA_bins[:-1], envBAA_comp_hist_norm_single, width=width, alpha=0.6, label='envBAA Comp Data', color='green', bottom=envBAA_thrs_hist_norm_single)
    plt.bar(envBAA_bins[:-1], remaining_envBAA_hist_norm_single, width=width, alpha=0.6, label='Remaining envBAA Layer Data', color='blue', bottom=envBAA_thrs_hist_norm_single+envBAA_comp_hist_norm_single)
    plt.xlim(left=0, right=35)
    plt.ylim(bottom=0, top=1.1)
    plt.title(set_name + ' - Bedload intensity distribution of frequency.\nTotal bedload, compensation and threshold ')
    plt.xlabel('Bedload intensity')
    plt.ylabel('Normalized Frequency')
    plt.legend()
    plt.text(0.5, -0.1, f"Generated by: {os.path.basename(__file__)}",
             transform=plt.gca().transAxes, ha='center', va='center', fontsize=10)
    plt.savefig(os.path.join(plot_dir_out, set_name + '_comp_thrs_bedload_intensity_bars.pdf'), dpi=800, bbox_inches='tight')
    plt.show()

# =============================================================================
# COMPUTE DISTRIBUTION OF FREQUENCY BARS PLOT - BEDLOAD ACT PERIODS
# =============================================================================
# TODO FINISH ACTG PERIODS HISTOGRAM
    num_bins = 30 # Define the number of bins for the histograms
    value_range = (0, 30) # Define the range of values for the histograms
    # Compute histograms using np.histogram
    envBAA_periods_comp_combined_hist, envBAA_periods_comp_combined_hist_bins = np.histogram(envBAA_periods_comp_combined, bins=envBAA_bins, range=value_range)
    envBAA_periods_thrs_combined_hist, envBAA_periods_thrs_combined_hist_bins = np.histogram(envBAA_periods_thrs_combined, bins=envBAA_bins, range=value_range)
    # Normalize histogram values
    envBAA_periods_comp_norm_single = envBAA_periods_comp_combined_hist/envBAA_hist
    envBAA_periods_thrs_norm_single = envBAA_periods_thrs_combined_hist/envBAA_hist