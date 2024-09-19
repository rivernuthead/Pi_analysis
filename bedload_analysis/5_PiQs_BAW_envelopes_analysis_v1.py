#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 15:44:32 2023

@author: erri

The aims of this script are:
    1. Considering the stack input mode, it takes as input activity maps from
    PiQs script as stack of activity map from "to regime" runs or merged single
    runs activity maps. Then it use the stack to compute all the available
    envelopes at increasing timescale. Timescales are given by
    the env_tscale_array at power of two: for example overy 2, 4, 8, 16 and so
    on minutues...
    2. It takes as input DoDs from the DoD_analysis repository and extract DoDs
    at the finest discretization. The script use the DoDs to compute all
    available envelopes increasing the number of DoDs involved.
    3. The script plots the comparison between envBAW and envMAW with the
    following strategy:
        - stack_input = 'merged_single_run' produces plots for q07, q10 and q20
          where envBAW comes from merged single runs envelopes.
        - stack_input = 'regime_run' produces plot for q15 because it is the only
          run where there is no correspondence between envMAW (q15_3) and envBAW
          (q15rgm2).
    4. The script produces plots of the instantaneous BAW and the comulative
    envBAW over time.
    
"""

# IMPORT PACKAGES ------------------------------------------------------------#
import numpy as np
import matplotlib.pyplot as plt
import os
from fractions import Fraction

# FUNCTIONS--------------------------------------------------------------------

def compute_multiples_of_2(start, end):
    multiples = []
    current = start
    while current < end:
        multiples.append(current)
        current *= 2
    return multiples

def stack_to_envelopes_stack(stack, time_discretization):
    """
    Reshapes and processes a stack of data into envelopes based on time discretization.

    Parameters:
    stack (np.ndarray): The input stack with shape (time, height, width).
    time_discretization (int): The number of time frames to group together.

    Returns:
    tuple: A tuple containing:
        - reshaped_stack (np.ndarray): The reshaped stack with shape (n_groups, time_discretization, height, width).
        - reshaped_sum_stack (np.ndarray): The sum of each group in the reshaped stack with shape (n_groups, height, width).
    
    Raises:
    TypeError: If the number of groups (n_groups) is less than or equal to zero.
    """
    
    import numpy as np
    
    # Calculate the number of groups we can form
    n_groups = stack.shape[0] // time_discretization
    
    # Check if the number of groups is valid
    if n_groups <= 0:
        raise TypeError("Stack time dimension has to be greater than or at least equal to the time discretization")

    else:
        # Trim the stack to make its total number of frames a multiple of time_discretization
        trimmed_stack = stack[:n_groups * time_discretization, :, :]
        
        # Reshape the trimmed stack into n_groups of stacks, each with time_discretization maps
        reshaped_stack = trimmed_stack.reshape(
            n_groups, time_discretization, stack.shape[1], stack.shape[2])
        
        # Sum each stack along the time_discretization dimension to obtain envelopes
        reshaped_sum_stack = np.nansum(reshaped_stack, axis=1)
        
        return reshaped_stack, reshaped_sum_stack


# DEFINE FOLDERS --------------------------------------------------------------
folder_home = os.getcwd()  # Setup home folder

# SELECT RUN ------------------------------------------------------------------

runs = ['05', '07', '10', '15', '20']

runs = ['05']


# ANALYSIS PARAMETERS ---------------------------------------------------------

# SCRIPT MODE -----------------------------------------------------------------
plt_show = 1
# stack_input = 'merged_single_run'
# stack_input = 'regime_run'
stack_input = 'single run'

# SCRIPT PARAMETERS -----------------------------------------------------------
downsampling_dim = 5

# EXNER TIME ------------------------------------------------------------------
Txnr_array = [20, 30, 48, 74, 100]

# RUN TIME --------------------------------------------------------------------
run_time_array = [16, 21, 28, 47, 81]  # Run duration in minutes

# RUN DURATION ----------------------------------------------------------------
run_frames = [10, 15, 20, 40, 88]  # Run duration in number of frames


# TIMESCALE -------------------------------------------------------------------
# env_tscale_array = [5,5,5,5]
# env_tscale_array = [1,1,1,1]
env_tscale_array = [2, 3, 4, 9, 22]  # 1/4 of the frame number (~1/8*Txnr)
# env_tscale_array = [3,4,6,9]        # 1/8 Txnr
# env_tscale_array = [5,7,12,18]      # 1/4 Txnr
# env_tscale_array = [10,15,24,37]    # 1/2 Txnr
# env_tscale_array = [15,22,36,55]    # 3/4 Txnr
# env_tscale_array = [20,30,48,74]    # 1 Txnr
# env_tscale_array = [25,37,60,92]    #5/4 Txnr
# env_tscale_array = [30,45,72,110]   #6/4 Txnr
# env_tscale_array = [35,52,84,129]   #7/4 Txnr
# env_tscale_array = [40,59,96,147]
# env_tscale_array = [46, 67, 108, 166]
# env_tscale_array = [51,74,120,184]
# env_tscale_array = [56,82,132,202]
# env_tscale_array = [61,89,144,221]
# env_tscale_array = [66,97,156,239]


# =============================================================================
# LOOP OVER THE RUNS
# =============================================================================
for run in runs:
    if '05' in run:
        skip       = 1
        env_tscale = env_tscale_array[4]
        Txnr       = Txnr_array[4]
        run_time   = run_time_array[4]
        run_frm    = run_frames[4]
        set_name = 'q05_1'
        rgm_code = 'q05rgm4'
        rgm_skip_frame = 0
        MAW_set_name = 'q05_1'
    
    if '07' in run:
        skip       = 1
        env_tscale = env_tscale_array[3]
        Txnr       = Txnr_array[3]
        run_time   = run_time_array[3]
        run_frm    = run_frames[3]
        set_name = 'q07_1'
        rgm_code = 'q07rgm'
        rgm_skip_frame = 0
        MAW_set_name = 'q07_1'
        
    if '10' in run:
        skip       = 1
        env_tscale = env_tscale_array[2]
        Txnr       = Txnr_array[2]
        run_time   = run_time_array[2]
        run_frm    = run_frames[2]
        set_name = 'q10_2'
        rgm_code = 'q10rgm2'
        rgm_skip_frame = 0
        MAW_set_name = 'q10_2'
        
    if '15' in run:
        skip       = 1
        env_tscale = env_tscale_array[1]
        Txnr       = Txnr_array[1]
        run_time   = run_time_array[1]
        run_frm    = run_frames[1]
        set_name = 'q15_2'
        rgm_code = 'q15rgm2'
        rgm_skip_frame = 5
        MAW_set_name = 'q15_3'
        
    if '20' in run:
        skip       = 1
        env_tscale = env_tscale_array[0]
        Txnr       = Txnr_array[0]
        run_time   = run_time_array[0]
        run_frm    = run_frames[0]
        set_name = 'q20_2'
        rgm_code = 'q20rgm2'
        rgm_skip_frame = 0
        MAW_set_name = 'q20_2'


    print('****************')
    print(run, '  Timescale: ', env_tscale)
    print('****************')
    
    # =============================================================================
    # DEFINE INPUT FILE
    # =============================================================================
    if 'merged_single_run' in stack_input:
        stack_path = os.path.join(os.getcwd(), 'activity_stack/activity_stack_cleaned', set_name +
                '_single_run_merged_envelope_BAA_stack_LR' + str(downsampling_dim) + '_cld.npy')

    if 'regime_run' in stack_input:
        stack_path = os.path.join(os.getcwd(), 'activity_stack/activity_stack_cleaned', rgm_code +
                '_envelope_BAA_stack_LR' + str(downsampling_dim) + '_cld.npy')

    if 'single_run' in stack_input:
        for run in runs:
            stack_path = os.path.join(os.getcwd(), 'activity_stack/activity_stack_cleaned', set_name +
                    '_envelope_BAA_stack_LR' + str(downsampling_dim) + '_cld.npy')













    
    # CREATE FOLDERS ---------------------------------------------------------#
    # Create the directory for single sums if it does not exist
    plots_dir = os.path.join(folder_home, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    # Load the numpy stack
    print('REMEMBER TO CHECK TE INPUT DATA PATH \n \n')

    if stack_input == 'merged_single_run':
        # 2D + t filtered stack as sum of all the single runs
        stack_path = os.path.join(os.getcwd(), 'activity_stack/activity_stack_cleaned', set_name + '_single_run_merged_BAA_stack_LR' + str(downsampling_dim) + '_cld.npy')
        stack = np.load(stack_path)
    elif stack_input == 'regime_run':
        # Stack from "to regime" run
        stack_path = os.path.join(os.getcwd(), 'activity_stack/activity_stack_cleaned', rgm_code + '_BAA_stack_LR' + str(downsampling_dim) + '_cld.npy')
        stack = np.load(stack_path)
        stack = stack[rgm_skip_frame:,:,:]
    elif stack_input == 'single_runs':
    # CONVERT STACK IN BOOL
    stack_bool = (stack>0)*1
    

    # Check the shape of the loaded array
    print("Original shape:", stack_bool.shape)
    print()

    # TEST -------------------------------------------------------------------#
    # stack_test = np.nansum(stack_bool[0:256, :, :], axis=0)
    # stack_test_bool = np.where(stack_test > 0, 1, 0)
    # BAA = np.nansum(stack_test_bool)
    # BAW_test = BAA/np.nansum(np.logical_not(np.isnan(stack_bool[0, :, :])))

    

    # PRELIMINARY ANALYSIS ---------------------------------------------------#
    # Compute the BAW
    BAA_ist = np.nansum(stack_bool, axis=1)/120
    BAW_ist_mean_array = np.nanmean(BAA_ist, axis=1) # Array of instantaneous BAW over time
    # Save instantaneous BAW over time
    if stack_input == 'merged_single_run':
        np.savetxt(os.path.join(plots_dir, run+'_single_run_merged_instantaneous_BAW.txt'), BAW_ist_mean_array, delimiter=',', fmt='%.4f')
    elif stack_input == 'regime_run':
        np.savetxt(os.path.join(plots_dir, run+'_regime_instantaneous_BAW.txt'), BAW_ist_mean_array, delimiter=',', fmt='%.4f')
    BAW_ist_mean = np.mean(BAW_ist_mean_array)
    BAW_ist_std = np.std(BAW_ist_mean_array)

    # Compute the envBAW
    cumulative_env_stack_cld = np.cumsum(stack_bool, axis=0)
    cumulative_env_stack_bool_cld = np.where(cumulative_env_stack_cld > 0, 1, 0)
    cumulative_envBAW = np.nansum(cumulative_env_stack_bool_cld, axis=1)/120
    cumulative_envBAW = np.nanmean(cumulative_envBAW, axis=1)
    
    # Save instantaneous envBAW over time
    if stack_input == 'merged_single_run':
        np.savetxt(os.path.join(plots_dir, run+'_single_run_merged_instantaneous_envBAW.txt'), cumulative_envBAW, delimiter=',', fmt='%.4f')
    elif stack_input == 'regime_run':
        np.savetxt(os.path.join(plots_dir, run+'_regime_instantaneous_envBAW.txt'), cumulative_envBAW, delimiter=',', fmt='%.4f')

    env_tscale_values = compute_multiples_of_2(env_tscale, stack_bool.shape[0])

    # BAW_x_values = np.array(env_tscale_x_values)/min(np.array(env_tscale_x_values))*1/8
    BAW_x_values = np.insert(np.array(env_tscale_values), 0, 0)/Txnr

    BAW_mean = [BAW_ist_mean]
    BAW_std = [BAW_ist_std]

    # Convert mean_x values to fractions
    fractions_x = [Fraction.from_float(
        val).limit_denominator() for val in BAW_x_values]
    # Format fractions as strings for x-tick labels
    xtick_labels = [f"{frac.numerator}/{frac.denominator}" if frac.denominator !=
                    1 else str(frac.numerator) for frac in fractions_x]

    for env_tscale, frac in zip(env_tscale_values, xtick_labels[1:]):

        print('Number of shots: ', env_tscale, '    ', frac, 'Exner Time')

        # Reshape the array to group every env_tscale matrices (140x1256)
        # The new shape will be (n_groups, Tscale, 140, 1256)
        n_groups = stack_bool.shape[0] // env_tscale
        print('Number of groups:', n_groups)
        # Trim the stack if the total number of frames exceeds n_groups*Tscale
        trimmed_stack = stack_bool[:n_groups * env_tscale,:,:]
        
        # Reshape the stack creating n_group number of stacks with a number of
        # Tscale maps each
        reshaped_stack = trimmed_stack.reshape(
            n_groups, env_tscale, stack_bool.shape[1], stack_bool.shape[2])

        # Sum every single stack to obtain n_group-number of envelopes
        envBAA_stack = np.nansum(reshaped_stack, axis=1)
        # Convert stack in bool matrix
        envBAA_stack_bool = np.where(envBAA_stack > 0, 1, 0)

        # Check the shape of the resulting array
        print("New shape:", envBAA_stack.shape)
        print()

        # Save the resulting stack
        envBAA_stack_filename = 'envBAA_stack_LR5_cld_'+str(env_tscale)+'.npy'
        np.save(envBAA_stack_filename, envBAA_stack)
        # print(f'Saved summed stack as {envBAA_stack_filename}')

        # # Create the directory for single sums if it does not exist
        # single_envelopes_dir = 'single_envelopes/' + \
        #     'env_tascle_'+str(env_tscale)+'/'+run
        # os.makedirs(single_envelopes_dir, exist_ok=True)

        # # Save each summed matrix individually
        # for i, single_sum in enumerate(envBAA_stack):
        #     single_filename = os.path.join(single_envelopes_dir, f'{i+1}_envBAA_stack_LR5_cld.npy')
        #     np.save(single_filename, single_sum)
        #     # print(f'Saved single summed matrix {i+1} as {single_filename}')

        # =====================================================================
        #     ANALYSIS
        # =====================================================================

        envBAW_array = np.nansum(envBAA_stack_bool, axis=1)/120
        envBAW_array_mean = np.mean(envBAW_array, axis=1)
        envBAW_mean = np.mean(envBAW_array_mean)
        envBAW_std = np.std(envBAW_array_mean)


        BAW_mean.append(envBAW_mean)
        BAW_std.append(envBAW_std)

    # SAVE envBAW MEAN AND STDEV IN A .txt FILE -------------------------------
    title = int(run)*0.1
    np.savetxt(os.path.join(plots_dir, run+'_single_run_merged_envBAW_mean_stdev.txt'), np.stack((BAW_x_values,BAW_mean, BAW_std)), fmt='%.4f', delimiter=',',
               header='Txnr, envBAW mean and stdev, ' + 'Q = ' + f"{title:.1f}" + ' l/s' + ' - Generated by: ' + os.path.basename(__file__))
    
    
    # =========================================================================
    #     BAW ENVELOPES PLOT
    # =========================================================================

    plt.figure(figsize=(20, 4))
    plt.errorbar(BAW_x_values, BAW_mean, yerr=BAW_std, fmt='o')
    plt.xlabel('Exner Time [-]')
    plt.ylabel('BAW [-]')
    title = int(run)*0.1
    plt.title('Q = ' + f"{title:.1f}" + ' l/s')
    plt.xticks(BAW_x_values, xtick_labels, fontsize=8)
    plt.tick_params(axis='x', pad=5)  # Adjust padding between ticks and labels
    plt.grid(True)
    plt.tight_layout()
    plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
              transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)

    # Save the plot as a PDF file
    plt.savefig(os.path.join(plots_dir, run +
                '_BAW_at_diff_envelopes.pdf'), dpi=800, bbox_inches='tight')
    plt.show()
    plt.close()


# =============================================================================
#     MAW AND BAW PLOT
# =============================================================================
    # PREPARE THE DATA
    MAW_data = np.loadtxt(os.path.join('/home/erri/Documents/PhD/Research/5_research_repos/DoD_analysis/output/report_' +
                          MAW_set_name+'/envelopes/'+MAW_set_name+'_envelope_timespan0_starting_points_raw.txt'), skiprows=2, delimiter=',')
    MAW_data = np.where(MAW_data == 0, np.nan, MAW_data)
    MAW_data_mean = np.nanmean(MAW_data, axis=1)
    MAW_data_std = np.nanstd(MAW_data, axis=1)
    MAW_x_values = np.arange(1, 10, 1)*run_time/Txnr

    plt.figure(figsize=(20, 4))
    if '15' in set_name:
        plt.errorbar(BAW_x_values, BAW_mean, yerr=BAW_std, fmt='o', label='envBAW - ' + rgm_code)
    else:
        plt.errorbar(BAW_x_values, BAW_mean, yerr=BAW_std, fmt='o', label='envBAW - ' + set_name)    

    plt.errorbar(MAW_x_values, MAW_data_mean,
                  yerr=MAW_data_std, fmt='o', label='envMAW - ' + MAW_set_name)
    # plt.errorbar(np.arange(len(cumulative_envBAW))/Txnr,
    #              cumulative_envBAW, fmt='.', label='raw envBAW\nsingle runs')
    plt.xlabel('Exner Time [-]')
    plt.ylabel('MAW & BAW [-]')
    title = int(run)*0.1
    plt.title('Q = ' + f"{title:.1f}" + ' l/s')
    # plt.xticks(BAW_x_values, xtick_labels, fontsize=8)
    plt.tick_params(axis='x', pad=5)  # Adjust padding between ticks and labels
    plt.grid(True)
    plt.tight_layout()
    # plt.xscale('log', base=2)
    plt.ylim((0.2, 1))
    plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
              transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
    plt.legend()
    # Save the plot as a PDF file
    plt.savefig(os.path.join(plots_dir, run +
                '_MAW_BAW_at_diff_envelopes.pdf'), dpi=800, bbox_inches='tight')
    plt.show()
    plt.close()

# =============================================================================
# PLOT INSTANTANEOUS BAW
# =============================================================================
# Data preparation
BAW_figure, BAW_ax = plt.subplots(figsize=(14, 6))
for run in runs:

    if '07' in run:
        set_name = 'q07_1'
        rgm_code = 'q07rgm'
    if '10' in run:
        set_name = 'q10_2'
        rgm_code = 'q10rgm2'
    if '15' in run:
        set_name = 'q15_2'
        rgm_code = 'q15rgm2'
    if '20' in run:
        set_name = 'q20_2'
        rgm_code = 'q20rgm2'
    
    
    if stack_input == 'merged_single_run':
        BAW_ist_mean_array = np.loadtxt(os.path.join(plots_dir, run+'_single_run_merged_instantaneous_BAW.txt'), delimiter=',')
        BAW_ax.set_title('Instantaneous BAW - Merged single runs')
        plot_path_out = os.path.join(plots_dir, 'single_run_merged_instantaneous_BAW.pdf')
        BAW_ax.plot(np.linspace(0, len(BAW_ist_mean_array), len(BAW_ist_mean_array)), BAW_ist_mean_array, label=set_name)
    elif stack_input == 'regime_run':
        BAW_ist_mean_array = np.loadtxt(os.path.join(plots_dir, run+'_regime_instantaneous_BAW.txt'), delimiter=',')
        BAW_ax.set_title('Instantaneous BAW - Regime run')
        plot_path_out = os.path.join(plots_dir, 'regime_run_instantaneous_BAW.pdf')
        BAW_ax.plot(np.linspace(0, len(BAW_ist_mean_array), len(BAW_ist_mean_array)), BAW_ist_mean_array, label=rgm_code)

    
BAW_ax.set_ylim(0.1,1)
BAW_ax.set_ylabel('BAW [-]')
BAW_ax.set_xlabel('time [Exner]')
BAW_ax.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
          transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
BAW_ax.grid(True)
BAW_ax.legend()

plt.savefig(plot_path_out, dpi=800, bbox_inches='tight')
plt.show()


# =============================================================================
# PLOT INSTANTANEOUS BAW
# =============================================================================
# Data preparation
envBAW_figure, envBAW_ax = plt.subplots(figsize=(14, 6))
for run in runs:
    
    if '07' in run:
        set_name = 'q07_1'
        rgm_code = 'q07rgm'
    if '10' in run:
        set_name = 'q10_2'
        rgm_code = 'q10rgm2'
    if '15' in run:
        set_name = 'q15_2'
        rgm_code = 'q15rgm2'
    if '20' in run:
        set_name = 'q20_2'
        rgm_code = 'q20rgm2'
    
    if stack_input == 'merged_single_run':
        cumulative_envBAW = np.loadtxt(os.path.join(plots_dir, run+'_single_run_merged_instantaneous_envBAW.txt'), delimiter=',')
        envBAW_ax.set_title('Instantaneous envBAW - Merged single runs')
        plot_path_out = os.path.join(plots_dir, 'single_run_merged_instantaneous_envBAW.pdf')
        envBAW_ax.plot(np.linspace(0, len(cumulative_envBAW), len(cumulative_envBAW)), cumulative_envBAW, label=set_name)
    elif stack_input == 'regime_run':
        cumulative_envBAW = np.loadtxt(os.path.join(plots_dir, run+'_regime_instantaneous_envBAW.txt'), delimiter=',')
        envBAW_ax.set_title('Instantaneous envBAW - Regime run')
        plot_path_out = os.path.join(plots_dir, 'regime_run_instantaneous_envBAW.pdf')
        envBAW_ax.plot(np.linspace(0, len(cumulative_envBAW), len(cumulative_envBAW)), cumulative_envBAW, label=rgm_code)

    
envBAW_ax.set_ylim(0.1,1)
envBAW_ax.set_ylabel('BAW [-]')
envBAW_ax.set_xlabel('time [Exner]')
envBAW_ax.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
          transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
envBAW_ax.grid(True)
envBAW_ax.legend()

plt.savefig(plot_path_out, dpi=800, bbox_inches='tight')
plt.show()