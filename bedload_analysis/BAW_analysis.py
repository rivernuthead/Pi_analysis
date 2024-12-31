# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 22:55:12 2024

@author: erri
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# RUNS
run_names = ['q05rgm4', 'q07rgm', 'q10rgm2', 'q15rgm2', 'q20rgm2']
# run_names = ['q05_rgm4']



w_dir = 'e:\\repos\\pi_analysis'
cleaned_stack_folder = os.path.join('e:\\repos\\pi_analysis', 'bedload_analysis', 'output_data', 'activity_stack', 'activity_stack_cleaned')
run_parameters = pd.read_csv(os.path.join('e:\\repos\\pi_analysis', 'Run_parameters.csv'))

import matplotlib.pyplot as plt
import numpy as np

# Example discharge values for the x-axis of the second plot
discharge_values = [0.5, 0.7, 1.0, 1.5, 2.0]

# Lists to store mean and std values
means = []
stds = []

# Create the figure and subplots with shared y-axis and reduced width for the right plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4), sharey=True, 
                                gridspec_kw={'width_ratios': [3, 1]})  # Adjusted width ratio

# Left plot: Time series
for run_name, discharge in zip(run_names, discharge_values):  # Assuming run_names match the discharge values
    print(run_name, 'is running...')
    
    # EXTRACT RUN PARAMETERS --------------------------------------------------
    selected_row = run_parameters.loc[run_parameters['run_code'] == run_name]
    if not selected_row.empty:
        disch_label = selected_row['disch label'].values[0]  # Discharge label
        exner_time = selected_row['Exner time [min]'].values[0]  # Exner time

    stack_cleaned = np.load(os.path.join(cleaned_stack_folder, run_name + '_BAA_stack_LR5_cld.npy'))
    stack_cleaned_bool = np.where(stack_cleaned > 0, 1, 0)
    
    BAW_array = np.sum((np.sum(stack_cleaned_bool, axis=1)), axis=1) / (120 * stack_cleaned.shape[2])
    BAW_array_time = np.linspace(0, stack_cleaned.shape[0], stack_cleaned.shape[0]) / exner_time
    
    # Filter for the last 8.5 minutes
    time_threshold = max(BAW_array_time) - 8.5
    valid_indices = BAW_array_time >= time_threshold
    filtered_time = BAW_array_time[valid_indices]
    filtered_BAW = BAW_array[valid_indices]
    filtered_time = filtered_time - np.min(filtered_time)
    print(f"Plotting from t = {min(filtered_time):.2f} to t = {max(filtered_time):.2f}")
    
    # Compute mean and std
    mean_value = np.mean(filtered_BAW)
    std_value = np.std(filtered_BAW)
    means.append(mean_value)
    stds.append(std_value)
    
    # Plot on the left subplot
    ax1.plot(filtered_time, filtered_BAW, linestyle='--', label=disch_label, linewidth=1.5, markersize=4)

# Labels, title, legend, and customization for the time series
ax1.set_xlabel('t / Exner time [-]', fontsize=10)
ax1.set_ylabel('Normalized active width [-]', fontsize=10)
ax1.set_title('Instantaneous active width', fontsize=12)
ax1.legend(bbox_to_anchor=(0, 1), loc='upper left', fontsize=8)
ax1.set_ylim(0.0, 1.0)
ax1.grid(True, linestyle='--', alpha=0.7)
ax1.tick_params(axis='both', which='major', labelsize=10)

# Right plot: Mean and Std Dev
ax2.errorbar(discharge_values, means, yerr=stds, fmt='o', linestyle='--', color='grey', ecolor='black', 
             capsize=5, label='Mean ± Std Dev')
ax2.set_xlabel('Discharge [-]', fontsize=10)
ax2.set_title('Bedload active width', fontsize=12)
ax2.grid(True, linestyle='--', alpha=0.7)
ax2.tick_params(axis='both', which='major', labelsize=10)
ax2.legend(fontsize=8)

# Adjust layout and show the combined plots
plt.tight_layout()
plt.show()
