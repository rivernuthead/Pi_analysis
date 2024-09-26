#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 10:05:33 2024

@author: erri
"""

# =============================================================================
# IMPORT PACKAGES 
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import os

# =============================================================================
# FUNCTIONS
# =============================================================================
# Task 1: Function to calculate consecutive 1, 0, -1 periods and store them in separate stacks.
def calculate_consecutive_periods(arr, t, y, x):
    """
    Calculate consecutive periods of 1s, 0s, and -1s in a 3D array and store 
    the lengths in separate stacks.

    Parameters:
    arr (ndarray): 3D array of shape (t, y, x) containing the time series.
    t (int): Number of time steps.
    y (int): Height dimension.
    x (int): Width dimension.

    Returns:
    tuple: Three 3D arrays storing lengths of consecutive 1s, 0s, and -1s.
    """
    ones_stack = np.zeros((t, y, x), dtype=int)
    zeros_stack = np.zeros((t, y, x), dtype=int)
    neg_ones_stack = np.zeros((t, y, x), dtype=int)

    for i in range(y):
        for j in range(x):
            series = arr[:, i, j]
            count = 1
            idx = 0
            for k in range(1, t):
                if series[k] == series[k - 1]:
                    count += 1
                else:
                    if series[k - 1] == 1:
                        ones_stack[idx, i, j] = count
                    elif series[k - 1] == 0:
                        zeros_stack[idx, i, j] = count
                    elif series[k - 1] == -1:
                        neg_ones_stack[idx, i, j] = count
                    count = 1
                    idx += 1

            # Assign the last counted period
            if series[-1] == 1:
                ones_stack[idx, i, j] = count
            elif series[-1] == 0:
                zeros_stack[idx, i, j] = count
            elif series[-1] == -1:
                neg_ones_stack[idx, i, j] = count

    return ones_stack, zeros_stack, neg_ones_stack

# Task 2: Function to calculate the distance between changes in sign
def calculate_distance_between_switches(arr, t, y, x):
    """
    Calculate the distances between changes in sign (1 to -1 or -1 to 1) in the 
    time series and store them in a 3D stack.

    Parameters:
    arr (ndarray): 3D array of shape (t, y, x) containing the time series.
    t (int): Number of time steps.
    y (int): Height dimension.
    x (int): Width dimension.

    Returns:
    ndarray: 3D array of distances between changes in sign.
    """
    switch_distance_stack = np.zeros((t, y, x), dtype=int)

    for i in range(y):
        for j in range(x):
            series = arr[:, i, j]
            changes = []
            start_found = False
            for k in range(t):
                if not start_found and series[k] != 0:
                    start_found = True
                    changes.append(k)
                elif start_found and (series[k] == 1 or series[k] == -1):
                    if series[k] != series[changes[-1]]:
                        changes.append(k)

            # Calculate distances
            if len(changes) > 1:
                for c in range(1, len(changes)):
                    distance = changes[c] - changes[c - 1]
                    switch_distance_stack[c - 1, i, j] = distance

    return switch_distance_stack


# =============================================================================
# MAIN SCRIPT
# =============================================================================
if __name__ == "__main__":
    
    run_names = ['q05_1', 'q07_1', 'q10_2', 'q15_2', 'q15_3', 'q20_2']
    
    plot_dir_path = os.path.join(os.getcwd(), 'plots', 'DoD_switch_analysis')
    if not os.path.exists(plot_dir_path):
        os.makedirs(plot_dir_path)
    
    # Common figure size for histograms
    fig_size = (10, 6)
    
    for run_name in run_names:
        # Load the data
        data = np.load(os.path.join(os.getcwd(), 'output/DoDs/DoDs_stack','DoD_stack_bool_' + run_name + '.npy' ))
        data = data[:,:,:,0]  # Select the first channel of the data
        t, y, x = data.shape[0], data.shape[1], data.shape[2]
        
        # Task 1: Calculate consecutive periods for 1, 0, and -1
        ones_stack, zeros_stack, neg_ones_stack = calculate_consecutive_periods(data, t, y, x)
    
        # Task 2: Calculate distance between changes in sign
        switch_distance_stack = calculate_distance_between_switches(data, t, y, x)
    
        # Task 3: Plot distribution of consecutive periods (-1, 0, 1)
        # Flatten stacks to get all period lengths
        ones_lengths = ones_stack.flatten()
        zeros_lengths = zeros_stack.flatten()
        neg_ones_lengths = neg_ones_stack.flatten()
    
        # Compute histograms using np.histogram (without density)
        ones_hist, ones_bins = np.histogram(ones_lengths[ones_lengths > 0], bins=range(1, t + 2))
        zeros_hist, zeros_bins = np.histogram(zeros_lengths[zeros_lengths > 0], bins=range(1, t + 2))
        neg_ones_hist, neg_ones_bins = np.histogram(neg_ones_lengths[neg_ones_lengths > 0], bins=range(1, t + 2))
    
        # Normalize the histograms by dividing by the total number of elements
        ones_hist_normalized = ones_hist / np.sum(ones_hist)
        zeros_hist_normalized = zeros_hist / np.sum(zeros_hist)
        neg_ones_hist_normalized = neg_ones_hist / np.sum(neg_ones_hist)
    
        # Plot for 1s
        plt.figure(figsize=fig_size)
        plt.bar(ones_bins[:-1], ones_hist_normalized, width=1, alpha=0.7, color='green', label='Consecutive 1s')
        plt.xlabel('Length of consecutive 1s')
        plt.ylabel('Normalized frequency')
        plt.title(run_name + ' Normalized distribution of consecutive 1s')
        plt.legend()
        plt.grid(True, axis='y')
        plt.ylim(0, 1)
        plt.xticks(range(1, t + 1))  # Set x-axis ticks at each unit
        # Add text to the plot with the script name
        plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
                 transform=plt.gca().transAxes, ha='center', va='center', fontsize=12)
        plt.savefig(os.path.join(plot_dir_path, run_name + '_consecutive_1s.pdf'))
    
        # Plot for 0s
        plt.figure(figsize=fig_size)
        plt.bar(zeros_bins[:-1], zeros_hist_normalized, width=1, alpha=0.7, color='blue', label='Consecutive 0s')
        plt.xlabel('Length of consecutive 0s')
        plt.ylabel('Normalized frequency')
        plt.title(run_name + ' Normalized distribution of consecutive 0s')
        plt.legend()
        plt.grid(True, axis='y')
        plt.ylim(0, 0.6)
        plt.xticks(range(1, t + 1))  # Set x-axis ticks at each unit
        # Add text to the plot with the script name
        plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
                 transform=plt.gca().transAxes, ha='center', va='center', fontsize=12)
        plt.savefig(os.path.join(plot_dir_path, run_name + '_consecutive_0s.pdf'))
    
        # Plot for -1s
        plt.figure(figsize=fig_size)
        plt.bar(neg_ones_bins[:-1], neg_ones_hist_normalized, width=1, alpha=0.7, color='red', label='Consecutive -1s')
        plt.xlabel('Length of consecutive -1s')
        plt.ylabel('Normalized frequency')
        plt.title(run_name + ' Normalized distribution of consecutive -1s')
        plt.legend()
        plt.grid(True, axis='y')
        plt.ylim(0, 1)
        plt.xticks(range(1, t + 1))  # Set x-axis ticks at each unit
        # Add text to the plot with the script name
        plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
                 transform=plt.gca().transAxes, ha='center', va='center', fontsize=12)
        plt.savefig(os.path.join(plot_dir_path, run_name + '_consecutive_neg_1s.pdf'))
    
        # Task 4: Plot distribution of switch distances
        switch_distances = switch_distance_stack.flatten()
        switch_hist, switch_bins = np.histogram(switch_distances[switch_distances > 0], bins=range(1, t + 2))
        switch_hist_normalized = switch_hist / np.sum(switch_hist)
    
        plt.figure(figsize=fig_size)
        plt.bar(switch_bins[:-1], switch_hist_normalized, width=1, alpha=0.7, color='purple', label='Switch distances')
        plt.xlabel('Distance between switches')
        plt.ylabel('Normalized frequency')
        plt.title(run_name + ' Normalized distribution of distances between switches')
        plt.legend()
        plt.grid(True, axis='y')
        plt.ylim(0, 0.5)
        plt.xticks(range(1, t + 1))  # Set x-axis ticks at each unit
        # Add text to the plot with the script name
        plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
                 transform=plt.gca().transAxes, ha='center', va='center', fontsize=12)
        plt.savefig(os.path.join(plot_dir_path, run_name + '_switch_distances.pdf'))
    
        # Task 5: Plot a map of the number of switches in each (y,x) position
        switch_count_map = np.count_nonzero(switch_distance_stack, axis=0)
    
        plt.figure(figsize=(14, 8))
        plt.imshow(switch_count_map, cmap='viridis', interpolation='nearest')
        plt.colorbar(label='Number of switches')
        plt.title(run_name + ' Map of number of switches per cell (y,x)')
        # Add text to the plot with the script name
        plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
                 transform=plt.gca().transAxes, ha='center', va='center', fontsize=12)
        plt.savefig(os.path.join(plot_dir_path, run_name + '_switch_map.pdf'))
        
        # Task 6: Plot histogram of number of switches per (y,x) cell
        switch_counts = switch_count_map.flatten()
        switch_count_hist, switch_count_bins = np.histogram(switch_counts[switch_counts > 0], bins=range(1, np.max(switch_counts) + 2))
        switch_count_hist_normalized = switch_count_hist / np.sum(switch_count_hist)

        plt.figure(figsize=fig_size)
        plt.bar(switch_count_bins[:-1], switch_count_hist_normalized, width=1, alpha=0.7, color='orange', label='Number of switches per cell')
        plt.xlabel('Number of switches per cell')
        plt.ylabel('Normalized frequency')
        plt.title(run_name + ' Normalized distribution of number of switches per cell')
        plt.legend()
        plt.grid(True, axis='y')
        plt.ylim(0, 1)
        plt.xlim(0, 8)
        plt.xticks(range(1, np.max(switch_counts) + 1))  # Set x-axis ticks for number of switches
        # Add text to the plot with the script name
        plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}",
                 transform=plt.gca().transAxes, ha='center', va='center', fontsize=12)
        plt.savefig(os.path.join(plot_dir_path, run_name + '_switch_count_histogram.pdf'))
        
        # Show all plots
        plt.show()
