#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 19:37:56 2024

@author: erri
"""

import numpy as np
import os
import matplotlib.pyplot as plt

# =============================================================================
# DEFINE FOLDERS PATH
# =============================================================================
plot_dir = os.path.join(os.getcwd(), 'plots', 'DEM_value_destination')
if not(os.path.exists(plot_dir)):
        os.mkdir(plot_dir)


#%%============================================================================
# DEM VALUE DESTINATION FOR DIFFERENT JUMP LENGTH
# =============================================================================
'''
The aim of this section 
'''
# -----------------------------------------------------------------------------
# RUN SELECTION
set_names = ['q05_1','q07_1', 'q10_2', 'q15_3', 'q20_2']


# Plot mode
plot_mode = [
              'overall_class_distribution_bars'
             ]
analysis_mode = [
    'computed_DEM'  # Use the DoD to compute DEMs after the first one
    ]
# ANALYSIS PARAMETERS
jumps = np.array([1,2,3,4,5,6,7,8])
# jumps = np.array([1])
# jumps = np.array([2])

# Initialize array
bar_top_persistence = np.zeros((len(set_names),len(jumps)))
pit_persistence = np.zeros((len(set_names),len(jumps)))


for set_name in set_names:
    print(set_name, ' is running...')
    
    home_dir = os.getcwd()
    output_dir = os.path.join(home_dir, 'output', set_name, 'DEM_analysis')
    if not(os.path.exists(output_dir)):
        os.mkdir(output_dir)
        
    # output_dir = os.path.join(home_dir, 'output', set_name, 'DEM_analysis', 'distribution_analysis')
    # if not(os.path.exists(output_dir)):
    #     os.mkdir(output_dir)
    

    DEM_path = os.path.join(home_dir, 'surveys',set_name)
    stack_path = os.path.join(DEM_path, set_name + '_DEM_stack.npy')
    DEM_data_raw = np.load(stack_path)
    DEM_data_raw = np.where(DEM_data_raw==-999, np.nan, DEM_data_raw)
    dim_t, dim_y, dim_x = DEM_data_raw.shape
    
    DEM_data = np.copy(DEM_data_raw)
    
    # IMPORT DoD DATA:
    DoD_stack_path = os.path.join(home_dir,'output','DoDs', 'DoDs_stack', 'DoD_stack_'+set_name+'.npy')
    DoD_raw = np.load(DoD_stack_path)
    DoD_data = DoD_raw[:,:,:,0]
    
    # COMPUTE THE DEM2 AS DEM1 + DoD2-1
    # DoD_data = np.where(DoD_data==0, np.nan, DoD_data)
    computed_DEM_data = DEM_data[:-1,:,:] + DoD_data[:,:,:]
    
    for jump in jumps:
        
        
        # DISTRIBUTION PARAMETERS
        num_bins = 7
        # DEM_bin_edges = [-50.0,-11.52,-4.59,-0.58,2.345,4.535,7.09,20.0]
        # DEM_bin_edges = [-60.0,-12.62,-5.69,-1.67,1.25,3.805,6.36,60.0]
        
        # BIN EDGES UPLOADED FROM ./DoD_analysis/array_values_equal_class_subdidision.py
        DEM_bin_edges = np.loadtxt(os.path.join(home_dir, 'output', 'distribution_edges', set_name + '_DEM_distribution_edges.txt'), delimiter=',')
        
        
        DEM_class_labels = [f"{DEM_bin_edges[i]},{DEM_bin_edges[i+1]}" for i in range(len(DEM_bin_edges) - 1)]
        hist_range = (-50, 50)  # Range of values to include in the histogram
        
        print()
        print()
        print('Jump: ', jump)
        print()
        
        
        for t in range(0,DEM_data.shape[0]-1-jump):
        # for t in range(0,1):
            print('t: ', t)
            print('t-1+jump: ', t-1+jump)
        
            # This is the array where the overall distribution will be stored
            dist_array_class1 = []
            dist_array_class2 = []
            dist_array_class3 = []
            dist_array_class4 = []
            dist_array_class5 = []
            dist_array_class6 = []
            dist_array_class7 = []
            
            for x in range(DEM_data.shape[1]):
                for y in range(DEM_data.shape[2]):
                    if 'computed_DEM' in analysis_mode:
                        if DEM_bin_edges[0]<=DEM_data[t,x,y]<DEM_bin_edges[1]:
                            dist_array_class1 = np.append(dist_array_class1, computed_DEM_data[t-1+jump,x,y])
                        if DEM_bin_edges[1]<=DEM_data[t,x,y]<DEM_bin_edges[2]:
                            dist_array_class2 = np.append(dist_array_class2, computed_DEM_data[t-1+jump,x,y])
                        if DEM_bin_edges[2]<=DEM_data[t,x,y]<DEM_bin_edges[3]:
                            dist_array_class3 = np.append(dist_array_class3, computed_DEM_data[t-1+jump,x,y])
                        if DEM_bin_edges[3]<=DEM_data[t,x,y]<DEM_bin_edges[4]:
                            dist_array_class4 = np.append(dist_array_class4, computed_DEM_data[t-1+jump,x,y])
                        if DEM_bin_edges[4]<=DEM_data[t,x,y]<DEM_bin_edges[5]:
                            dist_array_class5 = np.append(dist_array_class5, computed_DEM_data[t-1+jump,x,y])
                        if DEM_bin_edges[5]<DEM_data[t,x,y]<DEM_bin_edges[6]:
                            dist_array_class6 = np.append(dist_array_class6, computed_DEM_data[t-1+jump,x,y])
                        if DEM_bin_edges[6]<=DEM_data[t,x,y]<DEM_bin_edges[7]:
                            dist_array_class7 = np.append(dist_array_class7, computed_DEM_data[t-1+jump,x,y])
                    else:
                        if DEM_bin_edges[0]<=DEM_data[t,x,y]<DEM_bin_edges[1]:
                            dist_array_class1 = np.append(dist_array_class1, DEM_data[t+jump,x,y])
                        if DEM_bin_edges[1]<=DEM_data[t,x,y]<DEM_bin_edges[2]:
                            dist_array_class2 = np.append(dist_array_class2, DEM_data[t+jump,x,y])
                        if DEM_bin_edges[2]<=DEM_data[t,x,y]<DEM_bin_edges[3]:
                            dist_array_class3 = np.append(dist_array_class3, DEM_data[t+jump,x,y])
                        if DEM_bin_edges[3]<=DEM_data[t,x,y]<DEM_bin_edges[4]:
                            dist_array_class4 = np.append(dist_array_class4, DEM_data[t+jump,x,y])
                        if DEM_bin_edges[4]<=DEM_data[t,x,y]<DEM_bin_edges[5]:
                            dist_array_class5 = np.append(dist_array_class5, DEM_data[t+jump,x,y])
                        if DEM_bin_edges[5]<DEM_data[t,x,y]<DEM_bin_edges[6]:
                            dist_array_class6 = np.append(dist_array_class6, DEM_data[t+jump,x,y])
                        if DEM_bin_edges[6]<=DEM_data[t,x,y]<DEM_bin_edges[7]:
                            dist_array_class7 = np.append(dist_array_class7, DEM_data[t+jump,x,y])
        
        
        # OVERALL DISTRIBUTION
        overall_dist_matrix = np.zeros((num_bins,num_bins))
        
        hist_array_1 = np.copy(dist_array_class1)
        hist_array_1 = hist_array_1[~np.isnan(hist_array_1)] # Trim np.nan
        hist_1, DEM_bin_edges1 = np.histogram(hist_array_1, bins=DEM_bin_edges, range=hist_range)
        hist_1 = hist_1/np.nansum(hist_1)
        np.savetxt(os.path.join(output_dir, 'jump' + str(jump) + '_class1_distribution_DEM2DEM.txt'), np.round(hist_1, decimals=2))
        overall_dist_matrix[0,:] = hist_1
        
        hist_array_2 = np.copy(dist_array_class2)
        hist_array_2 = hist_array_2[~np.isnan(hist_array_2)] # Trimnp.nan
        hist_2, DEM_bin_edges2 = np.histogram(hist_array_2, bins=DEM_bin_edges, range=hist_range)
        hist_2 = hist_2/np.nansum(hist_2)
        np.savetxt(os.path.join(output_dir, 'jump' + str(jump) + '_class2_distribution_DEM2DEM.txt'), np.round(hist_2, decimals=2))
        overall_dist_matrix[1,:] = hist_2
        
        hist_array_3 = np.copy(dist_array_class3)
        hist_array_3 = hist_array_3[~np.isnan(hist_array_3)] # Trim np.nan
        hist_3, DEM_bin_edges3 = np.histogram(hist_array_3, bins=DEM_bin_edges, range=hist_range)
        hist_3 = hist_3/np.nansum(hist_3)
        np.savetxt(os.path.join(output_dir, 'jump' + str(jump) + '_class3_distribution_DEM2DEM.txt'), np.round(hist_3, decimals=2))
        overall_dist_matrix[2,:] = hist_3
        
        hist_array_4 = np.copy(dist_array_class4)
        hist_array_4 = hist_array_4[~np.isnan(hist_array_4)] # Trim np.nan
        hist_4, DEM_bin_edges4 = np.histogram(hist_array_4, bins=DEM_bin_edges, range=hist_range)
        hist_4 = hist_4/np.nansum(hist_4)
        np.savetxt(os.path.join(output_dir, 'jump' + str(jump) + '_class4_distribution_DEM2DEM.txt'), np.round(hist_4, decimals=2))
        overall_dist_matrix[3,:] = hist_4
        
        hist_array_5 = np.copy(dist_array_class5)
        hist_array_5 = hist_array_5[~np.isnan(hist_array_5)] # Trim np.nan
        hist_5, DEM_bin_edge5 = np.histogram(hist_array_5, bins=DEM_bin_edges, range=hist_range)
        hist_5 = hist_5/np.nansum(hist_5)
        np.savetxt(os.path.join(output_dir, 'jump' + str(jump) + '_class5_distribution_DEM2DEM.txt'), np.round(hist_5, decimals=2))
        overall_dist_matrix[4,:] = hist_5
        
        hist_array_6 = np.copy(dist_array_class6)
        hist_array_6 = hist_array_6[~np.isnan(hist_array_6)] # Trim np.nan
        hist_6, DEM_bin_edges6 = np.histogram(hist_array_6, bins=DEM_bin_edges, range=hist_range)
        hist_6 = hist_6/np.nansum(hist_6)
        np.savetxt(os.path.join(output_dir, 'jump' + str(jump) + '_class7_distribution_DEM2DEM.txt'), np.round(hist_6, decimals=2))
        overall_dist_matrix[5,:] = hist_6
        
        hist_array_7 = np.copy(dist_array_class7)
        hist_array_7 = hist_array_7[~np.isnan(hist_array_7)] # Trim np.nan
        hist_7, DEM_bin_edges7 = np.histogram(hist_array_7, bins=DEM_bin_edges, range=hist_range)
        hist_7 = hist_7/np.nansum(hist_7)
        np.savetxt(os.path.join(output_dir, 'jump' + str(jump) + '_class7_distribution_DEM2DEM.txt'), np.round(hist_7, decimals=2))
        overall_dist_matrix[6,:] = hist_7
        
        # Fill matrices:
        # hill_persistence[set_names.index(set_name),jump-1] = overall_dist_matrix[-1,-2] # Class 1
        # pit_persistence[set_names.index(set_name),jump-1] = overall_dist_matrix[0,1] # Class 7
        
        bar_top_persistence[set_names.index(set_name),jump-1] = overall_dist_matrix[-1,-2] # Class 2
        pit_persistence[set_names.index(set_name),jump-1] = overall_dist_matrix[0,1] # Class 6
        
        
        np.savetxt(os.path.join(output_dir, set_name + '_jump' + str(jump)+'_overall_distribution_matrix_DEM2DEM.txt'), np.round(overall_dist_matrix, decimals=6))

            #%%        
        if 'overall_class_distribution_bars' in plot_mode:
            
            np.loadtxt(os.path.join(output_dir, set_name + '_jump' + str(jump)+'_overall_distribution_matrix_DEM2DEM.txt'))
            # Your data
            data1 = np.copy(overall_dist_matrix[:,0])
            data2 = np.copy(overall_dist_matrix[:,1])
            data3 = np.copy(overall_dist_matrix[:,2])
            data4 = np.copy(overall_dist_matrix[:,3])
            data5 = np.copy(overall_dist_matrix[:,4])
            data6 = np.copy(overall_dist_matrix[:,5])
            data7 = np.copy(overall_dist_matrix[:,6])
            
            # Bar width and alpha
            bar_width = 0.6
            alpha = 0.6
            
            classes = DEM_class_labels
            class_distance = 1
            index = np.arange(len(classes)) * (bar_width * 7 + class_distance)
            
            colors = ['#440154', '#443983', '#31688e', '#21918c', '#35b779', '#90d743', '#fde725']
            
            plt.figure(figsize=(10, 6))
            
            # Create the bar chart
            plt.bar(index - 3*bar_width, data1, width=bar_width, label=DEM_class_labels[0], color=colors[0])
            plt.bar(index - 2*bar_width, data2, width=bar_width, label=DEM_class_labels[1], color=colors[1])
            plt.bar(index - bar_width, data3, width=bar_width, label=DEM_class_labels[2], color=colors[2])
            plt.bar(index, data4, width=bar_width, label=DEM_class_labels[3], color=colors[3])
            plt.bar(index + bar_width, data5, width=bar_width, label=DEM_class_labels[4], color=colors[4])
            plt.bar(index + 2*bar_width, data6, width=bar_width, label=DEM_class_labels[5], color=colors[5])
            plt.bar(index + 3*bar_width, data7, width=bar_width, label=DEM_class_labels[6], color=colors[6])
            
            # Configure the chart
            plt.xlabel('Origin DEMs value classes')
            plt.ylabel('Frequency ')
            plt.title(set_name + ' - Jump: ' + str(jump))
            # Set x-axis labels
            plt.xticks(index, DEM_class_labels, rotation=15)  # Rotate x-axis labels for readability
            
            # Add a legend
            plt.legend(bbox_to_anchor=(1.001, 1), loc='upper left')
            
            # Set the y limit
            # plt.ylim(0,0.50)
            plt.ylim(0,1.0)
            
            # Add text to the plot with the script name
            plt.text(0.5, -0.2, f"Generated by: {os.path.basename(__file__)}", transform=plt.gca().transAxes, ha='center', va='center', fontsize=12)
            
            # Save image and report
            plot_path = os.path.join(plot_dir, set_name + '_jump' + str(jump) + '_DEM2DEM_overall_distribution_sep.pdf')
            plt.savefig(plot_path, dpi=400, bbox_inches='tight')
            
            # Show the plot
            plt.show()
    
    header = 'rows are runs, columns at different jumps \n' + str(set_names)

    np.savetxt(os.path.join(home_dir, 'output', 'DEM2DEM_bar_top_persitence.txt'), bar_top_persistence, header=header, delimiter='\t')
    np.savetxt(os.path.join(home_dir, 'output', 'DEM2DEM_pit_persitence.txt'), pit_persistence, header=header, delimiter='\t')

#%%############################################################################
# ALGOROTHM CHECK
###############################################################################
set_name = 'q07_1'
home_dir = os.getcwd()
output_dir = os.path.join(home_dir, 'output', set_name, 'DEM_analysis')

# IMPORT DATA
DEM_path = os.path.join(home_dir, 'surveys',set_name)
stack_path = os.path.join(DEM_path, set_name + '_DEM_stack.npy')
DEM_data_raw = np.load(stack_path)
DEM_data_raw = np.where(DEM_data_raw==-999, np.nan, DEM_data_raw)
dim_t, dim_y, dim_x = DEM_data_raw.shape

DEM_data = np.copy(DEM_data_raw)

# IMPORT DoD DATA:
DoD_stack_path = os.path.join(home_dir,'output','DoDs', 'DoDs_stack', 'DoD_stack_'+set_name+'.npy')
DoD_raw = np.load(DoD_stack_path)
DoD_data = DoD_raw[:,:,:,0]

# COMPUTE THE DEM2 AS DEM1 + DoD2-1
computed_DEM_data = DEM_data[:-1,:,:] + DoD_data[:,:,:]
# computed_DEM_data = DEM_data

# IMPORT DISTRIBUTION
# overall_dist_matrix = np.loadtxt(os.path.join(output_dir, set_name + '_overall_distribution_matrix_DEM2DEM.txt'))


#%% COMPUTING...

# DEM0 = DEM_data[0,:,:] # DEM0
# # DEM1 = DEM_data[1,:,:] # DEM1
# DoD1_0 = DoD_data[0,:,:] # DoD 1-0
# DEM1 = DEM0 + DoD1_0 # Computed DEM1

# ## DEM2DEM Extract values that are in class7 both for DEM0 and DEM1
# # DEM0_c7 = np.where(DEM0>6.36,1,np.nan)
# # DEM1_c7 = np.where(DEM1>6.36,1,np.nan)

# DEM0_c7 = np.where(DEM0>6.36,1,0)
# DEM1_c7 = np.where(DEM1>6.36,1,0)

# DEM0_c7_to_DEM1_c7 = DEM0_c7*DEM1_c7 # DEM_c7 persistence
# percent_persistence = np.nansum(DEM0_c7_to_DEM1_c7)/np.nansum(DEM0_c7)
# print('\n \n')
# print('Number of points in DEM0 class7: ', np.nansum(DEM0_c7))
# print('Number of points in DEM1 class7: ', np.nansum(DEM1_c7))
# print('Number of points that both was in DEM0 class 7 and are in DEM1 class7: ', np.nansum(DEM0_c7_to_DEM1_c7))

# print('Percentage of persistence DEM classe 7: ', percent_persistence)

# print('Percentage of persistence DEM classe 7 from histogram: ', overall_dist_matrix[6,6])

