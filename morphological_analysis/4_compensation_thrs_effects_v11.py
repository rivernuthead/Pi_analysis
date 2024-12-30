import os
import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.cm as cm

start = time.time() # Set initial time

# =============================================================================
# SCRIPT MODE
# =============================================================================
'''
set_name mode:
    0 = set_names in the set_names list
    1 = one set_name at time
    2 = bath process 'all the set_names in the folder'
plot_mode:
    ==0 plot mode OFF
    ==1 plot mode ON
'''
set_name_mode = 1
plot_mode = 1

# =============================================================================
# ARRAY OF set_names
# =============================================================================
set_names = ['q07_1', 'q10_2', 'q15_3', 'q20_2']
# set_names = ['q07_1', 'q10_2', 'q10_3', 'q10_4','q15_2', 'q15_3', 'q20_2']
# set_names = ['q10_3', 'q10_4']

set_names = ['q05_1']
# set_names = ['q07_1']
# set_names = ['q10_2']
# set_names = ['q15_2']
# set_names = ['q15_3']
# set_names = ['q20_2']

# home_dir = os.path.join(os.getcwd(), 'morphological_analysis') # Home directory
home_dir = os.path.join(os.getcwd(), 'morphological_analysis') # Home directory
output_dir = os.path.join(home_dir, 'output_data','comp_thrs_analysis')
if not(os.path.exists(output_dir)):
        os.mkdir(output_dir)
report_dir = os.path.join(home_dir, 'output_data', 'report')

# INITIALIZE THE FIGURE -------------------------------------------------------
for set_name in set_names:
    # FOLDER SETUP ------------------------------------------------------------
    
    report_dir = os.path.join(report_dir, 'report_' + set_name)
    if not(os.path.exists(report_dir)):
        os.mkdir(report_dir)
    set_name_dir = os.path.join(home_dir, 'surveys')
    DoDs_folder = os.path.join(home_dir, 'output_data', 'DoDs', 'DoDs_stack') # Input folder

    print(set_name, ' is running...')
    print()
    
    # SET EXNER TIME ----------------------------------------------------------
    if '05' in set_name:
        Txnr = 100
        Txnr_scale = 0.81
    elif '07' in set_name:
        Txnr_scale = 0.64
        Txnr = 74
    elif '10' in set_name:
        Txnr_scale = 0.58
        Txnr = 48
    elif '15' in set_name:
        Txnr_scale = 0.70
        Txnr = 30
    elif '20' in set_name:
        Txnr_scale = 0.80
        Txnr = 20

    # IMPORT STACK ------------------------------------------------------------
    stack_dir = os.path.join(home_dir, 'output_data', 'DoDs', 'DoDs_stack') # Define the stack directory
    stack=np.load(os.path.join(stack_dir, 'DoD_stack_'+set_name+'.npy')) # Load the stack
    stack_bool=np.load(os.path.join(stack_dir, 'DoD_stack_bool_'+set_name+'.npy'))
    
    # ADD MAW CURVE -----------------------------------------------------------
    # Compute MAW
    MAW_stack    = np.nansum(abs(stack_bool), axis=1)/120
    MAW_matrix   = np.nanmean(MAW_stack, axis=1) # MAW matrix
    MAW_matrix   = np.where(MAW_matrix==0, np.nan, MAW_matrix)
    MAW_mean_std = np.stack((np.nanmean(MAW_matrix, axis=0), np.nanstd(MAW_matrix, axis=0)))
    
    # COMPUTE THE envMAW-------------------------------------------------------
    envMAW_stack      = np.copy(abs(stack_bool))
    envMAW_stack      = envMAW_stack[:,:,:,0]
    envMAW_stack      = np.cumsum(envMAW_stack, axis=0)
    envMAW_bool_stack = np.where(envMAW_stack>0,1,envMAW_stack)
    envMAW_array      = np.nansum(envMAW_bool_stack, axis=1)/120
    envMAW            = np.nanmean(envMAW_array, axis=1)
    
    
    # Define the report stack in which all the difference matrix will be stored
    report_stack            = np.zeros((stack_bool.shape[0]-1,stack_bool.shape[1],stack_bool.shape[2],stack_bool.shape[3]-1))
    report_stack            = np.where(report_stack==0, np.nan, np.nan)
    compensation_stack      = np.copy(report_stack)
    compensation_pure_stack = np.copy(report_stack)
    threshold_pure_stack    = np.copy(report_stack)
    
    index=0
    for i in range(0,stack.shape[3]-1): # For each delta
        
        for t in range(0,stack.shape[0]-1-index):
            
            # Extract the MAA matrices
            envMAA = stack_bool[t:t+2+i,:,:,0] # envMAA
            
            # Compute the DoD as boolean map
            DoD = abs(stack_bool[t,:,:,i+1])
            
            for x in range(0,stack_bool.shape[2]):
                for y in range(0,stack_bool.shape[1]):
                    
                    envMA_array = envMAA[:,y,x]
                    DoD_array   = DoD[y,x]
                    
                    compensation      = np.any(abs(envMA_array) == 1) & (DoD_array == 0)
                    compensation_pure = np.any(envMA_array == -1) and np.any(envMA_array == 1) and np.all(DoD_array == 0)
                    threshold_pure    = np.all((envMA_array == 0) & (DoD_array == 1))
                    
                    compensation_stack[t,y,x,i]      = compensation*1
                    compensation_pure_stack[t,y,x,i] = compensation_pure*1
                    threshold_pure_stack[t,y,x,i]    = threshold_pure*1
        

# =============================================================================
#   SAVE MAPS STACK AS NPY FILE
# =============================================================================
    np.save(os.path.join(stack_dir, set_name + '_compensation_map_stack.npy'), compensation_stack)
    np.save(os.path.join(stack_dir, set_name + '_compensation_pure_map_stack.npy'), compensation_pure_stack)   
    np.save(os.path.join(stack_dir, set_name + '_threshold_pure_map_stack.npy'), threshold_pure_stack)

# =============================================================================
# SAVE MAPS AS SINGLE MAP     
# =============================================================================
    for i in range(0,compensation_stack.shape[0]):
            comp_map      = np.save(os.path.join(stack_dir, set_name + '_' + str(i+1) + '_compensation_map.npy'), compensation_stack[i,:,:,0])
            pure_comp_map = np.save(os.path.join(stack_dir, set_name + '_' + str(i+1) + '_compensation_pure_map.npy'), compensation_pure_stack[i,:,:,0])
            thrs_pure_map = np.save(os.path.join(stack_dir, set_name + '_' + str(i+1) + '_threshold_pure_map.npy'), threshold_pure_stack[i,:,:,0])


# =============================================================================
#   BUILD REPORT AND MEAN - STDEV ARRAY FOR PLOTS
# =============================================================================
    
    report_compensation_stack = np.nansum(compensation_stack, axis=2)
    report_compensation_stack = np.nansum(report_compensation_stack, axis=1)
    report_compensation_stack = report_compensation_stack/stack_bool.shape[2]/120
    report_compensation_stack = np.where(report_compensation_stack==0, np.nan, report_compensation_stack)
    compensation_mean_array   = np.nanmean(report_compensation_stack, axis=0)
    compensation_std_array    = np.nanstd(report_compensation_stack, axis=0)
    
    report_compensation_pure_stack = np.nansum(compensation_pure_stack, axis=2)
    report_compensation_pure_stack = np.nansum(report_compensation_pure_stack, axis=1)
    report_compensation_pure_stack = report_compensation_pure_stack/stack_bool.shape[2]/120
    report_compensation_pure_stack = np.where(report_compensation_pure_stack==0, np.nan, report_compensation_pure_stack)
    compensation_pure_mean_array   = np.nanmean(report_compensation_pure_stack, axis=0)
    compensation_pure_std_array    = np.nanstd(report_compensation_pure_stack, axis=0)
    
    report_threshold_pure_stack = np.nansum(threshold_pure_stack, axis=2)
    report_threshold_pure_stack = np.nansum(report_threshold_pure_stack, axis=1)
    report_threshold_pure_stack = report_threshold_pure_stack/stack_bool.shape[2]/120
    report_threshold_pure_stack = np.where(report_threshold_pure_stack==0, np.nan, report_threshold_pure_stack)
    threshold_pure_mean_array   = np.nanmean(report_threshold_pure_stack, axis=0)
    threshold_pure_std_array    = np.nanstd(report_threshold_pure_stack, axis=0)
    
    np.savetxt(os.path.join(output_dir, set_name + '_compensation_dimless_width.txt'), report_compensation_pure_stack, delimiter=',', fmt = "%.4f")
    np.savetxt(os.path.join(output_dir, set_name + '_under_thrs_dimless_width.txt'), report_threshold_pure_stack, delimiter=',', fmt = "%.4f")
    
    # Save mean and std arrays for plotting
    np.savetxt(os.path.join(output_dir, set_name + '_compensation_mean_array.txt'), compensation_mean_array, delimiter=',', fmt = "%.4f")
    np.savetxt(os.path.join(output_dir, set_name + '_compensation_std_array.txt'), compensation_std_array, delimiter=',', fmt = "%.4f")
    np.savetxt(os.path.join(output_dir, set_name + '_compensation_pure_mean_array.txt'), compensation_pure_mean_array, delimiter=',', fmt = "%.4f")
    np.savetxt(os.path.join(output_dir, set_name + '_compensation_pure_std_array.txt'), compensation_pure_std_array, delimiter=',', fmt = "%.4f")
    np.savetxt(os.path.join(output_dir, set_name + '_threshold_pure_mean_array.txt'), threshold_pure_mean_array, delimiter=',', fmt = "%.4f")
    np.savetxt(os.path.join(output_dir, set_name + '_threshold_pure_std_array.txt'), threshold_pure_std_array, delimiter=',', fmt = "%.4f")
    np.savetxt(os.path.join(output_dir, set_name + '_MAW_mean_std.txt'), MAW_mean_std, delimiter=',', fmt = "%.4f")
    # np.savetxt(os.path.join(output_dir, set_name + '_envMAW.txt'), envMAW, delimiter=',', fmt = "%.4f")
    

print('End of data saving phase')

# #%%============================================================================
# # PLOTS
# # =============================================================================

# set_names = ['q05_1', 'q07_1', 'q10_2', 'q15_3', 'q20_2']

# set_names = ['q05_1']
# # set_names = ['q07_1']
# # set_names = ['q10_2']



# # Function to extract the last 80% of a colormap
# def extract_colormap(cmap_name, num_colors):
#     cmap = cm.get_cmap(cmap_name)
#     colors = cmap(np.linspace(0.4, 1, num_colors))
#     return colors

# # Generate colormaps for each type of plot
# thrs_colors = extract_colormap('Purples', len(set_names))
# comp_colors = extract_colormap('Greens', len(set_names))
# maw_colors = extract_colormap('Greys', len(set_names))
# envmaw_colors = extract_colormap('Blues', len(set_names))
# baw_colors = extract_colormap('Reds', len(set_names))


# # =============================================================================
# # # PLOT FOR EACH SET_NAME
# # =============================================================================
# for i, set_name in enumerate(set_names):
#     plt.figure(figsize=(10, 8))
    
#     if set_name == 'q15_3':
#         BAW_set_name = 'q15_2'
#     else:
#         BAW_set_name = set_name
        
#     # SET EXNER TIME ----------------------------------------------------------
#     if '05' in set_name:
#         Txnr = 100
#         Txnr_scale = 0.81
#     elif '07' in set_name:
#         Txnr_scale = 0.64
#         Txnr = 74
#     elif '10' in set_name:
#         Txnr_scale = 0.58
#         Txnr = 48
#     elif '15' in set_name:
#         Txnr_scale = 0.70
#         Txnr = 30
#     elif '20' in set_name:
#         Txnr_scale = 0.80
#         Txnr = 20
    
#     # IMPORT STACK ------------------------------------------------------------
#     stack_dir = os.path.join(home_dir, 'output', 'DoDs', 'DoDs_stack') # Define the stack directory
#     stack=np.load(os.path.join(stack_dir, 'DoD_stack_'+set_name+'.npy')) # Load the stack
#     stack_bool=np.load(os.path.join(stack_dir, 'DoD_stack_bool_'+set_name+'.npy'))

#     # COMPUTE X-VALUES ARRAYS------------------------------------------------------
#     # Compensation and threshold
#     comp_thrs_x_values = np.linspace(2, stack.shape[0], stack.shape[0]-1)/2 * Txnr_scale

#     # MAW
#     MAW_x_values = np.linspace(1, stack.shape[0], stack.shape[0])/2*Txnr_scale
    
#     # Load the mean and std arrays from text files
#     output_dir = os.path.join(home_dir, 'output', 'report_' + set_name, 'comp_thrs_analysis')
#     compensation_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_mean_array.txt'), delimiter=',')
#     compensation_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_std_array.txt'), delimiter=',')
#     compensation_pure_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_pure_mean_array.txt'), delimiter=',')
#     compensation_pure_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_pure_std_array.txt'), delimiter=',')
#     threshold_pure_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_threshold_pure_mean_array.txt'), delimiter=',')
#     threshold_pure_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_threshold_pure_std_array.txt'), delimiter=',')
#     MAW_mean_std = np.loadtxt(os.path.join(output_dir, set_name + '_MAW_mean_std.txt'), delimiter=',')

#     envBAW_mean_std = np.loadtxt(os.path.join('/home/erri/Documents/PhD/Research/5_research_repos/PiQs_analysis/output_report/', BAW_set_name + '_envBAW_merged_single_runs.txt'), delimiter=',', skiprows=0)
#     envMAW = np.loadtxt(os.path.join(home_dir,'output','report_'+ set_name, 'envelopes', set_name + '_timespan0_envMAW_report.txt'), delimiter=',')
    
#     # BAW_x_values = envBAW_mean_std[:,0]*np.max(MAW_x_values)/np.max(envBAW_mean_std[:,0])
    
#     # Plot the datasets
#     plt.errorbar(envBAW_mean_std[0,:]/Txnr, envBAW_mean_std[1,:], yerr=envBAW_mean_std[1,:], fmt='s', linestyle='--', ecolor='black', color='red', capsize=5, label='envBAW Mean ± Std')
#     plt.errorbar(MAW_x_values, envMAW[:,1], yerr=envMAW[:,2], label='envMAW', color='blue', marker='o')
#     plt.errorbar(MAW_x_values, MAW_mean_std[0,:], yerr=MAW_mean_std[1,:], fmt='s', linestyle='--', ecolor='black', color='gray', capsize=5, label='MAW Mean ± Std')
#     plt.errorbar(comp_thrs_x_values, compensation_mean_array, yerr=compensation_std_array, fmt='o', linestyle='-', ecolor='black', color='green', capsize=5, label='Compensation\n'+set_name)
#     plt.errorbar(comp_thrs_x_values, threshold_pure_mean_array, yerr=threshold_pure_std_array, fmt='o', linestyle='-', ecolor='black', color='purple', capsize=5, label='Thrs and filt\n'+set_name)
    
    
#     # Adding labels and title
#     plt.xlabel('Exner time [-]')
#     plt.ylabel('Width %')
#     plt.title(f'envMAW, envBAW, MAW, thrs and comp widths - {set_name}')
#     plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
#     plt.ylim(0.0,1.0)
#     # Save the plot as a PDF file
#     plt.savefig(os.path.join(home_dir, 'output', f'REPORT_comp_thrs_filt_plot_{set_name}.pdf'), dpi=800, bbox_inches='tight')
#     plt.show()


# #%%
# # =============================================================================
# # 
# # =============================================================================
# for i, set_name in enumerate(set_names):
#     plt.figure(figsize=(10, 8))
    
#     if set_name == 'q15_3':
#         BAW_set_name = 'q15_2'
#     else:
#         BAW_set_name = set_name
    
#     # Load the mean and std arrays from text files
#     output_dir = os.path.join(home_dir, 'output', 'report_' + set_name, 'comp_thrs_analysis')
#     compensation_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_mean_array.txt'), delimiter=',')
#     compensation_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_std_array.txt'), delimiter=',')
#     compensation_pure_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_pure_mean_array.txt'), delimiter=',')
#     compensation_pure_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_pure_std_array.txt'), delimiter=',')
#     threshold_pure_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_threshold_pure_mean_array.txt'), delimiter=',')
#     threshold_pure_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_threshold_pure_std_array.txt'), delimiter=',')
#     MAW_mean_std = np.loadtxt(os.path.join(output_dir, set_name + '_MAW_mean_std.txt'), delimiter=',')

#     BAW_mean_std = np.loadtxt(os.path.join('/home/erri/Documents/PhD/Research/5_research_repos/PiQs_analysis/output_report/', BAW_set_name + '_envBAW_groups_report_merged_single_runs.txt'), delimiter=',', skiprows=0)
#     envMAW = np.loadtxt(os.path.join(home_dir,'output','report_'+ set_name, 'envelopes', set_name + '_timespan0_envMAW_report.txt'), delimiter=',')
    
#     BAW_x_values = BAW_mean_std[:,0]*np.max(MAW_x_values)/np.max(BAW_mean_std[:,0])
    
#     # Plot the datasets
#     plt.errorbar(envBAW_mean_std[0,:]/Txnr, envBAW_mean_std[1,:], yerr=envBAW_mean_std[1,:], fmt='s', linestyle='--', ecolor='black', color='red', capsize=5, label='envBAW Mean ± Std')
#     plt.errorbar(MAW_x_values, envMAW[:,1], yerr=envMAW[:,2], label='envMAW', color='blue', marker='o')
#     plt.errorbar(MAW_x_values, MAW_mean_std[0,:], yerr=MAW_mean_std[1,:], fmt='s', linestyle='--', ecolor='black', color='gray', capsize=5, label='MAW Mean ± Std')
#     plt.errorbar(comp_thrs_x_values, compensation_mean_array, yerr=compensation_std_array, fmt='o', linestyle='-', ecolor='black', color='green', capsize=5, label='Compensation\n'+set_name)
#     plt.errorbar(comp_thrs_x_values, threshold_pure_mean_array, yerr=threshold_pure_std_array, fmt='o', linestyle='-', ecolor='black', color='purple', capsize=5, label='Thrs and filt\n'+set_name)
    
    
#     # Plot purple and green areas
#     # plt.fill_between(MAW_x_values[1:], MAW_mean_std[0,1:], MAW_mean_std[0,1:] + threshold_pure_mean_array, color='purple', alpha=0.5, label='Thrs and filt')
#     plt.fill_between(MAW_x_values[1:], MAW_mean_std[0,1:], MAW_mean_std[0,1:] - threshold_pure_mean_array + compensation_mean_array, color='green', alpha=0.5, label='Compensation - Thrs')

    
#     # Adding labels and title
#     plt.xlabel('Exner time [-]')
#     plt.ylabel('Width %')
#     plt.title(f'envMAW, envBAW, MAW, thrs and comp widths - {set_name}')
#     plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
#     plt.ylim(0.0,1.0)
#     # Save the plot as a PDF file
#     # plt.savefig(os.path.join(home_dir, 'output', f'REPORT_comp_thrs_filt_plot_{set_name}.pdf'), dpi=800, bbox_inches='tight')
#     plt.show()



# #%%
# # =============================================================================
# # # OVERALL PLOT WITH ALL SET_NAMES
# # =============================================================================
# plt.figure(figsize=(10, 8))

# for i, set_name in enumerate(set_names):
    
#     if set_name == 'q15_3':
#         BAW_set_name = 'q15_2'
#     else:
#         BAW_set_name = set_name
        
        
#     # Load the mean and std arrays from text files
#     output_dir = os.path.join(home_dir, 'output', 'report_' + set_name, 'comp_thrs_analysis')
#     compensation_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_mean_array.txt'), delimiter=',')
#     compensation_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_std_array.txt'), delimiter=',')
#     threshold_pure_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_threshold_pure_mean_array.txt'), delimiter=',')
#     threshold_pure_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_threshold_pure_std_array.txt'), delimiter=',')
#     MAW_mean_std = np.loadtxt(os.path.join(output_dir, set_name + '_MAW_mean_std.txt'), delimiter=',')
#     envBAW_mean_std = np.loadtxt(os.path.join('/home/erri/Documents/PhD/Research/5_research_repos/PiQs_analysis/output_report/', BAW_set_name + '_envBAW_merged_single_runs.txt'), delimiter=',', skiprows=0)
#     envMAW = np.loadtxt(os.path.join(home_dir,'output','report_'+ set_name, 'envelopes', set_name + '_timespan0_envMAW_report.txt'), delimiter=',')
    
#     # BAW_x_values = envBAW_mean_std[:,0]*np.max(MAW_x_values)/np.max(envBAW_mean_std[:,0])
    
#     # Assign colors from the extracted colormaps
#     thrs_color = thrs_colors[i]
#     comp_color = comp_colors[i]
#     maw_color = maw_colors[i]
#     envmaw_color = envmaw_colors[i]
#     baw_color = baw_colors[i]
    
#     # Plot the datasets
#     plt.errorbar(BAW_x_values, envBAW_mean_std[:,1], yerr=envBAW_mean_std[:,2]/2, fmt='s', linestyle='--', ecolor='black', color='red', capsize=5, label='envBAW Mean ± Std')
#     plt.errorbar(MAW_x_values, MAW_mean_std[0,:], yerr=MAW_mean_std[1,:]/2, fmt='s', linestyle='--', ecolor='black', color=maw_color, capsize=5, label=f'MAW {set_name} Mean ± Std')
#     plt.errorbar(comp_thrs_x_values, compensation_mean_array, yerr=compensation_std_array, fmt='o', linestyle='-', ecolor='black', color=comp_color, capsize=5, label=f'Compensation {set_name}')
#     plt.errorbar(comp_thrs_x_values, threshold_pure_mean_array, yerr=threshold_pure_std_array, fmt='o', linestyle='-', ecolor='black', color=thrs_color, capsize=5, label=f'Thrs and filt {set_name}')
#     plt.errorbar(MAW_x_values, envMAW[:,1], yerr=envMAW[:,2], fmt='o', linestyle='-', ecolor='black', color=envmaw_color, capsize=5, label=f'Thrs and filt {set_name}')
# # Adding labels and title
# plt.xlabel('Exner time [-]')
# plt.ylabel('Width %')
# plt.title('envBAW, MAW, thrs and comp widths - All datasets')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
# plt.ylim(0.0,1.0)
# # Save the overall plot as a PDF file
# plt.savefig(os.path.join(home_dir, 'output', 'REPORT_comp_thrs_filt_plot_overall.pdf'), dpi=800, bbox_inches='tight')
# plt.show()