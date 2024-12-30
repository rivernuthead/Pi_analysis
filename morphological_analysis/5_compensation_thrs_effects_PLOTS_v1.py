# IMPORT LIBRARIES
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# =============================================================================
# SCRIPT PARAMETERS
# =============================================================================

# EXNER TIME ------------------------------------------------------------------
Txnr_array = [20, 30, 48, 74, 100]

# RUN DURATION [MINS] ---------------------------------------------------------
run_time_array = [15, 19, 25, 47, 81]  # Run duration in minutes

# RUN DURATION [FRAMES] -------------------------------------------------------
run_frames = [10, 15, 20, 40, 88]  # Run duration in number of frames

# EXNER SCALE -----------------------------------------------------------------
Txnr_scale_array = np.array(run_time_array)/np.array(Txnr_array)

# TIMESCALE -------------------------------------------------------------------
env_tscale_array = [2, 3, 4, 9, 22]  # 1/4 of the frame number (~1/8*Txnr)

#==============================================================================
# PLOTS
# =============================================================================

set_names = ['q05_1', 'q07_1', 'q10_2', 'q15_3', 'q20_2']
# set_names = ['q05_1']

# set_names = ['q05_1']
# set_names = ['q07_1']
# set_names = ['q10_2']
# set_names = ['q20_2']

# FOLDER SETUP ----------------------------------------------------------------
home_dir = os.path.join(os.getcwd(), 'morphological_analysis') # Home directory

# FUNCTIONS -------------------------------------------------------------------

# Function to extract the last 80% of a colormap
def extract_colormap(cmap_name, num_colors):
    cmap = cm.get_cmap(cmap_name)
    colors = cmap(np.linspace(0.4, 1, num_colors))
    return colors

# Generate colormaps for each type of plot
thrs_colors = extract_colormap('Purples', len(set_names))
comp_colors = extract_colormap('Greens', len(set_names))
maw_colors = extract_colormap('Greys', len(set_names))
envmaw_colors = extract_colormap('Blues', len(set_names))
baw_colors = extract_colormap('Reds', len(set_names))


# =============================================================================
# # PLOT FOR EACH SET_NAME
# =============================================================================
for i, set_name in enumerate(set_names):
    
    if set_name == 'q15_3':
        BAW_set_name = 'q15_2'
    else:
        BAW_set_name = set_name
        

    if '05' in set_name:
        skip = 1
        env_tscale = env_tscale_array[4]
        Txnr = Txnr_array[4]
        run_time = run_time_array[4]
        run_frm = run_frames[4]
        Txnr_scale = Txnr_scale_array[4]
        rgm_code = 'q05rgm4'
        rgm_skip_frame = 0
        MAW_set_name = 'q05_1'

    if '07' in set_name:
        skip = 1
        env_tscale = env_tscale_array[3]
        Txnr = Txnr_array[3]
        run_time = run_time_array[3]
        run_frm = run_frames[3]
        Txnr_scale = Txnr_scale_array[3]
        rgm_code = 'q07rgm'
        rgm_skip_frame = 0
        MAW_set_name = 'q07_1'

    if '10' in set_name:
        skip = 1
        env_tscale = env_tscale_array[2]
        Txnr = Txnr_array[2]
        run_time = run_time_array[2]
        run_frm = run_frames[2]
        Txnr_scale = Txnr_scale_array[2]
        rgm_code = 'q10rgm2'
        rgm_skip_frame = 0
        MAW_set_name = 'q10_2'

    if '15' in set_name:
        skip = 1
        env_tscale = env_tscale_array[1]
        Txnr = Txnr_array[1]
        run_time = run_time_array[1]
        run_frm = run_frames[1]
        Txnr_scale = Txnr_scale_array[1]
        rgm_code = 'q15rgm2'
        rgm_skip_frame = 5
        MAW_set_name = 'q15_3'

    if '20' in set_name:
        skip = 1
        env_tscale = env_tscale_array[0]
        Txnr = Txnr_array[0]
        run_time = run_time_array[0]
        run_frm = run_frames[0]
        Txnr_scale = Txnr_scale_array[0]
        rgm_code = 'q20rgm2'
        rgm_skip_frame = 0
        MAW_set_name = 'q20_2'
    
    # SET FOLDERS
    env_dir = os.path.join(home_dir,'output_data', 'envelopes')
    if not(os.path.exists(env_dir)):
        os.mkdir(env_dir)
    # IMPORT STACK ------------------------------------------------------------
    stack_dir = os.path.join(home_dir, 'output_data', 'DoDs', 'DoDs_stack') # Define the stack directory
    stack=np.load(os.path.join(stack_dir, 'DoD_stack_'+set_name+'.npy')) # Load the stack
    stack_bool=np.load(os.path.join(stack_dir, 'DoD_stack_bool_'+set_name+'.npy'))

    # COMPUTE X-VALUES ARRAYS------------------------------------------------------
    # Compensation and threshold
    comp_thrs_x_values = np.linspace(2, stack.shape[0], stack.shape[0]-1) * Txnr_scale

    # MAW
    MAW_x_values = np.linspace(1, stack.shape[0], stack.shape[0])*Txnr_scale
    
    # Load the mean and std arrays from text files
    output_dir = os.path.join(home_dir, 'output_data', 'comp_thrs_analysis')
    compensation_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_mean_array.txt'), delimiter=',')
    compensation_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_std_array.txt'), delimiter=',')
    compensation_pure_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_pure_mean_array.txt'), delimiter=',')
    compensation_pure_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_pure_std_array.txt'), delimiter=',')
    threshold_pure_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_threshold_pure_mean_array.txt'), delimiter=',')
    threshold_pure_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_threshold_pure_std_array.txt'), delimiter=',')
    MAW_mean_std = np.loadtxt(os.path.join(output_dir, set_name + '_MAW_mean_std.txt'), delimiter=',')

    envBAW_mean_std = np.loadtxt(os.path.join(os.getcwd(), 'bedload_analysis','output_data', 'output_report', BAW_set_name + '_envBAW_merged_single_runs.txt'), delimiter=',', skiprows=0)
    envMAW = np.loadtxt(os.path.join(env_dir, set_name + '_timespan0_envMAW_report.txt'), delimiter=',')
    
    BAW_x_values = np.linspace(0, np.max(MAW_x_values), len(envBAW_mean_std[0,:]))
    
    # Resample envBAW to match MAW_x_values
    envBAW_y_resampled = np.interp(MAW_x_values, envBAW_mean_std[0,:], envBAW_mean_std[1,:])
    envBAW_std_resampled = np.interp(MAW_x_values, envBAW_mean_std[0,:], envBAW_mean_std[2,:])
    plt.figure(figsize=(20, 20))
    # Plot the datasets
    plt.errorbar(envBAW_mean_std[0,:], envBAW_mean_std[1,:], yerr=envBAW_mean_std[2,:], fmt='s', linestyle='--', ecolor='black', color='red', capsize=5, label='envBAW Mean ± Std')
    plt.errorbar(MAW_x_values, envMAW[:,1], yerr=envMAW[:,2], label='envMAW', color='blue', marker='o')
    plt.errorbar(MAW_x_values, MAW_mean_std[0,:], yerr=MAW_mean_std[1,:], fmt='s', linestyle='--', ecolor='black', color='gray', capsize=5, label='MAW Mean ± Std')
    plt.errorbar(comp_thrs_x_values, compensation_mean_array, yerr=compensation_std_array, fmt='o', linestyle='-', ecolor='black', color='green', capsize=5, label='Compensation\n'+set_name)
    plt.errorbar(comp_thrs_x_values, threshold_pure_mean_array, yerr=threshold_pure_std_array, fmt='o', linestyle='-', ecolor='black', color='purple', capsize=5, label='Thrs and filt\n'+set_name)
    
    plt.xlabel('Exner time [-]')
    plt.ylabel('Width %')
    plt.title(f'envMAW, envBAW, MAW, thrs and comp widths - {set_name}')
    plt.legend(bbox_to_anchor=(0.85, 0.4), loc='upper left')
    plt.ylim(0.0,1.0)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.text(0.5, -0.12, f"Generated by: {os.path.basename(__file__)}", transform=plt.gca().transAxes, ha='center', va='center', fontsize=12)
    plt.savefig(os.path.join(home_dir, 'output_data', f'comp_thrs_filt_plot_{set_name}.pdf'), dpi=800, bbox_inches='tight')
    plt.show()
    
    # Prepare data for report
    report_lines = []
    
    # Adding each curve data to the report
    curves = [
        ("envBAW Mean ± Std", envBAW_mean_std[0,:], envBAW_mean_std[1,:], envBAW_mean_std[2,:]),
        ("envBAW (Resampled) Mean ± Std", MAW_x_values, envBAW_y_resampled, envBAW_std_resampled),
        ("envMAW", MAW_x_values, envMAW[:,1], envMAW[:,2]),
        ("MAW Mean ± Std", MAW_x_values, MAW_mean_std[0,:], MAW_mean_std[1,:]),
        (f"Compensation {set_name}", comp_thrs_x_values, compensation_mean_array, compensation_std_array),
        (f"Thrs and filt {set_name}", comp_thrs_x_values, threshold_pure_mean_array, threshold_pure_std_array)
    ]
    
    for name, x_values, y_values, std_values in curves:
        report_lines.append(f"Curve: {name}")
        report_lines.append("X Values, Y Values, Std Dev")
        for x, y, std in zip(x_values, y_values, std_values):
            report_lines.append(f"{x:.4f}, {y:.4f}, {std:.4f}")
        report_lines.append("\n")  # Add a blank line between different curves
    
    # Write the report to a txt file
    report_file_path = os.path.join(home_dir, 'output_data', f'curve_report_{set_name}.txt')
    with open(report_file_path, 'w') as report_file:
        report_file.write("\n".join(report_lines))
    
    print(f"Report saved to {report_file_path}")

#%%
# =============================================================================
# 
# =============================================================================
# for i, set_name in enumerate(set_names):    
#     if set_name == 'q15_3':
#         BAW_set_name = 'q15_2'
#     else:
#         BAW_set_name = set_name
    
#     # Load the mean and std arrays from text files
#     output_dir = os.path.join(home_dir, 'output_data', 'comp_thrs_analysis')
#     compensation_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_mean_array.txt'), delimiter=',')
#     compensation_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_std_array.txt'), delimiter=',')
#     compensation_pure_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_pure_mean_array.txt'), delimiter=',')
#     compensation_pure_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_pure_std_array.txt'), delimiter=',')
#     threshold_pure_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_threshold_pure_mean_array.txt'), delimiter=',')
#     threshold_pure_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_threshold_pure_std_array.txt'), delimiter=',')
#     MAW_mean_std = np.loadtxt(os.path.join(output_dir, set_name + '_MAW_mean_std.txt'), delimiter=',')

#     BAW_mean_std = np.loadtxt(os.path.join('/home/erri/Documents/PhD/Research/5_research_repos/PiQs_analysis/output_report/', BAW_set_name + '_envBAW_groups_report_merged_single_runs.txt'), delimiter=',', skiprows=0)
#     envMAW = np.loadtxt(os.path.join(home_dir,'output_data', 'envelopes', set_name + '_timespan0_envMAW_report.txt'), delimiter=',')
    
#     BAW_x_values = np.linspace(0, np.max(MAW_x_values), len(envBAW_mean_std[0,:]))
#     plt.figure(figsize=(20, 20))
#     # Plot the datasets
#     plt.errorbar(envBAW_mean_std[0,:]/Txnr, envBAW_mean_std[1,:], yerr=envBAW_mean_std[2,:], fmt='s', linestyle='--', ecolor='black', color='red', capsize=5, label='envBAW Mean ± Std')
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
#     plt.legend(bbox_to_anchor=(0.85, 0.4), loc='upper left')
#     plt.ylim(0.0,1.0)
#     plt.text(0.5, -0.8, f"Generated by: {os.path.basename(__file__)}", transform=plt.gca().transAxes, ha='center', va='center', fontsize=6)
#     # plt.savefig(os.path.join(home_dir, 'output_data', f'REPORT_comp_thrs_filt_plot_{set_name}.pdf'), dpi=800, bbox_inches='tight')
#     plt.show()



#%%
# =============================================================================
# # OVERALL PLOT WITH ALL SET_NAMES
# =============================================================================
plt.figure(figsize=(20, 20))

for i, set_name in enumerate(set_names):
    
    if set_name == 'q15_3':
        BAW_set_name = 'q15_2'
    else:
        BAW_set_name = set_name
        
        
    # Load the mean and std arrays from text files
    output_dir = os.path.join(home_dir, 'output_data', 'comp_thrs_analysis')
    compensation_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_mean_array.txt'), delimiter=',')
    compensation_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_compensation_std_array.txt'), delimiter=',')
    threshold_pure_mean_array = np.loadtxt(os.path.join(output_dir, set_name + '_threshold_pure_mean_array.txt'), delimiter=',')
    threshold_pure_std_array = np.loadtxt(os.path.join(output_dir, set_name + '_threshold_pure_std_array.txt'), delimiter=',')
    MAW_mean_std = np.loadtxt(os.path.join(output_dir, set_name + '_MAW_mean_std.txt'), delimiter=',')
    envBAW_mean_std = np.loadtxt(os.path.join(os.getcwd(), 'bedload_analysis', 'output_data', 'output_report', BAW_set_name + '_envBAW_merged_single_runs.txt'), delimiter=',', skiprows=0)
    envMAW = np.loadtxt(os.path.join(home_dir,'output_data', 'envelopes', set_name + '_timespan0_envMAW_report.txt'), delimiter=',')
    
    BAW_x_values = np.linspace(0, np.max(MAW_x_values), len(envBAW_mean_std[0,:]))
    
    # Assign colors from the extracted colormaps
    thrs_color = thrs_colors[i]
    comp_color = comp_colors[i]
    maw_color = maw_colors[i]
    envmaw_color = envmaw_colors[i]
    baw_color = baw_colors[i]
    
    # Plot the datasets
    plt.errorbar(BAW_x_values, envBAW_mean_std[1,:], yerr=envBAW_mean_std[2,:], fmt='s', linestyle='--', ecolor='black', color='red', capsize=5, label='envBAW Mean ± Std')
    plt.errorbar(MAW_x_values, MAW_mean_std[0,:], yerr=MAW_mean_std[1,:]/2, fmt='s', linestyle='--', ecolor='black', color=maw_color, capsize=5, label=f'MAW {set_name} Mean ± Std')
    plt.errorbar(comp_thrs_x_values, compensation_mean_array, yerr=compensation_std_array, fmt='o', linestyle='-', ecolor='black', color=comp_color, capsize=5, label=f'Compensation {set_name}')
    plt.errorbar(comp_thrs_x_values, threshold_pure_mean_array, yerr=threshold_pure_std_array, fmt='o', linestyle='-', ecolor='black', color=thrs_color, capsize=5, label=f'Thrs and filt {set_name}')
    plt.errorbar(MAW_x_values, envMAW[:,1], yerr=envMAW[:,2], fmt='o', linestyle='-', ecolor='black', color=envmaw_color, capsize=5, label=f'Thrs and filt {set_name}')
# Adding labels and title
plt.xlabel('Exner time [-]')
plt.ylabel('Width %')
plt.title('envBAW, MAW, thrs and comp widths - All datasets')
plt.legend(bbox_to_anchor=(0.85, 0.4), loc='upper left')
plt.ylim(0.0,1.0)
# Save the overall plot as a PDF file
plt.savefig(os.path.join(home_dir, 'output_data', 'REPORT_comp_thrs_filt_plot_overall.pdf'), dpi=800, bbox_inches='tight')
plt.show()