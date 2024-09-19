'''
This script processes image stacks from various experimental runs.
It performs the following operations:

1. Imports necessary libraries and custom functions.
2. Defines and processes a list of runs:
    - Loads a mask image and applies downsampling.
    - Loads data stacks corresponding to each run.
    - Converts the stacks to boolean arrays.
    - Applies a domain mask to the stacks.
    - Performs spatial and temporal filtering on the stacks.
    - Saves the processed stacks to specified directories.

3. Defines and processes a list of merged runs:
    - Determines a code based on the set name.
    - Loads a mask image and applies downsampling.
    - Loads merged data stacks.
    - Converts the stacks to boolean arrays.
    - Applies a domain mask to the stacks.
    - Performs spatial and temporal filtering on the stacks.
    - Saves the processed stacks to specified directories.

The script is divided into two main sections:
    "SINGLE RUN SECTION" for individual runs and "MERGED SINGLE RUN SECTION"
    for merged runs.

Parameters:
- `runs`: List of runs to process individually.
- `set_names`: List of sets to process as merged runs.
- `downsampling_dim`: Dimension for linear downsampling of the mask.
- `thrs`: Threshold value for spatial and temporal filtering.

The script assumes the presence of certain directories and files such as
'mask.tif' and '.npy' files corresponding to the activity stacks for each run.

Modify the lists of runs and set names, as well as other parameters to
customize the processing.
'''


# =============================================================================
# SINGLE RUN SECTION
# =============================================================================

# IMPORT LIBRARIES ------------------------------------------------------------
import os
import numpy as np
from PIL import Image
from PiQs_BAW_func_v1 import *

# RUNS ------------------------------------------------------------------------
# runs = ['q05r1','q05r2','q05r3','q05r4','q05r5','q05r6','q05r7','q05r8','q05r9'
#         ,'q07r1','q07r2','q07r3','q07r4','q07r5','q07r6','q07r7','q07r8','q07r9'
#         ,'q10r1','q10r2','q10r3','q10r4','q10r5','q10r6','q10r7','q10r8','q10r9'
#         ,'q15r1','q15r2','q15r3','q15r4','q15r5','q15r6','q15r7','q15r8','q15r9'
#         ,'q20r1','q20r2','q20r3','q20r4','q20r5','q20r6','q20r7','q20r8','q20r9']

runs = ['q05r1','q05r2','q05r3','q05r4','q05r5','q05r6','q05r7','q05r8','q05r9']
# runs = ['q07r1','q07r2','q07r3','q07r4','q07r5','q07r6','q07r7','q07r8','q07r9']
# runs = ['q10r1','q10r2','q10r3','q10r4','q10r5','q10r6','q10r7','q10r8','q10r9']
# runs = ['q15r1','q15r2','q15r3','q15r4','q15r5','q15r6','q15r7','q15r8','q15r9']
# runs = ['q20r1','q20r2','q20r3','q20r4','q20r5','q20r6','q20r7','q20r8','q20r9']

# runs = ['q07rgm', 'q10rgm2', 'q15rgm2', 'q20rgm2']
# runs = ['q10rgm2', 'q15rgm2', 'q20rgm2']
# runs = ['q07rgm', 'q10rgm2', 'q15rgm2', 'q20rgm2']
# runs = ['q07rgm', 'q10rgm2', 'q15rgm2']
# runs = ['q07rgm']
# runs = ['q10rgm2']
# runs = ['q20rgm2']

# runs = ['q20rgm2']
# runs = ['q20r9']

# SETUP FOLDERS----------------------------------------------------------------
folder_home = os.getcwd() # Setup home folder

# SCRIPT PARAMETERS------------------------------------------------------------
downsampling_dim = 5

for run in runs:

    print('****************')
    print(run)
    print('****************')
    
    # LOAD MASK ---------------------------------------------------------------
    mask_path = os.path.join(folder_home, 'mask.tif') # Define image path
    mask = Image.open(mask_path) # Open image as image
    mask_arr = np.array(mask)
    mask_arr_LR = non_overlapping_average(mask_arr, kernel_size=downsampling_dim) # perform linear downsampling

    # DEFINE PATHS ------------------------------------------------------------
    path_report = os.path.join(folder_home, 'output_report', run)

    if not(os.path.exists(path_report)):
        os.mkdir(path_report)
        
    # DEFINE STACK PATH--------------------------------------------------------
    stack_path = os.path.join(folder_home, 'activity_stack', run + '_BAA_stack_LR' + str(downsampling_dim) +  '.npy')
    
    # LOAD THE DATA STACK------------------------------------------------------
    stack = np.load(stack_path)
    
    # CONVERT THE STACK TO BOOLEAN---------------------------------------------
    stack_bool = np.where(stack>0, 1, stack)

    # APPLY DOMAIN MASK--------------------------------------------------------
    dim_t, dim_y, dim_x = stack_bool.shape # Define dimension

    # RESHAPE MASK-------------------------------------------------------------
    mask_arr_LR = mask_arr_LR[:dim_y,:dim_x] # Reshape
    mask_arr_LR = np.where(mask_arr_LR>0, 1, np.nan) # Set np.nan outside the domain

    # APPLY DOMAIN MASK--------------------------------------------------------
    stack_bool[:] = stack_bool[:]*mask_arr_LR
    stack_bool = np.where(stack_bool<0, np.nan, stack_bool)
        
    stack_bool_raw = np.copy(stack_bool)

    # =========================================================================
    # PERFORM SPATIAL AND TEMPORAL FILTERING
    # =========================================================================

    thrs = 0.4
    stack_bool_cld = spatial_temporal_activation_filtering(stack_bool_raw, (3,3,3), thrs)
    np.save(os.path.join(os.path.join(folder_home, 'activity_stack/activity_stack_cleaned'), run + '_BAA_stack_LR' + str(downsampling_dim) + '_cld.npy'),np.around(stack_bool_cld, decimals=0))

#%%
# =============================================================================
# MERGED SINGLE RUN SECTION
# =============================================================================

# IMPORT LIBRARIES ------------------------------------------------------------
import os
import numpy as np
from PIL import Image
from PiQs_BAW_func_v1 import *

# set_names = ['q07_1', 'q10_2', 'q15_2', 'q20_2']
set_names = ['q07_1']

for set_name in set_names:
    if '07' in set_name:
        code = '07'
    if '10' in set_name:
        code = '10'
    if '15' in set_name:
        code = '15'
    if '20' in set_name:
        code = '20'

    print('SINGLE RUN MERGED STACK FILTERING...')
    print(set_name)
    print()
    
    # SETUP FOLDERS------------------------------------------------------------
    folder_home = os.getcwd() # Setup home folder
    
    # SCRIPT PARAMETERS--------------------------------------------------------
    downsampling_dim = 5
    
    # LOAD MASK ---------------------------------------------------------------
    mask_path = os.path.join(folder_home, 'mask.tif') # Define image path
    mask = Image.open(mask_path) # Open image as image
    mask_arr = np.array(mask)
    mask_arr_LR = non_overlapping_average(mask_arr, kernel_size=downsampling_dim) # perform linear downsampling

    # DEFINE PATHS ------------------------------------------------------------
    path_report = os.path.join(folder_home, 'output_report', set_name)
    
    # IMPORT STACK ------------------------------------------------------------
    stack_path = os.path.join(os.getcwd(), 'activity_stack', set_name +
            '_single_run_merged_envelope_BAA_stack_LR' + str(downsampling_dim) + '.npy')
    stack = np.load(stack_path)
    
    # CONVERT THE STACK TO BOOLEAN --------------------------------------------
    stack_bool = np.where(stack>0, 1, stack)

    # APPLY DOMAIN MASK -------------------------------------------------------
    dim_t, dim_y, dim_x = stack_bool.shape # Define dimension

    # RESHAPE MASK ------------------------------------------------------------
    mask_arr_LR = mask_arr_LR[:dim_y,:dim_x] # Reshape
    mask_arr_LR = np.where(mask_arr_LR>0, 1, np.nan) # Set np.nan outside the domain

    # APPLY DOMAIN MASK -------------------------------------------------------
    stack_bool[:] = stack_bool[:]*mask_arr_LR
    stack_bool = np.where(stack_bool<0, np.nan, stack_bool)
        
    stack_bool_raw = np.copy(stack_bool)
        
    # =========================================================================
    # PERFORM SPATIAL AND TEMPORAL FILTERING
    # =========================================================================
    thrs = 0.4
    stack_bool_cld = spatial_temporal_activation_filtering(stack_bool_raw, (3,3,3), thrs)
    
    # SAVE STACK DATA ---------------------------------------------------------
    np.save(os.path.join(folder_home, 'activity_stack/activity_stack_cleaned', set_name + '_single_run_merged_BAA_stack_LR' + str(downsampling_dim) + '_cld.npy'), np.around(stack_bool_cld, decimals=0))
    # np.save(os.path.join(folder_home, 'activity_stack/activity_stack_cleaned', set_name + '_single_run_merged_envBAA_stack_LR' + str(downsampling_dim) + '_cld.npy'),np.around(env_stack_bool_cld, decimals=0))