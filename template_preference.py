####################################################################
# Copy this file to snapshot_preference.py and modify as needed

# A string containing the name of the paper.  Multiple lines are fine.
paper_name='Hindman, B.W., et al., 2020, ApJ, 898, 120\nMorphological Classification of the Convective Regimes in Rotating Stars'

magnetic = False   # Set to true for magnetic runs if you wish to plot magnetic quantities (if output).

##############################################################################################
# Decide which data types we would like to post-processes and save.
process_G_Avgs        = True # Concatenate G_Avgs (or not)
process_Shell_Avgs    = True # Time-Average Shell_Avgs (or not)
process_AZ_Avgs       = True # Time-Average AZ_Avgs (or not)
process_Shell_Spectra = True # Time-Average Shell_Spectra (or not)
process_Shell_Slices  = True # Backup last Shell_Slices file (or not)
process_Checkpoints   = True # Backup last Checkpoint and main_input (or not)

##############################################################################################
# We can optionally plot each of the processed files,
# creating a multi-page pdf in the process.
plot_G_Avgs        = True 
plot_Shell_Avgs    = True 
plot_AZ_Avgs       = True # Not yet implemented
plot_Shell_Spectra = True 
plot_Shell_Slices  = False

##############################################################################################
#  Decide how much data to process.
#  For time-averaged outputs, averaging is carried out over the last num_X files found in each subdirectory.

max_iter = 900000000    # 900 million;  maximum iteration number to search for. 
                        # This is unused for checkpoints (last checkpoint always taken).
num_Shell_Avgs    = 15  # Number of Shell_Avgs files to time-average over.
num_AZ_Avgs       = 15  # Number of AZ_Avgs files to time-average over.
num_Shell_Spectra = 15  # Number of Shell_Spectra files to time-average oer
num_G_Avgs        = -1  # Number of G_Avgs files to concatenate (-1 indicates 'use all available files')



velocity = [1,2,3]
flux_bal = [1433, 1455, 1470, 1923, 1935]
thermal = [501]
ke = [401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412]

plot_breakdown = True  # Set to true to plot the different contributions to kinetic and magnetic energy

# Decide which values will be extracted from the different outputs and stored
# in the post-processed data.  Non-existent values will be stored as zero.
# All values here are required to make the full demonstration pdf.

azavg_values = [3, 201, 202, 501]  # vphi, r/theta mass flux, entropy/temperature
shellavg_values = velocity+ke+thermal+flux_bal  # velocity, KE, entropy/temperature, Energy Fluxes -- add rheating flux, get sign right on visc
shell_spectra_values = [1,2,3,501] # velocity and entropy/temperature
globalavg_values = [401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412]  # all flavors of KE

# Add magnetic quantity codes
if (magnetic):
    m_values = [1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110, 1111, 1112]  # all flavors of ME
    globalavg_values += m_values
    shellavg_values += m_values
    shell_spectra_values += [801, 802, 803]
    
    
# Units
# For nondimensional runs, set these to None (the NoneType, not the string 'None').
time_units = 's'  # time
energy_density_units = 'erg cm$^{-3}$'
length_units = 'cm'
me_denom = r'8\,\pi\,'  #latex codes OK here -- omit $-signs

