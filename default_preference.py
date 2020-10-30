#####################################################
# Copy this file as needed
# Each variable defined below must be defined

paper_name='Hindman, B.W., et al., 2020, ApJ, 898, 120\nMorphological Classification of the Convective Regimes in Rotating Stars'

magnetic = True

velocity = [1,2,3]
flux_bal = [1433, 1455, 1470, 1923, 1935]
thermal = [501]
ke = [401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412]


azavg_values = [3, 201, 202, 501]  # vphi, r/theta mass flux, entropy/temperature
shellavg_values = velocity+ke+thermal+flux_bal  # velocity, KE, entropy/temperature, Energy Fluxes -- add rheating flux, get sign right on visc
shell_spectra_values = [1,2,3,501] # velocity and entropy/temperature
globalavg_values = [401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412]  # all flavors of KE


max_iter = 900000000 # 900 million;  maximum timestep to search for.  Should be plenty..
num_shellavg = 15   # Number of Shell_Avgs files to time-average over.
num_azavg = 15      # Number of AZ_Avgs files to time-average over.
num_shellspectra = 15 # Number of Shell_Spectra files to time-average oer
num_gavg = -1      # Number of G_Avgs file to averages over (set to -1 to explicitly say 'all'; -1 only needed at the moment)

#In all cases, averaging is done over the last num_X files in the directory
#Each averaging routine also accepts a dt keyword, which confines averaging to the last dt's worth of time
#If both are set, which ever criteria is met first, dt or nfiles, controls the halt the averaging.

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
#me_denom = r'\mathrm{E}\,\mathrm{Pr}\,'
