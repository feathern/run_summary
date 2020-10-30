# Summarization Script
from rayleigh_diagnostics import build_file_list, TimeAvg_AZAverages, AZ_Avgs, TimeAvg_ShellAverages, Shell_Avgs, Compile_GlobalAverages
from rayleigh_diagnostics import Shell_Spectra, TimeAvg_ShellSpectra
from snapshot_g_avgs import plot_energy_v_time
from snapshot_shell_avgs import plot_energy_v_radius, plot_energy_flux
from snapshot_shell_spectra import plot_spectra
import sys
import importlib
import os

from title_page import title_page

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib


pfile = 'default_preference'

prefs = importlib.import_module(pfile)

mnames = []
indirs = ['/home/feathern/runs/model_69_2']
outdirs = ['/home/feathern/runs/snapshot/model_69_2']


print(prefs.azavg_values)

have_azavg = len(prefs.azavg_values) > 0
have_gavg = len(prefs.globalavg_values) > 0
have_shellavg = len(prefs.shellavg_values) > 0
have_spectra = len(prefs.shell_spectra_values) > 0
debug = True

titles=[]
for i,indir in enumerate(indirs):
    input_dir= indir
    output_dir = outdirs[i]
    snapdir = output_dir+'/snapshot_data'
    snapfile = output_dir+'/snapshot.pdf'
    mtemp = output_dir.split('/')
    model_name = mtemp[len(mtemp)-1]
    os.makedirs(snapdir, exist_ok=True)

    # First, create the snapshot data
    if (have_gavg):
        gafiles = build_file_list(0,prefs.max_iter,path=input_dir+'/G_Avgs')
        ga_ofile=snapdir+'/G_Avgs.dat'
        # <------------ put some plot_breakdown logic here for the energy densities.
        titles += ['Kinetic Energy vs. Time']+prefs.magnetic*['Magnetic Energy vs. Time']
        if (not debug):
            Compile_GlobalAverages(gafiles,ga_ofile, nfiles=prefs.num_gavg, qcodes = prefs.globalavg_values)    
   
    #azfiles = build_file_list(0,9000000000,path=input_dir+'/AZ_Avgs')
    #ofile = 'aztest.dat'
    #TimeAvg_AZAverages(azfiles,ofile, qcodes=azavg_values, dt = 1e8, nfiles = num_azavg)
    #atest = AZ_Avgs(ofile,path='./', time_avg = True)
    #print(atest.time)
    #print(atest.iters)
    #ofile ='satest.dat'

    if (have_shellavg):
        safiles = build_file_list(0,prefs.max_iter,path=input_dir+'/Shell_Avgs')
        sa_ofile=snapdir+'/Shell_Avgs.dat'
        titles += ['Kinetic Energy vs. Radius']+prefs.magnetic*['Magnetic Energy vs. Radius']        
        if (not debug):
            TimeAvg_ShellAverages(safiles,sa_ofile, qcodes=prefs.shellavg_values, dt = 1e8,  nfiles = prefs.num_shellavg)


    if (have_spectra):
        spfiles = build_file_list(0,prefs.max_iter,path=input_dir+'/Shell_Spectra')
        sp_ofile=snapdir+'/Shell_Spectra.dat'
        titles += ['Kinetic Energy Spectrum']+prefs.magnetic*['Magnetic Energy Spectrum']
        if (not debug):
            TimeAvg_ShellSpectra(spfiles,sp_ofile, qcodes=prefs.shell_spectra_values, nfiles = prefs.num_shellspectra)

    
    #Now, begin to plot it
    
    with PdfPages(snapfile) as pdf:

        title_page(plt,titles, name = model_name, paper = prefs.paper_name)
        pdf.savefig()
        plt.close()
        plot_energy_v_time(ga_ofile, pdf,magnetic=prefs.magnetic, tunits=prefs.time_units, eunits=prefs.energy_density_units, 
                           me_denom = prefs.me_denom, plot_breakdown=True)
        plot_energy_v_radius(sa_ofile, pdf, magnetic=prefs.magnetic, lunits=prefs.length_units,eunits=prefs.energy_density_units,
                              me_denom=prefs.me_denom, tunits=prefs.time_units,plot_breakdown=True)       
        plot_energy_flux(sa_ofile, pdf, magnetic=False, funits=None)                    
        plot_spectra(sp_ofile,velocity=True,magnetic=prefs.magnetic,thermal=False, decades=10, pdf = pdf)
        #pdf.savefig()
        #plt.close()

