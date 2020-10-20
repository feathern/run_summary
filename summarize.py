# Summarization Script
from rayleigh_diagnostics import build_file_list, TimeAvg_AZAverages, AZ_Avgs, TimeAvg_ShellAverages, Shell_Avgs, Compile_GlobalAverages
from rayleigh_diagnostics import Shell_Spectra, TimeAvg_ShellSpectra
from snapshot_g_avgs import plot_energy_v_time
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
        TimeAvg_ShellAverages(safiles,sa_ofile, qcodes=prefs.shellavg_values, dt = 1e8,  nfiles = prefs.num_shellavg)
        #atest = Shell_Avgs(ofile,path='./', time_avg = True)

        
    #ofile='sptest.dat'
    #spfiles = build_file_list(0,9000000000,path=input_dir+'/Shell_Spectra')
    #TimeAvg_ShellSpectra(spfiles,ofile, qcodes=shell_spectra_values, nfiles = num_shellspectra)
    #atest = Shell_Spectra(ofile,path='./', time_avg = True)
    #print(atest.time)
    #print(atest.iters)
    
    #Now, begin to plot it
    titles = ['Kinetic Energy vs. Time']+prefs.magnetic*['Magnetic Energy vs. Time']
    with PdfPages(snapfile) as pdf:

        title_page(plt,titles, name = model_name, paper = prefs.paper_name)
        pdf.savefig()
        plt.close()
        plot_energy_v_time(ga_ofile, velocity=True, decades=0,magnetic=prefs.magnetic, pdf=pdf)
        #pdf.savefig()
        #plt.close()

