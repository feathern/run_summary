# Summarization Script
from rayleigh_diagnostics import build_file_list, TimeAvg_AZAverages, AZ_Avgs, TimeAvg_ShellAverages, Shell_Avgs, Compile_GlobalAverages
from rayleigh_diagnostics import Shell_Spectra, TimeAvg_ShellSpectra
from snapshot_g_avgs import plot_energy_v_time
from snapshot_shell_avgs import plot_energy_v_radius, plot_energy_flux
from snapshot_shell_spectra import plot_spectra
from snapshot_shell_slices import plot_shell_slice
import sys
import importlib
import os

from title_page import title_page

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib

import subprocess as sp

pfile = 'template_preference'

prefs = importlib.import_module(pfile)

mnames = []
indirs = ['/home/feathern/runs/model_69_2']
outdirs = ['/home/feathern/runs/snapshot/model_69_2']

twait = 1000 # max time to wait on any command (seconds)




process_G_Avgs        = prefs.process_G_Avgs
process_Shell_Avgs    = prefs.process_Shell_Avgs 
process_AZ_Avgs       = prefs.process_AZ_Avgs 
process_Shell_Spectra = prefs.process_Shell_Spectra 
process_Shell_Slices  = prefs.process_Shell_Slices 
process_Checkpoints   = prefs.process_Checkpoints 

plot_G_Avgs        = prefs.plot_G_Avgs
plot_Shell_Avgs    = prefs.plot_Shell_Avgs
plot_AZ_Avgs       = prefs.plot_AZ_Avgs
plot_Shell_Spectra = prefs.plot_Shell_Spectra
plot_Shell_Slices  = prefs.plot_Shell_Slices




plot_data = plot_G_Avgs or plot_Shell_Avgs or plot_AZ_Avgs or plot_Shell_Spectra or plot_Shell_Slices
process_data = process_G_Avgs or process_Shell_Avgs or process_AZ_Avgs or process_Shell_Spectra or process_Shell_Slices


titles=[]
for i,indir in enumerate(indirs):
    input_dir= indir
    output_dir = outdirs[i]
    creation_cmd='mkdir -p '+output_dir
    s=sp.Popen(creation_cmd,shell=True)
    s.wait(timeout=twait)    
    
    snapdir = output_dir+'/snapshot_data'
    snapfile = output_dir+'/snapshot.pdf'
    mtemp = output_dir.split('/')
    model_name = mtemp[len(mtemp)-1]
    os.makedirs(snapdir, exist_ok=True)

    print('\n-------- Snapshotting Model: '+model_name+' ----------')
    print('---------- input dir : '+indir)
    print('---------- output dir: '+output_dir)


    if (process_data):
        # First, create the snapshot data

        if (process_G_Avgs):
            print('\n           Concatenating G_Avgs...')

            try:
                gafiles = build_file_list(0,prefs.max_iter,path=input_dir+'/G_Avgs')
            except:
                print('                 Error creating list of G_Avgs')
                gafiles = []
                
            if (len(gafiles) > 0):
                have_G_Avgs = True
                ga_ofile=snapdir+'/G_Avgs.dat'
                Compile_GlobalAverages(gafiles,ga_ofile, nfiles=prefs.num_G_Avgs, qcodes = prefs.globalavg_values) 
                if (plot_G_Avgs):
                    titles += ['Kinetic Energy vs. Time']+prefs.magnetic*['Magnetic Energy vs. Time']
            else:
                have_G_Avgs = False
                
            

        if (process_Shell_Avgs):
            print('           Time-Averaging Shell_Avgs...')    
            
            try:
                safiles = build_file_list(0,prefs.max_iter,path=input_dir+'/Shell_Avgs')
            except:
                print('                 Error creating list of Shell_Avgs')
                safiles = []
                
            if (len(safiles) > 0):     
                have_Shell_Avgs = True
                sa_ofile=snapdir+'/Shell_Avgs.dat'
                TimeAvg_ShellAverages(safiles,sa_ofile, qcodes=prefs.shellavg_values, dt = 1e8,  nfiles = prefs.num_Shell_Avgs)
                if (plot_Shell_Avgs):
                    titles += ['Kinetic Energy vs. Radius']+prefs.magnetic*['Magnetic Energy vs. Radius']    
            else:
                have_Shell_Avgs = False

        if (process_AZ_Avgs):
            print('           Time-Averaging AZ_Avgs...')    
            
            try:
                azfiles = build_file_list(0,prefs.max_iter,path=input_dir+'/AZ_Avgs')
            except:
                print('                 Error creating list of AZ_Avgs')
                azfiles = []
             
            if (len(azfiles) > 0):
                have_AZ_Avgs = True
                az_ofile=snapdir+'/AZ_Avgs.dat'
                TimeAvg_AZAverages(azfiles,az_ofile, qcodes=prefs.azavg_values, dt = 1e8,  nfiles = prefs.num_AZ_Avgs)                
                if (plot_AZ_Avgs):
                    titles += ['Kinetic Energy vs. Radius']+prefs.magnetic*['Magnetic Energy vs. Radius']        
            else:
                have_AZ_Avgs = False
                
        if (process_Shell_Spectra):
            print('           Time-Averaging Shell_Spectra...')           
            
            try:    
                spfiles = build_file_list(0,prefs.max_iter,path=input_dir+'/Shell_Spectra')
            except:
                print('                 Error creating list of Shell_Spectra')
                spfiles = []
            if (len(spfiles) > 0):
                have_Shell_Spectra = True
                sp_ofile=snapdir+'/Shell_Spectra.dat'
                TimeAvg_ShellSpectra(spfiles,sp_ofile, qcodes=prefs.shell_spectra_values, nfiles = prefs.num_Shell_Spectra)
                if (plot_Shell_Spectra):
                    titles += ['Kinetic Energy Spectrum']+prefs.magnetic*['Magnetic Energy Spectrum']               
            else:
                have_Shell_Spectra = False
                
        if (process_Shell_Slices):
            print('           Saving last Shell_Slices file...')
            try:
                shell_files = build_file_list(0,prefs.max_iter,path=input_dir+'/Shell_Slices')
            except:
                print('Error creating list of Shell_Slices')
                shell_files = []
            if (len(shell_files) > 0):
                have_Shell_Slices = True
                nshell = len(shell_files)
                last_shell = shell_files[nshell-1]
                shell_file = snapdir+'/Shell_Slices.dat'
                creation_cmd='cp '+last_shell+' '+shell_file
                s=sp.Popen(creation_cmd,shell=True)
                s.wait(timeout=twait)
                if (plot_Shell_Slices):
                    titles += ['Shell_Slice Snapshot']
            else:
                have_Shell_Slices = False
                
        if (process_Checkpoints):
            #Note:  need to add bette error checking into here.
            print('           Saving last checkpoint...')
            try:
                checkdir= output_dir+'/Checkpoints'
                os.makedirs(checkdir, exist_ok=True)            
                check_file = input_dir+'/Checkpoints/last_checkpoint'
                cfile = open(check_file,'r')
                eof = False
                clines=[]
                while (not eof):                            # Keep reading forever
                    line = cfile.readline()
                    ilen = len(line)
                    if (ilen == 0):
                        eof = True
                    else:
                        clines.append(line)
                if (len(clines) == 1):
                    check_pref = input_dir+'/Checkpoints/'+clines[0].split()[0]
                    check_pref2 = clines[0].split()[0]
                else:
                    check_pref = input_dir+'/Checkpoints/quicksave_'+clines[1].split()[0]
                    check_pref2 = 'quicksave_'+clines[1].split()[0]
                print('              ...Saving '+check_pref2)
                copycmd = 'cp '+check_file+' '+checkdir+'/.'
                s=sp.Popen(copycmd,shell=True)
                s.wait(timeout=twait)    
                
                copycmd = 'cp -r '+check_pref+' '+checkdir+'/.'
                s=sp.Popen(copycmd,shell=True)
                s.wait(timeout=twait)               

                copycmd = 'cp '+input_dir+'/main_input '+output_dir+'/.'
                s=sp.Popen(copycmd,shell=True)
                s.wait(timeout=twait)         
            except:
                print('                 ... Error... no checkpoint backed up.')

            
    #Now, make the plot if desired
    if (plot_data):
        with PdfPages(snapfile) as pdf:

            title_page(plt,titles, name = model_name, paper = prefs.paper_name)
            pdf.savefig()
            plt.close()
            
            if (have_G_Avgs and plot_G_Avgs):
                plot_energy_v_time(ga_ofile, pdf,magnetic=prefs.magnetic, tunits=prefs.time_units, eunits=prefs.energy_density_units, 
                                   me_denom = prefs.me_denom, plot_breakdown=prefs.plot_breakdown)
            if (have_Shell_Avgs and plot_Shell_Avgs):
                plot_energy_v_radius(sa_ofile, pdf, magnetic=prefs.magnetic, lunits=prefs.length_units,eunits=prefs.energy_density_units,
                                      me_denom=prefs.me_denom, tunits=prefs.time_units,plot_breakdown=prefs.plot_breakdown)       
                plot_energy_flux(sa_ofile, pdf, magnetic=False, funits=None)  
                           
            if (have_Shell_Spectra and plot_Shell_Spectra):       
                plot_spectra(sp_ofile,velocity=True,magnetic=prefs.magnetic,thermal=False, decades=4, pdf = pdf)
                
            if (have_Shell_Slices and plot_Shell_Slices):
                plot_shell_slice(shell_file,pdf, magnetic=prefs.magnetic)


