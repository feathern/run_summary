def plot_energy_v_radius(infile, pdf, magnetic=False, min_decades=2, max_decades=5, plot_breakdown=True, lunits = None, eunits = None, 
                          tunits = None, me_denom=''):
    """
    
    plot_energy_v_radius:  Plots kinetic and magnetic energy densities vs. radius from a time-averaged Shell_Avgs file.
                         Plots are saved to a pdf.
                         
    Input Parameters:
            infile         -- the name of the Shell_Avgs file to process.
            pdf            -- the PdfPages object corresponding to the pdf file to plot to
            magnetic       -- (optional; default = False) Set to True to plot magnetic energy densities as well
            min_decades    -- (optional; default = 2) minimum number of decades along the y-axis.
            max_decades    -- (optional; default = 4) maximum number of decades along the y-axis
            plot_breakdown -- (optional; default = True) If true, fluctuating and mean energy densities, 
                              along with their r, theta, phi breakdown are plotted.  If False, only total
                              energy density, fluctuating energy density, and mean energy density are shown
                              (single plot, no r,theta,phi breakdown).
            lunits         -- (optional; default = None) length units for dimensional runs (e.g., cm)
            tunits         -- (optional; default = None) time units for dimensional runs (e.g., s)            
            eunits         -- (optional; default = None) energy density units for dimensional runs      
                   
    """
    from matplotlib import pyplot as plt
    import numpy as np
    import numpy
    from rayleigh_diagnostics import Shell_Avgs
    from matplotlib import rcParams 

    sa= Shell_Avgs(infile,path='',time_avg=True)
    radius = sa.radius

    # Page size and margins
    pxsize = 8.5
    pysize = 11.0
    ymar = 1.0/pysize
    xmar = 1.0/pxsize
    
    yspace = ymar*0.8  # vertical spacing between plots
    xsize = (1.0-2*xmar)*0.65 # length of plots
    # Note that the plot height is determined further down 
    # based on how many things are being plotted.

    ######################################################
    # Plot labelling
    labels = []
    qcodes = [] 
    ptitles = []
    

    # Kinetic Energy Codes and Labels
    ke_labels = []
    ke_codes = []
    ke_titles = []
    
    tlabels = [r'KE  = $\frac{1}{2V}\int_V\,\,\rho\left|\mathbf{u}\right|^2dV$ ', 
               r'MKE = $\frac{1}{2V}\int_V\,\,\rho\left|\overline{\mathbf{u}}\right|^2dV$',
               r'FKE = $\frac{1}{2V}\int_V\,\,\rho\left|\mathbf{u}-\overline{\mathbf{u}}\right|^2dV$']

    tqcodes = [401, 405, 409]
    

    
    ke_labels.append(tlabels)
    ke_codes.append(tqcodes)
    ke_titles.append('Kinetic Energy Density')
    

    if (plot_breakdown):
        mlabels = [r'MKE  = $\frac{1}{2V}\int_V\,\,\rho\left|\overline{\mathbf{u}}\right|^2dV$ ', 
        r'MRKE = $\frac{1}{2V}\int_V\,\,\rho\overline{u_r}^2dV$', 
        r'MTKE = $\frac{1}{2V}\int_V\,\,\rho\overline{u_\theta}^2dV$', 
        r'MPKE = $\frac{1}{2V}\int_V\,\,\rho\overline{u_\phi}^2dV$']
        mqcodes = [405, 406, 407, 408]
        ke_labels.append(mlabels)
        ke_codes.append(mqcodes)
        ke_titles.append('Mean Kinetic Energy Density')
   
        flabels = [r'FKE  = $\frac{1}{2V}\int_V\,\,\rho\left|\mathbf{u}-\overline{\mathbf{u}}\right|^2dV$ ', 
                   r'FRKE = $\frac{1}{2V}\int_V\,\,\rho\left(u_r-\overline{u_r}\right)^2dV$', 
                   r'FTKE = $\frac{1}{2V}\int_V\,\,\rho\left(u_\theta-\overline{u_\theta}\right)^2dV$', 
                   r'FPKE = $\frac{1}{2V}\int_V\,\,\rho\left(u_\phi-\overline{u_\phi}\right)^2dV$']        
        fqcodes = [409, 410, 411, 412]
        ke_labels.append(flabels)
        ke_codes.append(fqcodes)
        ke_titles.append('Fluctuating Kinetic Energy Density')

    labels.append(ke_labels)
    qcodes.append(ke_codes)
    ptitles.append(ke_titles)

    #Magnetic energy codes and labels.
    if (magnetic):
        me_labels = []
        me_codes  = []
        me_titles = []

        tlabels = [r'ME  = $\frac{1}{'+me_denom+r' V}\int_V\,\,\left|\mathbf{B}\right|^2dV$ ', 
           r'MME = $\frac{1}{'+me_denom+r' V}\int_V\,\,\left|\overline{\mathbf{B}}\right|^2dV$',
           r'FME = $\frac{1}{'+me_denom+r' V}\int_V\,\,\left|\mathbf{B}-\overline{\mathbf{B}}\right|^2dV$']

        me_labels.append(tlabels)
        me_codes.append([1101,1105, 1109])
        me_titles.append('Magnetic Energy Density')

        if (plot_breakdown):
            mlabels = [r'MME  = $\frac{1}{'+me_denom+r' V}\int_V\,\,\left|\overline{\mathbf{B}}\right|^2dV$ ', 
                       r'MRME = $\frac{1}{'+me_denom+r' V}\int_V\,\,\overline{B_r}^2dV$', 
                       r'MTME = $\frac{1}{'+me_denom+r' V}\int_V\,\,\overline{B_\theta}^2dV$', 
                       r'MPME = $\frac{1}{'+me_denom+r' V}\int_V\,\,\overline{B_\phi}^2dV$']
            mqcodes = [1105, 1106, 1107, 1108]
            me_labels.append(mlabels)
            me_codes.append(mqcodes)
            me_titles.append('Mean Magnetic Energy Density')
       

            flabels = [r'FME  = $\frac{1}{'+me_denom+r' V}\int_V\,\,\left|\mathbf{B}-\overline{\mathbf{B}}\right|^2dV$ ', 
                       r'FRME = $\frac{1}{'+me_denom+r' V}\int_V\,\,\left(B_r-\overline{B_r}\right)^2dV$', 
                       r'FTME = $\frac{1}{'+me_denom+r' V}\int_V\,\,\left(B_\theta-\overline{B_\theta}\right)^2dV$', 
                       r'FPME = $\frac{1}{'+me_denom+r' V}\int_V\,\,\left(B_\phi-\overline{B_\phi}\right)^2dV$']    
            fqcodes = [1109, 1110, 1111, 1112]
            me_labels.append(flabels)
            me_codes.append(fqcodes)
            me_titles.append('Fluctuating Magnetic Energy Density')    
            
        labels.append(me_labels)
        qcodes.append(me_codes)
        ptitles.append(me_titles)        


    # X- and Y- Axis Labels
    xlabel = r'Radius'
    if (lunits != None):
        xlabel = xlabel+' ('+lunits+')'
    
    ylabel = r'Energy Density'
    if (eunits != None):
        ylabel = ylabel+' ('+eunits+')'

    if (tunits != None):
        tstring = ' '+tunits
    else:
        tstring=''

    #########################################################
    # Create a message regarding temporal averaging
    txt="Time-averaging information:\n\n"
    txt=txt+"  first timestep: "+"{:,}".format(sa.iters[0])+'\n'
    txt=txt+"  last timestep : "+"{:,}".format(sa.iters[1])+'\n\n'
    txt=txt+"  start time    : "+"{:.4e}".format(sa.time[0])+tstring+'\n'
    txt=txt+"  end time      : "+"{:.4e}".format(sa.time[1])+tstring+'\n\n'    
    txt=txt+"  delta time    : "+"{:.4e}".format(sa.time[1]-sa.time[0])+tstring    
    
    ##################################################################################
    # Plotting
    for k,qtemp in enumerate(qcodes):   # Iterate over Kinetic and Magnetic energy densities
        newfig = True
        if (magnetic and not plot_breakdown):
            if (k > 0):
                newfig = False    
        if (newfig):
            fig = plt.figure(figsize=(pxsize,pysize))     
            time_info_shown = False  

        nplot = len(qtemp)
        if (nplot < 2):   # In the case of a single plot, prevent it from filling the full vertical page height.
            nplot = 2
        ysize = (1.0-2*ymar-(nplot-1)*yspace)/nplot 


        for j, qcode in enumerate(qtemp):

            yone = 1-ymar-(j+1)*ysize-j*yspace
            if (not newfig): # We will put ME and KE on the same page -- shift plot down
                yone = yone-yspace-ysize
            ax1 = fig.add_axes((xmar, yone, xsize,ysize))

            ax1.set_xlabel(xlabel)
            ax1.set_ylabel(ylabel)
            ax1.set_title(ptitles[k][j])
            
            maxvals = []
            minvals = []
            rcParams.update({'font.size': 10})     
            for i in range(len(qcode)):
                ax1.plot(radius, sa.vals[:,0,sa.lut[qcode[i]],0], label=labels[k][j][i])
                maxvals.append(numpy.max(sa.vals[:,0,sa.lut[qcode[i]],0]))   
                minvals.append(numpy.min(sa.vals[:,0,sa.lut[qcode[i]],0]))
                       
            ax1.set_yscale('log')


            # Attemp to do some sensible scaling on the y-axis
            ymax = numpy.max(maxvals)
            ymin = numpy.min(minvals)

            ylogtop = int(numpy.log10(ymax))+1
            ybottemp = int(numpy.log10(ymin))

            if (min_decades > 0):
                ylogbot = ylogtop-min_decades
                
                if (ylogbot > ybottemp):
                    ylogbot = ybottemp
            else:
                ylogbot = ybottemp

            ydiff = ylogtop-ylogbot
            if (ydiff > max_decades):
                ylogbot = ylogtop-max_decades

            ax1.set_ylim([10**ylogbot,10**ylogtop])

            ax1.legend(bbox_to_anchor=(1,1.05),loc='upper left', prop={'size' : 12})
            if (j == 0 and not time_info_shown):
                fig.text(xmar+xsize*1.05, 1-ymar-ysize*1.05, txt)
                time_info_shown = True
                
        # Finally, if we're not plotting the full energy breakdown, but we are plotting 
        # magnetism, we want to save this plot to the same page.
        savepage = True
        if (magnetic and not plot_breakdown):
            if (k == 0):
                savepage = False
        if (savepage):
            pdf.savefig()
            plt.close()

def plot_energy_flux(infile, pdf, magnetic=False, funits=None):
    """
    
    plot_energy_v_radius:  Plots kinetic and magnetic energy densities vs. radius from a time-averaged Shell_Avgs file.
                         Plots are saved to a pdf.
                         
    Input Parameters:
            infile         -- the name of the Shell_Avgs file to process.
            pdf            -- the PdfPages object corresponding to the pdf file to plot to
            magnetic       -- (optional; default = False) Set to True to plot magnetic energy densities as well
            min_decades    -- (optional; default = 2) minimum number of decades along the y-axis.
            max_decades    -- (optional; default = 4) maximum number of decades along the y-axis
            plot_breakdown -- (optional; default = True) If true, fluctuating and mean energy densities, 
                              along with their r, theta, phi breakdown are plotted.  If False, only total
                              energy density, fluctuating energy density, and mean energy density are shown
                              (single plot, no r,theta,phi breakdown).
            lunits         -- (optional; default = None) length units for dimensional runs (e.g., cm)
            tunits         -- (optional; default = None) time units for dimensional runs (e.g., s)            
            eunits         -- (optional; default = None) energy density units for dimensional runs      
                   
    """
    from matplotlib import pyplot as plt
    import numpy as np
    import numpy
    from rayleigh_diagnostics import Shell_Avgs
    from matplotlib import rcParams 

    sa= Shell_Avgs(infile,path='',time_avg=True)
    radius = sa.radius

    # Page size and margins
    pxsize = 8.5
    pysize = 11.0
    ymar = 1.0/pysize
    xmar = 1.0/pxsize
    
    yspace = ymar*0.8  # vertical spacing between plots
    xsize = (1.0-2*xmar)*0.65 # length of plots
    ysize = 0.4
    yone = 1-ymar-ysize
    
    
    xlabel = 'Radius'
    ylabel = r'Flux $\times\,4\pi r^2$'
    ptitle = 'Energy Flux Balance'

    #########################################################
    # Create a message regarding temporal averaging
    tstring=''
    txt="Time-averaging information:\n\n"
    txt=txt+"  first timestep: "+"{:,}".format(sa.iters[0])+'\n'
    txt=txt+"  last timestep : "+"{:,}".format(sa.iters[1])+'\n\n'
    txt=txt+"  start time    : "+"{:.4e}".format(sa.time[0])+tstring+'\n'
    txt=txt+"  end time      : "+"{:.4e}".format(sa.time[1])+tstring+'\n\n'    
    txt=txt+"  delta time    : "+"{:.4e}".format(sa.time[1]-sa.time[0])+tstring    

    
    qcodes = [1433, 1455, 1470, 1923, 1935]
    fig = plt.figure(figsize=(pxsize,pysize))     
    ax1 = fig.add_axes((xmar, yone, xsize,ysize))

    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title(ptitle)
            
    maxvals = []
    minvals = []
    rcParams.update({'font.size': 10})     
    fpr = 4*numpy.pi*radius**2
    tot = numpy.zeros(sa.nr,dtype='float64')
    for i in range(len(qcodes)):
        val = sa.vals[:,0,sa.lut[qcodes[i]],0]*fpr
        tot+= val
        ax1.plot(radius, val, label=str(i))
        maxvals.append(numpy.max(val))   
        minvals.append(numpy.min(sa.vals[:,0,sa.lut[qcodes[i]],0]))
    ax1.plot(radius,tot,label='Total')

    # Attemp to do some sensible scaling on the y-axis
    ymax = numpy.max(maxvals)
    ymin = numpy.min(minvals)

    ytop = 1.1*ymax
    ybottom = -.1*ymax
    if (ymin < ybottom):
        ybottom = ymin

    ax1.set_ylim([ybottom, ytop])

    ax1.legend(bbox_to_anchor=(1,1.05),loc='upper left', prop={'size' : 12})

    fig.text(xmar+xsize*1.05, 1-ymar-ysize*1.05, txt)

    pdf.savefig()
    plt.close()


