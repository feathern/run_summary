def plot_energy_v_time(infile, velocity=True,magnetic=False, decades=5, pdf=None,plot_breakdown=True):
    # Provide full path to infile
    from matplotlib import pyplot as plt
    import numpy as np
    import numpy
    from rayleigh_diagnostics import G_Avgs

    if (magnetic):
        velocity = False

    ga= G_Avgs(infile,path='')
    time = ga.time
    nt = len(ga.time)
    energy = numpy.zeros((nt,3),dtype='float64')

    pxsize = 8.5
    pysize = 11.0
    ymar = 1.0/pysize
    xmar = 1.0/pxsize
    yspace = ymar*0.75


    # First, set up the labeling for the different plots
    labels = []
    qcodes = [] 
    ptitles = []
    ###########################
    # KE Codes

    ke_labels = []
    ke_codes = []
    ke_titles = []
    
    tlabels = ['Total KE', 'Axisymmetric KE', 'Non-Axisymmetric KE']
    tqcodes = [401, 405, 409]
    

    
    ke_labels.append(tlabels)
    ke_codes.append(tqcodes)
    ke_titles.append('Kinetic Energy Density')
    
    plot_breakdown = True

    if (plot_breakdown):
        mlabels = ['Mean KE', 'Mean Radial KE', r'Mean $\theta$ KE', r'Mean $\phi$ KE']
        mqcodes = [405, 406, 407, 408]
        ke_labels.append(mlabels)
        ke_codes.append(mqcodes)
        ke_titles.append('Mean Kinetic Energy Density')
   
        flabels = ['Fluctuating KE', 'Fluctuating Radial KE', r'Fluctuating $\theta$ KE', r'Fluctuating $\phi$ KE']
        fqcodes = [409, 410, 411, 412]
        ke_labels.append(flabels)
        ke_codes.append(fqcodes)
        ke_titles.append('Fluctuating Kinetic Energy Density')

    labels.append(ke_labels)
    qcodes.append(ke_codes)
    ptitles.append(ke_titles)

   
    if (magnetic):
        me_labels = []
        me_codes  = []
        me_titles = []
        
        me_labels.append(['Total ME', 'Axisymmetric ME', 'Non-Axisymmetric ME'])
        me_codes.append([1101,1105, 1109])
        me_titles.append('Magnetic Energy Density')

        if (plot_breakdown):
            mlabels = ['Mean ME', 'Mean Radial ME', r'Mean $\theta$ ME', r'Mean $\phi$ ME']
            mqcodes = [1105, 1106, 1107, 1108]
            me_labels.append(mlabels)
            me_codes.append(mqcodes)
            me_titles.append('Mean Magnetic Energy Density')
       
            flabels = ['Fluctuating ME', 'Fluctuating Radial ME', r'Fluctuating $\theta$ ME', r'Fluctuating $\phi$ ME']
            fqcodes = [1109, 1110, 1111, 1112]
            me_labels.append(flabels)
            me_codes.append(fqcodes)
            me_titles.append('Fluctuating Magnetic Energy Density')    

        labels.append(me_labels)
        qcodes.append(me_codes)
        ptitles.append(me_titles)        


        
    if (velocity):
        for i in range(3):
            energy[:,0]= ga.vals[:,ga.lut[401]]
            energy[:,1]= ga.vals[:,ga.lut[405]]
            energy[:,2]= ga.vals[:,ga.lut[409]]           

        txt=    "Volume-averaged kinetic energy density. Shown are:\n"
        ptitle="Kinetic Energy Density vs. Time"
    elif (magnetic):
        for i in range(3):
            energy[:,0]= ga.vals[:,ga.lut[1101]]
            energy[:,1]= ga.vals[:,ga.lut[1105]]
            energy[:,2]= ga.vals[:,ga.lut[1109]]            

        txt=    "Volume-averaged magnetic energy density. Shown are:\n"
        ptitle= "Magnetic Energy Density vs Time"
    ptitle=''
    txt=''
    txt=txt+"  (1) Total energy density.\n"
    txt=txt+"  (2) Mean energy density ( $m=0$ only).\n"
    txt=txt+"  (3) Fluctuating energy density ( $m\\ne0$ only).\n"

    for k,qtemp in enumerate(qcodes):   # Iterate over Kinetic and Magnetic energy densities
        

        nplot = len(qtemp)
        #fig, ax = plt.subplots(nrows=nplot,ncols=1,figsize=(8.5,11))
        
        fig = plt.figure(figsize=(pxsize,pysize))
        print(qtemp, nplot)
        ysize = (1.0-2*ymar-(nplot-1)*yspace)/nplot
        xsize = 1.0-2*xmar

        for j, qcode in enumerate(qtemp):
            #fig = plt.figure(figsize=(8.5,11))
            yone = 1-ymar-(j+1)*ysize-j*yspace
            ax1 = fig.add_axes((xmar, yone, xsize,ysize))
            #    ax1 = ax
            #else:
            #    ax1 = ax[j]
            ax1.set_title(ptitle)
            ax1.set_xlabel('Time')
            ax1.set_ylabel('Energy Density')
            ax1.set_title(ptitles[k][j])
            
            maxvals = []
            minvals = []
            for i in range(len(qcode)):
                ax1.plot(time, ga.vals[:,ga.lut[qcode[i]]], label=labels[k][j][i])
                maxvals.append(numpy.max(ga.vals[:,ga.lut[qcode[i]]]))   
                minvals.append(numpy.min(ga.vals[:,ga.lut[qcode[i]]]))
                       
            ax1.set_yscale('log')
            ax1.legend()

            # Attemp to do some sensible scaling
            
            ymax = numpy.max(maxvals)
            ymin = numpy.min(minvals)

            ylogtop = int(numpy.log10(ymax))+1
            ybottemp = int(numpy.log10(ymin))-1

            if (decades > 0):
                ylogbot = ylogtop-decades
                
                if (ylogbot < ybottemp):
                    ylogbot = ybottemp
            else:
                ylogbot = ybottemp

            ax1.set_ylim([10**ylogbot,10**ylogtop])

            #fig.legend(bbox_to_anchor=(0.9,0.2))
            #plt.tight_layout()
        if (pdf != None):
            pdf.savefig()
            plt.close()

            # make the edge colors match the facecolors
            # center text

            #fig.text(.1, .2, txt) #, ha='center')


            # resize the figure to match the aspect ratio of the Axes    
            #fig.set_size_inches(7, 8, forward=True)

            #plt.savefig(pfile)
        #plot_energy_v_time('gatest.dat','gplot.pdf')

