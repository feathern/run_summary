#plot_spectra(sp_ofile,velocity=True,magnetic=prefs.magnetic,thermal=False, decades=10, pdf = pdf)
def plot_shell_slice(shell_file,pdf,magnetic=False):
    from rayleigh_diagnostics import Shell_Slices
    from projection import plot_ortho
    from matplotlib import gridspec
    import numpy
    import matplotlib.pyplot as plt


    pdf_file = 'two_rows_three_columns.pdf'

    pltin = 3.5  # Size of each subimage in inches (will be square)

    nrow=3
    ncol=2
    fig = plt.figure(constrained_layout=False, figsize=(8.5,11) ) # figsize=(pltin*ncol,pltin*nrow*1.1))
    spec = gridspec.GridSpec(ncols=ncol, nrows=nrow*2, figure=fig, height_ratios=[1,.1]*nrow, width_ratios=[1]*ncol)


    plt.rcParams.update({'font.size': 12})

    qcodes = [1,2, 3]   # V_R and V_phi and Temperature/Entropy (Boussinesq/Anelastic)
    names = [r'v$_r$',r'v$_\theta$', r'v$_\phi$']
    if (magnetic):
        qcodes = qcodes+[801,802,803]
        names = names+[r'B$_r$',r'B$_\theta$', r'B$_\phi$']
    subtract_mean = [False, False, True]  # Temperature/Entropy renders better with mean subtracted
                                          # (This is because the color table runs from - to + by default.)
    
    lv = 0   # level to plot (uppermost in file)
    ts = 0   # timestep to plot (first timestep in file)
    
    gwidth=1 # width of grid lines for each image
    style = '-'  # line style for the grid lines
    # Width of the horizon line or each image (Default: 2)
    hwidth = 2

    # A color table for each image (Default: RdYlBu_r)
    cmap = "RdYlBu_r"

    latcen = 45 # Latitudinal center (degrees) of vantage point (Default: 45 N)
    loncen = 0 # Longitudinal center (degrees) of vantage point (Default: 0)



    ##########################################################
    # If the grid is plotted, the number of latitude lines
    # for the grid can be controlled via the nlats keyword.
    # Default: 9
    # Note that if nlats is even, the equator will not be drawn
    nlat = 9


    ##############################################################################
    # Similarly, the nlons keyword can be used to control longitude lines
    # More precisely, it controls the number of MERIDIANS (great circles) drawn
    # Default:  8

    nlon = 8

    #Longitude grid-lines can be drawn in one of two ways:
    # 1)  Completely to the pole (polar_style = 'polar')
    # 2)  Truncated at the last line of latitue drawn (polar_style = 'truncated')
    # Default:  "truncated"
    pstyle = 'truncated'


    ##############################################################
    # We can also control the way in which the image is saturated
    # via the scale_type keyword.   There are three possibilities:
    # 1) scale_type=['rms', a], where a is of type 'float'
    #    In this instance, the image bounds are -a*rms(data), +a*rms(data)
    # 2) scale_type = ['abs', a]
    #    In this instance, the image bounds are -a*abs(data), +a*abs(data)
    # 3) scale_type= ['force', [a,b]]
    #    In this instance, the image bounds are a,b
    # 4) scale_type = [None,None]
    #    In this instance, the image bounds are min(image), max(image)
    # Default:  [None,None]
    # Note that rms and abs are taken over projected image values, not input data
    # (you only see half the data in the image)
    scale_type = ['rms',2.0]


    # Number of pixels across each projected, interpolated image
    # 768 is the default and seems to do a reasonable job
    nyzi = 768

    s1=Shell_Slices(shell_file,path='')
    data = numpy.zeros((s1.nphi,s1.ntheta),dtype='float64')
    costheta = s1.costheta

    counter = 0
    for i in range(nrow):
        for j in range(ncol):
            if (counter < len(qcodes)):
                qcode = qcodes[counter]
                name = names[counter]
                row_ind = 2*i  # The color bars are actually on a row of their own.  Skip those rows with 2*i 
                col_ind = j

            
                data[:,:] = s1.vals[:,:,lv,s1.lut[qcode],ts]
                print('Plotting: ', name)   
                if (subtract_mean[j]):
                    data = data-numpy.mean(data)


                ax    = fig.add_subplot(spec[row_ind,col_ind])
                cspec = spec[row_ind+1,col_ind]
                caxis=None

                plot_ortho(data,s1.costheta,fig,ax,caxis, hwidth=hwidth, gridstyle=style, 
                           gridwidth=gwidth, nyz=nyzi, colormap=cmap, 
                           latcen=latcen,
                           pole_style=pstyle, nlats = nlat,scale_type=scale_type)
                ind = s1.inds[lv]    # The radial index within the full Rayleigh grid
                rval = s1.radius[lv] # The value of that radius
                ptitle = name
                if (j == 0 and i == 0):  # only print radius information in first column
                    ptitle=ptitle+"    (radius = "+"{:2e}".format(rval) +", index ="+str(ind)+")"
                ax.set_title(ptitle)
                counter+=1




    pdf.savefig()
    plt.close()

