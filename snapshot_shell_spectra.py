def plot_spectra(infile, velocity=True,magnetic=False,thermal=False, decades=10, pdf = None, time_avg=True, rind=None, fontsize=12):
    # Provide full path to infile
    from matplotlib import pyplot as plt
    from matplotlib import rcParams
    import numpy as np
    import numpy
    from rayleigh_diagnostics import Shell_Spectra

    if (magnetic):
        velocity = False
    if (thermal):
        velocity = False

    spec = Shell_Spectra(infile,path='',time_avg=time_avg)
    if (rind == None):
        rind = 0 #spec.nr//2+spec.nr%2
        #print(spec.nr, rind, spec.lpower.shape)
    nell = spec.nell
    ellp1=np.arange(1,nell+1,dtype='int32')   
    power = numpy.zeros((nell,3),dtype='float64')
    rcParams.update({'font.size': fontsize})



    vvec = [r'V$_r$', r'V$_\theta$', r'V$_\phi$']
    vtxt = 'Time-averaged velocity power spectrum at radius r = '+"{:.4e}".format(spec.radius[rind])+".\n\nShown are:\n"
    vtitle = 'Velocity Power Spectrum'
    vlabels = [r' $P_\ell(r)\equiv\,\,\sum_{m=0}^{\ell}\left(|v_{\ell,r}^{m}(r)|^2+|v_{\ell,\theta}^{m}(r)|^2+|v_{\ell,\phi}^{m}(r)|^2  \right)$', 
              r' $\overline{P}_\ell(r)\equiv\,\,|v_{\ell,r}^{0}(r)|^2+|v_{\ell,\theta}^{0}(r)|^2+|v_{\ell,\phi}^{0}(r)|^2 $  ', 
              r" $P^{\prime}_\ell\equiv\,\,P_\ell(r)-\overline{P}_\ell(r)$"]


    bvec = [r'B$_r$', r'B$_\theta$', r'B$_\phi$']
    btxt = 'Time-averaged velocity power spectrum at radius r = '+"{:.4e}".format(spec.radius[rind])+".\n\nShown are:\n"
    btitle= 'Magnetic Power Spectrum'
    blabels = [r' $P_\ell(r)\equiv\,\,\sum_{m=0}^{\ell}\left(|B_{\ell,r}^{m}(r)|^2+|B_{\ell,\theta}^{m}(r)|^2+|B_{\ell,\phi}^{m}(r)|^2  \right)$', 
              r' $\overline{P}_\ell(r)\equiv\,\,|B_{\ell,r}^{0}(r)|^2+|B_{\ell,\theta}^{0}(r)|^2+|B_{\ell,\phi}^{0}(r)|^2 $  ', 
              r" $P^{\prime}_\ell\equiv\,\,P_\ell(r)-\overline{P}_\ell(r)$"]    
    
    
    vec = [vvec, bvec]
    txts = [vtxt, btxt]
    titles = [vtitle, btitle]
    labels = [vlabels, blabels]

    pxsize = 8.5
    pysize = 11
    ymar = 1.0/pysize
    xmar = 1.0/pxsize
    xpsize = 1-2*xmar
    yone = 0.6
    ypsize = 0.3


    nplt = 1
    if (magnetic):    
        nplt=2
        
        
    for k in range(nplt):    
        if (k==0):
            for i in range(3):
                vr = spec.lpower[:,rind,spec.lut[1],0,i]
                vt = spec.lpower[:,rind,spec.lut[2],0,i]
                vp = spec.lpower[:,rind,spec.lut[3],0,i]
            
                power[:,i] = vr+vt+vp

        elif (k==1):
            for i in range(3):    
                br = spec.lpower[:,rind,spec.lut[801],0,i]
                bt = spec.lpower[:,rind,spec.lut[802],0,i]
                bp = spec.lpower[:,rind,spec.lut[803],0,i]
            
                power[:,i] = br+bt+bp
            txt=    "Time-averaged magnetic power spectrum at radius r = "+"{:.4e}".format(spec.radius[rind])+". Shown are:\n"            

        ptitle = titles[k]+' ( r = '+"{:.4e}".format(spec.radius[rind])+')'
        txt=txts[k]+"  (1) Total power, integrated over all $m$-values.\n"
        txt=txt+"  (2) Axisymmetric power, associated with $m=0$ only.\n"
        txt=txt+"  (3) Convective Power, integrated over all $m\\ne0$.\n"
        txt=txt+"\n"
        #txt=txt+"Radius:  "+"{:.4e}".format(spec.radius[rind])
        #txt=txt+"\n\n"
        txt=txt+"Time-averaging information:\n\n"
        txt=txt+"  Beginning timestep: "+"{:,}".format(spec.iters[0])+'\n'
        txt=txt+"  Ending timestep   : "+"{:,}".format(spec.iters[1])+'\n\n'
        txt=txt+"  Beginning time: "+"{:.4e}".format(spec.time[0])+'\n'
        txt=txt+"  Ending time   : "+"{:.4e}".format(spec.time[1])+'\n\n'    
        txt=txt+"  Delta time    : "+"{:.4e}".format(spec.time[1]-spec.time[0])    
        
        Pdef=r"P$_\ell\,\,\equiv\,\,\sum_{m=0}^{m=\ell}\left(v_r  \right)"

        fig = plt.figure(figsize=(pxsize,pysize))
        
        
        #fig = plt.figure()
        ax1 = fig.add_axes((xmar, yone, xpsize, ypsize))

        ax1.set_title(ptitle)
        ax1.set_xlabel('Degree $\ell$+1')
        ax1.set_ylabel('Power ')
        for i in range(3):
            ax1.plot(ellp1, power[:,i], label=labels[k][i])
        ax1.set_yscale('log')
        ax1.set_xscale('log')

        ymax = numpy.max(power)
        ylogtop = int(numpy.log10(ymax))+1
        ylogbot = ylogtop-decades
        #print(ymax,ylogtop,ylogbot)
        ax1.set_ylim([10**ylogbot,10**ylogtop])

        rcParams.update({'font.size': 16})        
        fig.legend(bbox_to_anchor=(xmar,0.55),loc='upper left')
        rcParams.update({'font.size': fontsize})
        fig.text(xmar, ymar, txt) #, ha='center')
        #fig.text(xmar, 0.2, Pdef)


        # resize the figure to match the aspect ratio of the Axes    
        #fig.set_size_inches(7, 8, forward=True)
        if (pdf != None):
            pdf.savefig()
            plt.close()
    #plt.savefig(pfile)

