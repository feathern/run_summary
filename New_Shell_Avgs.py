import numpy as np

maxq = 4000

def get_lut(quantities):
    """return the lookup table based on the quantity codes"""
    nq = len(quantities)
    lut = np.zeros(maxq) + maxq
    for i,q in enumerate(quantities):
        if ((0 <= q) and ( q <= maxq-1)): # quantity must be in [0, maxq-1]
            lut[q] = i
    return lut.astype('int32')


class Shell_Avgs:
    """Rayleigh Shell-Averaged Data Structure
    ----------------------------------

    self.niter                         : number of time steps
    self.nq                            : number of diagnostic quantities output
    self.nr                            : number of radial points
    self.qv[0:nq-1]                    : quantity codes for the diagnostics output
    self.radius[0:nr-1]                : radial grid

    For version 1 (earliest Rayleigh versions):
    self.vals[0:nr-1,0:nq-1,0:niter-1] : The spherically averaged diagnostics
                                             

    For versions > 2:
    self.vals[0:n-1,0:3,0:nq-1,0:niter-1] : The spherically averaged diagnostics
                                             0-3 refers to moments (index 0 is mean, index 3 is kurtosis)    
    self.iters[0:niter-1]              : The time step numbers stored in this output file
    self.time[0:niter-1]               : The simulation time corresponding to each time step
    self.version                       : The version code for this particular output (internal use)
    self.lut     


    Initialization Examples:
        (1):  Read in a single Shell_Avgs file
              a = Shell_Avgs('00000001',path='./Shell_Avgs/')
              
        (2):  Concatenate time-series data from multiple Shell_Avgs files:
              a = Shell_Spectra(['00000001','00000002'])

        (3):  Concatenate time-series data from multiple Shell_Avgs files + save data to new Shell_Avgs file:
              a = Shell_Spectra(['00000001','00000002'],ofile='my_save_file.dat')
              
        (4):  Concatenate time-series data from multiple Shell_Avgs files.
              Extract only quantity codes 401 and 402:
              a = Shell_Spectra(['00000001','00000002'], qcodes = [401,402])

    Additional Notes:
        For concatenated data:
        
            - Output files are written in identical format to a standard Shell_Avgs file. They may be
              read in using the same syntax described above.
              
    """

    def __init__(self,filename='none',path='Shell_Avgs/', ofile='none', qcodes=[],concatenate=False):

        """
           Initializer for the Shell_Avgs class.
           
           Input parameters:
               filename  : (a) {string} The Shell_Avgs file to read.
                         : (b) {list of strings} The Shell_Avgs files whose time-series
                               data you wish to concatenate or time-average.
               path      : The directory where the file is located (if full path not in filename)
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)
               ofile     : {optional; string} Filename to save time-averaged data to, if desired.
               concatenate : {optional; Boolean} Concatenate data from multiple files if a list is provided.
                             If a list of files is provided and this flag is not set, data will be time-averaged. 
        """
        # Check to see if we are compiling data from multiple files
        # This is true if (a) ofile is set or (b) filename is a list of files, rather than a single string
        # if (a) is false, but (b) is true, no output is written, but compiled output is returned        
        
        #(a)
        multiple_files = False
        mainfile=filename
        if (ofile != 'none'):
            multiple_files = True
            
        #(b)
        ltype = type([])        
        ftype=type(filename)
        if (ltype == ftype):
            multiple_files = True
        else:
            if (mainfile == 'none'):
                mainfile = path+'00000001'
            else:
                mainfile = path+mainfile             
        
        if (multiple_files):
            mainfile = filename[0]
        
        # Read the main file no matter what    
        self.read_dimensions(mainfile)   
        self.read_data(qcodes=qcodes)
        
        if (multiple_files):
            if (concatenate):
                self.compile_multiple_files(filename, qcodes = qcodes,path=path)
            else:
                self.time_average_files(filename, qcodes = qcodes,path=path)
                
        if (ofile != 'none'):
            self.write(ofile)
            
            
    def write(self,outfile):
        """
            Writing routine for Shell_Avgs class.
            
            Input parameters:
               outfile  : (a) {string} The Shell_Avgs file to be written.

            Calling Syntax:
                my_gavgs.write("my_file.dat")
                (writes all data contained in my_gavgs to my_file.dat, in standard Shell_Avgs format)

        """
        one_rec = np.dtype([('vals', np.float64, [self.nq,]), ('times',np.float64), ('iters', np.int32)  ])
        fstruct = np.dtype([ ('fdims', np.int32,4), ('qvals', np.int32,(self.nq)), ('fdata', one_rec, [self.niter,])  ])
        
        odata = np.zeros((1,),dtype=fstruct)
        print(odata['fdims'].shape)
        odata['fdims'][0,0]=314
        odata[0]['fdims'][1]=self.version
        odata[0]['fdims'][2]=self.niter
        odata[0]['fdims'][3]=self.nq
        odata[0]['qvals'][:]=self.qv[:]
        odata[0]['fdata']['times'][:]=self.time
        odata[0]['fdata']['iters'][:]=self.iters
        odata[0]['fdata']['vals'][:,:]=self.vals
        print(odata['qvals'])
        fd = open(outfile,'wb')
        odata.tofile(fd)
        fd.close()
        

    def compile_multiple_files(self,filelist,qcodes=[],path=''):
        """
           Time-series concatenation routine for the Shell_Avgs class.
           
           Input parameters:
               filelist  : {list of strings} The Shell_Avgs files to be concatenated.
               path      : The directory where the files are located (if full path not in filename)
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)

           Notes:
                - This routine is incompatibile with version=1 Shell_Avgs due to lack of 
                  moments output in that original version.  All other versions are compatible.
                  
                - This routine assumes radial resolution does not change across files in filelist.

        """
        if (self.version == 1):
            print('Error: This is a version=1 Shell_Avgs file.  Concatenation only works with version>1 Shell_Avgs.')
            return
            
        nfiles = len(filelist)
        new_nrec = self.niter*nfiles
        self.niter = new_nrec
        self.vals = np.zeros((self.nr,4,self.nq,self.niter),dtype='float64')
        self.iters = np.zeros(self.niter          ,dtype='int32')
        self.time  = np.zeros(self.niter          ,dtype='float64')        
        
        k = 0
        for i in range(nfiles):
            a = Shell_Avgs(filelist[i],qcodes=qcodes,path=path)

            #If the iteration count isn't constant throughout a run,
            #we may need to resize the arrays
            if (k+a.niter > self.niter):
                self.niter = k+a.niter
                self.time.resize((self.niter))
                self.iters.resize( (self.niter))
                # Note that using numpy's resize routine gets a little tricky here
                # due to the striping of the vals array.  Handle this 'manually'
                # for now
                vals = np.zeros((self.nr,4,self.nq,self.niter),dtype='float64')
                vals[:,:,:,0:k] = self.vals[:,:,:,0:k]
                self.vals = vals



            self.time[k:k+a.niter]   = a.time[:]
            self.iters[k:k+a.niter]  = a.iters[:]
            self.vals[:,:,:,k:k+a.niter] = a.vals[:,:,:,:]
            k+=a.niter
            
        # Trim off any excess zeros 
        # (in case last final was incomplete)
        self.niter = k
        self.time = self.time[0:self.niter]
        self.iters = self.iters[0:self.niter]
        self.vals = self.vals[:,:,:,0:self.niter]


    def time_average_files(self,filelist,qcodes=[],path=''):
        """
           Time-series concatenation routine for the Shell_Avgs class.
           
           Input parameters:
               filelist  : {list of strings} The Shell_Avgs files to be time-averaged.
               path      : The directory where the files are located (if full path not in filename)
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)

           Notes:
                - This routine is incompatibile with version=1 Shell_Avgs due to lack of 
                  moments output in that original version.  All other versions are compatible.
                  
                - This routine assumes radial resolution does not change across files in filelist.

        """
        if (self.version == 1):
            print('Error: This is a version=1 Shell_Avgs file.  Time-averaging only works with version>1 Shell_Avgs.')
            return
            
        nfiles = len(filelist)

        self.niter = 1
        self.vals = np.zeros((self.nr,4,self.nq,1),dtype='float64')
        self.iters = np.zeros(2          ,dtype='int32')
        self.time  = np.zeros(2          ,dtype='float64')        
        
        icount = 0
        for i in range(nfiles):
            a = Shell_Avgs(filelist[i],qcodes=qcodes,path=path)
            if (i == 0):
                self.iters[0] = a.iters[0]
                self.time[0] = a.time[0]


            for j in range(a.niter):
                self.vals[:,:,:,0]+=a.vals[:,:,:,j]
            icount+= a.niter  

        self.iters[1] = a.iters[a.niter-1]
        self.time[1]  = a.time[a.niter-1]
        self.vals = self.vals/icount

        
        
    def read_data(self,qcodes = []):
        """
           Data-reading routine for the Shell_Avgs class.
           
           Input parameters:
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)
           
           Notes:
                - This routine does not read header information.  Read_dimensions must be called first
                  so that dimentions and the file descriptor are initialized. 
        """
    
        self.qv    = np.zeros(self.nq             ,dtype='int32')
        self.vals  = np.zeros((self.nr, 4, self.nq, self.niter),dtype='float64')
        self.iters = np.zeros(self.niter          ,dtype='int32')
        self.time  = np.zeros(self.niter          ,dtype='float64')
        self.radius = np.zeros(self.nr, dtype='float64')

        # Set up the record structure
        if (self.version < 6):

            if (self.version ==1):
                one_rec = np.dtype([('vals', np.float64, [self.nq,self.nr]), ('times',np.float64), ('iters', np.int32)  ])
            else:
                one_rec = np.dtype([('vals', np.float64, [self.nq, 4, self.nr]), ('times',np.float64), ('iters', np.int32)  ])   
            
        else:
            # Things are a little more complicated following the parallel I/O redo.
            # Store all data values in a 1-D array and rerrange each record at the end
            one_rec = np.dtype([('vals', np.float64, [self.nr*4*self.nq]), ('times',np.float64), ('iters', np.int32)  ])

        fstruct = np.dtype([ ('qvals', np.int32,(self.nq)), ('radius',np.float64,(self.nr)), ('fdata', one_rec, [self.niter,])  ])

        
        if (self.byteswap):
                fdata = np.fromfile(self.fd,dtype=fstruct,count=1).byteswap()
        else:
                fdata = np.fromfile(self.fd,dtype=fstruct)
                
        self.fd.close()
        
        self.time[:] = fdata['fdata']['times'][0,:]
        self.iters[:] = fdata['fdata']['iters'][0,:]
        self.qv[:] = fdata['qvals'][0,:]
        self.radius[:] = fdata['radius'][0,:]
        
        ####################################################3
        # Reading in the values is a little more complicated since the Shell_Avgs
        # data structure has undergone a few different iterations.
        if (self.version >= 6):
            vals = np.zeros((self.nr,4,self.nq,self.niter),dtype='float64')
            for i in range(self.niter):
                rind=0
                kone = 0
                nr_base = self.nr//self.npcol
                nr_mod = self.nr % self.npcol
                for j in range(self.npcol):
                    nrout= nr_base
                    if (j < nr_mod) :
                        nrout=nrout+1
                    ktwo = kone+nrout*4*self.nq
                    tmp = np.reshape(fdata['fdata']['vals'][0,i,kone:ktwo], (nrout,4,self.nq), order = 'F'   )
                    vals[rind:rind+nrout,:,:,i] = tmp[:,:,:]
                    rind=rind+nrout            
                    kone = ktwo
        elif (self.version == 1):
            #vals = numpy.zeros((self.nr,self.nq,self.niter)
            vals = np.transpose(fdata['fdata']['vals'][0,:,:,:])

        else:
            #vals = numpy.zeros((self.niter,self.nq,4,self.nr))
            vals = np.transpose(fdata['fdata']['vals'][0,:,:,:,:])

        self.lut = get_lut(self.qv)  # Lookup table
        
        ########################################################
        # We may want to extract a subset of the quantity codes
        if (len(qcodes) == 0):
            self.vals = vals
        else:
            nqfile = self.nq        # number of quantity codes in the file
            qget = np.array(qcodes,dtype='int32')
            self.qv = qget  # Don't update the lookup table yet
            self.nq = len(self.qv)  # number of quantity codes we will extract
            
            if (self.version == 1):
                self.vals  = np.zeros((self.nr,self.nq, self.niter),dtype='float64')
            else:
                self.vals  = np.zeros((self.nr,4, self.nq, self.niter),dtype='float64')            
                
            for q in range(self.nq):
                qcheck = self.lut[qget[q]]
                if (qcheck < maxq):
                    if (self.version == 1):
                        self.vals[:,q,:] = vals[:,qcheck,:]
                    else:
                        self.vals[:,:,q,:] = vals[:,:,qcheck,:]
            self.lut = get_lut(self.qv)  # Rebuild the lookup table since qv has changed


            
    def read_dimensions(self,the_file,closefile=False):
        """
        
            Header-reading routine for the Shell_Avgs class.
            
            This routine initialized the file descriptor and reads the dimensions 
            and Endian information of a Shell_Avgs file only. It does not read the Shell_Avgs data itself.

           
            Input parameters:
               the_file    : {string} Shell_Avgs file whose header is to be read.
               closefile   : (Boolean; optional; default=False) Set to True to close the file 
                             after reading the header information. 
                   

        """
        self.fd = open(the_file,'rb')        
        specs = np.fromfile(self.fd,dtype='int32',count=6)
        bcheck = specs[0]       # If not 314, we need to swap the bytes
        self.byteswap = False
        if (bcheck != 314):
            specs.byteswap()
            self.byteswap = True
            
        self.version = specs[1]
        self.niter   = specs[2]
        self.nr      = specs[3]
        self.nq      = specs[4]
        if (self.version >= 6):
            self.npcol = specs[5]
        else:
            #versions < 6 do not have npcol stored.
            #rewind by 4 bytes
            self.fd.seek(-4,1)
        
        
        if (closefile):
            self.fd.close()
   
