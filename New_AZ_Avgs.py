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


class AZ_Avgs:
    """Rayleigh Shell-Averaged Data Structure
    ----------------------------------

    self.niter                         : number of time steps
    self.nq                            : number of diagnostic quantities output
    self.nr                            : number of radial points
    self.qv[0:nq-1]                    : quantity codes for the diagnostics output
    self.radius[0:nr-1]                : radial grid

    For version 1 (earliest Rayleigh versions):
    self.vals[0:nr-1,0:nq-1,0:niter-1] : The spherically averaged diagnostics
                                             

    For versions > 1:
    self.vals[0:n-1,0:3,0:nq-1,0:niter-1] : The spherically averaged diagnostics
                                            0-3 refers to moments taken over spherical shells
                                            (index 0 is mean, index 3 is kurtosis)    
                                            
    self.iters[0:niter-1]              : The time step numbers stored in this output file
    self.time[0:niter-1]               : The simulation time corresponding to each time step
    self.version                       : The version code for this particular output (internal use)
    self.lut     

    self.time_averaged                 : If True, data has been time averaged.
                                         In this case, iters and time have dimension 2, and contain the
                                         initial and final times and iterations of the averaged data.

    Initialization Examples:
        (1):  Read in a single AZ_Avgs file
              a = AZ_Avgs('00000001',path='./AZ_Avgs/')
              
        (2):  Time-average data from multiple AZ_Avgs files:
              a = AZ_Avgs(['00000001','00000002'])
              
        (2):  Concatenate time-series data from multiple AZ_Avgs files:
              a = AZ_Avgs(['00000001','00000002'],time_average=False)

        (3):  Time-averaged data from multiple AZ_Avgs files + save data to new AZ_Avgs file:
              a = AZ_Avgs(['00000001','00000002'],ofile='my_save_file.dat')
              
        (4):  Concatenate time-series data from multiple AZ_Avgs files.
              Extract only quantity codes 401 and 402:
              a = AZ_Avgs(['00000001','00000002'], qcodes = [401,402], time_average = False)

    Additional Notes:
        For concatenated or averaged data:
        
            - Output files are written in identical format to a standard AZ_Avgs file. They may be
              read in using the same syntax described above.
              
            - To save concatenated or averaged data to a new file after initialization:
              a.write('filename')  [ where 'a' is the AZ_Avgs object name]
              
    """

    def __init__(self,filename='none',path='AZ_Avgs/', ofile=None, qcodes=[],time_average=True):

        """
           Initializer for the AZ_Avgs class.
           
           Input parameters:
               filename     : (a) {string} The AZ_Avgs file to read.
                            : (b) {list of strings} The AZ_Avgs files whose time-series
                               data you wish to concatenate or time-average.
               path         : The directory where the file is located (if full path not in filename)
               qcodes       : {optional; list of ints} Quantity codes you wish to extract (default is to extract all)
               ofile        : {optional; string}; default = None Filename to save time-averaged data to, if desired.
               time_average : {optional; Boolean; default = True} Time-average data from multiple files if a list is provided.
                              If a list of files is provided and this flag is set to False, data will be concatenated instead. 
        """
        # Check to see if we are compiling data from multiple files
        # This is true if (a) ofile is set or (b) filename is a list of files, rather than a single string
        # if (a) is false, but (b) is true, no output is written, but compiled output is returned        
        
        #(a)
        multiple_files = False
        mainfile=filename
        if (ofile != None):
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
        self.read_data(qcodes=qcodes, ntheta=ntheta)
        
        if (multiple_files):
            if (not time_average):
                self.compile_multiple_files(filename, qcodes = qcodes,path=path,ntheta=ntheta)
            else:
                self.time_average_files(filename, qcodes = qcodes,path=path, ntheta=ntheta)
                
        if (ofile != None):
            self.write(ofile)
            
            
    def write(self,outfile):
        """
            Writing routine for AZ_Avgs class.
            
            Input parameters:
               outfile  : (a) {string} The AZ_Avgs file to be written.

            Calling Syntax:
                my_shellavgs.write("my_file.dat")
                (writes all data contained in my_shellavgs to my_file.dat, in standard AZ_Avgs format)

        """

        numdim = 6

        one_rec = np.dtype([('vals', np.float64, [self.nq, self.nr, self.ntheta]), ('times',np.float64), ('iters', np.int32)  ])   

        fstruct = np.dtype([ ('fdims', np.int32, numdim), ('qvals', np.int32,(self.nq)), ('radius',np.float64,(self.nr)), 
                             ('costheta',np.float64,(self.ntheta)), ('fdata', one_rec, [self.niter,])  ])
        
                        
        odata = np.zeros((1,),dtype=fstruct)
        odata['fdims'][0,0]=314
        odata[0]['fdims'][1]=self.version
        odata[0]['fdims'][2]=self.niter
        odata[0]['fdims'][3]=self.nr
        odata[0]['fdims'][4]=self.ntheta
        odata[0]['fdims'][5]=self.nq

	                
        odata[0]['qvals'][:]=self.qv[:]
        odata[0]['radius'][:]=self.radius[:]
        odata[0]['costheta'][:]=self.costheta[:]
        odata[0]['fdata']['times'][:]=self.time[0:self.niter]
        odata[0]['fdata']['iters'][:]=self.iters[0:self.niter]
        
        odata[0]['fdata']['vals'][:,:,:,:]=np.transpose(self.vals[:,:,:,:])
            
        fd = open(outfile,'wb')
        odata.tofile(fd)
        if (self.time_averaged):
            self.time[1].tofile(fd)
            self.iters[1].tofile(fd)
        fd.close()
        

    def compile_multiple_files(self,filelist,qcodes=[],path=''):
        """
           Time-series concatenation routine for the AZ_Avgs class.
           
           Input parameters:
               filelist  : {list of strings} The AZ_Avgs files to be concatenated.
               path      : The directory where the files are located (if full path not in filename)
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)
               ntheta    : {optional; int; default = 0} Set this value to correct the variance in 
                           version=2 AZ_Avgs (only mean is preserved otherwise for that version).
           Notes:
                - This routine is incompatibile with version=1 AZ_Avgs due to lack of 
                  moments output in that original version.  All other versions are compatible.
                  
                - This routine assumes radial resolution does not change across files in filelist.

        """

            
        nfiles = len(filelist)
        new_nrec = self.niter*nfiles
        self.niter = new_nrec
        self.vals = np.zeros((self.nr,4,self.nq,self.niter),dtype='float64')
        self.iters = np.zeros(self.niter          ,dtype='int32')
        self.time  = np.zeros(self.niter          ,dtype='float64')        
        
        k = 0
        for i in range(nfiles):

            a = AZ_Avgs(filelist[i],qcodes=qcodes,path=path)

            #If the iteration count isn't constant throughout a run,
            #we may need to resize the arrays
            if (k+a.niter > self.niter):
                self.niter = k+a.niter
                self.time.resize((self.niter))
                self.iters.resize( (self.niter))
                # Note that using numpy's resize routine gets a little tricky here
                # due to the striping of the vals array.  Handle this 'manually'
                # for now
                vals = np.zeros((self.ntheta,self.nr, self.nq,self.niter),dtype='float64')
                vals[:,:,:,0:k] = self.vals[:,:,:,0:k]
                self.vals = vals

            self.time[k:k+a.niter]   = a.time[:]
            self.iters[k:k+a.niter]  = a.iters[:]
            self.vals[:,:,:,k:k+a.niter] = a.vals[:,:,:,:]
            k+=a.niter
            
        # Trim off any excess zeros 
        # (in case last file was incomplete)
        self.niter = k
        self.time = self.time[0:self.niter]
        self.iters = self.iters[0:self.niter]
        self.vals = self.vals[:,:,:,0:self.niter]


    def time_average_files(self,filelist,qcodes=[],path=''):
        """
           Time-series concatenation routine for the AZ_Avgs class.
           
           Input parameters:
               filelist  : {list of strings} The AZ_Avgs files to be time-averaged.
               path      : The directory where the files are located (if full path not in filename)
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)
               ntheta    : {optional; int; default = 0} Set this value to correct the variance in 
                           version=2 AZ_Avgs (only mean is preserved otherwise for that version).               

           Notes:
                - This routine is incompatibile with version=1 AZ_Avgs due to lack of 
                  moments output in that original version.  All other versions are compatible.
                  
                - This routine assumes radial resolution does not change across files in filelist.

        """

            
        nfiles = len(filelist)

        self.niter = 1
        self.vals = np.zeros((self.ntheta,self.nr,self.nq,1),dtype='float64')
        self.iters = np.zeros(2          ,dtype='int32')
        self.time  = np.zeros(2          ,dtype='float64')        
        
        #Integrate in time using a midpoint method.
        last_dt = 0.0
        total_time = 0.0
        a = AZ_Avgs(filelist[0],qcodes=qcodes,path=path)
        for i in range(nfiles):
            weights = np.zeros(a.niter,dtype='float64')
            
            if (i != nfiles-1):
                b = AZ_Avgs(filelist[i+1],qcodes=qcodes,path=path)
                
            if (i == 0):
                self.iters[0] = a.iters[0]
                self.time[0] = a.time[0]
                
            for j in range(a.niter-1):
                weights[j] = 0.5*(last_dt + a.time[j+1]-a.time[j])
                
            if (i != nfiles-1):
                weights[a.niter-1] = 0.5*(last_dt+b.time[0]-a.time[a.niter-1])
            else:
                weights[a.niter-1] = 0.5*(a.time[a.niter-1]-a.time[a.niter-2])

            total_time = total_time+np.sum(weights)
            for j in range(a.niter):
                self.vals[:,:,:,0]+=a.vals[:,:,:,j]*weights[j]
            last_dt = a.time[a.niter-1]-a.time[a.niter-2]
            if (i != nfiles-1):
                a = b

        self.iters[1] = a.iters[a.niter-1]
        self.time[1]  = a.time[a.niter-1]
        self.vals = self.vals/total_time
        self.version=a.version+100  # 100+ version numbers indicate time-averaging
        self.time_averaged = True
        
        
    def read_data(self,qcodes = []):
        """
           Data-reading routine for the AZ_Avgs class.
           
           Input parameters:
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)
               ntheta    : {optional; int; default = 0} Set this value to correct the variance in 
                           version=2 AZ_Avgs (only mean is preserved otherwise for that version).
           
           Notes:
                - This routine does not read header information.  Read_dimensions must be called first
                  so that dimentions and the file descriptor are initialized. 
        """
    
        self.qv    = np.zeros(self.nq             ,dtype='int32')
        self.vals  = np.zeros((self.ntheta, self.nr, self.nq, self.niter),dtype='float64')
        self.iters = np.zeros(self.niter+self.time_averaged          ,dtype='int32')
        self.time  = np.zeros(self.niter+self.time_averaged          ,dtype='float64')
   
            
        self.radius = np.zeros(self.nr, dtype='float64')
        self.costheta = np.zeros(self.ntheta, dtype='float64')

        # Set up the record structure

        one_rec = np.dtype([('vals', np.float64, [self.nq, self.nr, self.ntheta]), ('times',np.float64), ('iters', np.int32)  ])   
            

        fstruct = np.dtype([ ('qvals', np.int32,(self.nq)), ('radius',np.float64,(self.nr)), 
                             ('costheta',np.float64,(self.ntheta)), ('fdata', one_rec, [self.niter,])  ])

        
        if (self.byteswap):
                fdata = np.fromfile(self.fd,dtype=fstruct,count=1).byteswap()
                if (self.time_averaged):
                    self.time[1] = np.fromfile(self.fd,dtype='float64',count=1).byteswap() 
                    self.iters[1] = np.fromfile(self.fd,dtype='int32',count=1).byteswap()    
        else:
                fdata = np.fromfile(self.fd,dtype=fstruct)
                if (self.time_averaged):
                    self.time[1] = np.fromfile(self.fd,dtype='float64',count=1) 
                    self.iters[1] = np.fromfile(self.fd,dtype='int32',count=1)                  

        self.fd.close()
        
        self.time[:] = fdata['fdata']['times'][0,:self.niter]
        self.iters[:] = fdata['fdata']['iters'][0,:self.niter]
        self.qv[:] = fdata['qvals'][0,:]
        self.radius[:] = fdata['radius'][0,:]
        self.costheta[:] = fdata['costheta'][0,:]
        


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
            

            self.vals  = np.zeros((self.ntheta,self.nr, self.nq, self.niter),dtype='float64')            
                
            for q in range(self.nq):
                qcheck = self.lut[qget[q]]
                if (qcheck < maxq):
                        self.vals[:,:,q,:] = vals[:,:,qcheck,:]
            self.lut = get_lut(self.qv)  # Rebuild the lookup table since qv has changed
            
    def read_dimensions(self,the_file,closefile=False):
        """
        
            Header-reading routine for the AZ_Avgs class.
            
            This routine initialized the file descriptor and reads the dimensions 
            and Endian information of a AZ_Avgs file only. It does not read the AZ_Avgs data itself.

           
            Input parameters:
               the_file    : {string} AZ_Avgs file whose header is to be read.
               closefile   : (Boolean; optional; default=False) Set to True to close the file 
                             after reading the header information. 
                   

        """
        print('the_file is: ', the_file)
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
        self.ntheta  = specs[4]
        self.nq      = specs[5]
        
        self.time_averaged = (self.version > 100)
        if (closefile):
            self.fd.close()
   
