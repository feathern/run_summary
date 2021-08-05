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


class G_Avgs:
    """Rayleigh GlobalAverage Structure
    ----------------------------------
    self.niter                  : number of time steps
    self.nq                     : number of diagnostic quantities output
    self.qv[0:nq-1]             : quantity codes for the diagnostics output
    self.vals[0:niter-1,0:nq-1] : The globally averaged diagnostics 
    self.iters[0:niter-1]       : The time step numbers stored in this output file
    self.time[0:niter-1]        : The simulation time corresponding to each time step
    self.version                : The version code for this particular output (internal use)
    self.lut                    : Lookup table for the different diagnostics output
    """

    def __init__(self,filename='none',path='G_Avgs/', ofile='none', qcodes=[]):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename)
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
            self.compile_multiple_files(filename, qcodes = qcodes,path=path)

        if (ofile != 'none'):
            self.write(ofile)
    def write(self,outfile):
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
        nfiles = len(filelist)
        new_nrec = self.niter*nfiles
        self.niter = new_nrec
        self.vals = np.zeros((self.niter,self.nq),dtype='float64')
        self.iters = np.zeros(self.niter          ,dtype='int32')
        self.time  = np.zeros(self.niter          ,dtype='float64')        
        
        k = 0
        for i in range(nfiles):
            a = G_Avgs(filelist[i],qcodes=qcodes,path=path)
            self.time[k:k+a.niter]   = a.time[:]
            self.iters[k:k+a.niter]  = a.iters[:]
            self.vals[k:k+a.niter,:] = a.vals[:]
            k+=a.niter
        
        
        
    def read_data(self,qcodes = []):
    
        self.qv    = np.zeros(self.nq             ,dtype='int32')
        self.vals  = np.zeros((self.niter,self.nq),dtype='float64')
        self.iters = np.zeros(self.niter          ,dtype='int32')
        self.time  = np.zeros(self.niter          ,dtype='float64')

        
        one_rec = np.dtype([('vals', np.float64, [self.nq,]), ('times',np.float64), ('iters', np.int32)  ])
        fstruct = np.dtype([ ('qvals', np.int32,(self.nq)), ('fdata', one_rec, [self.niter,])  ])
        
        if (self.byteswap):
                fdata = np.fromfile(self.fd,dtype=fstruct,count=1).byteswap()
        else:
                fdata = np.fromfile(self.fd,dtype=fstruct)
                
        self.time[:] = fdata['fdata']['times'][0,:]
        self.iters[:] = fdata['fdata']['iters'][0,:]
        self.qv[:] = fdata['qvals'][0,:]
   
        self.fd.close()
        
        self.lut = get_lut(self.qv)  # Lookup table
        
        ########################################################
        # We may want to extract a subset of the quantity codes
        if (len(qcodes) == 0):
            self.vals[:,:] = fdata['fdata']['vals'][0,:,:]
        else:
            nqfile = self.nq        # number of quantity codes in the file
            qget = np.array(qcodes,dtype='int32')
            self.qv = qget  # Don't update the lookup table yet
            self.nq = len(self.qv)  # number of quantity codes we will extract
            self.vals  = np.zeros((self.niter,self.nq),dtype='float64')
            for q in range(self.nq):
                qcheck = self.lut[qget[q]]
                if (qcheck < maxq):
                    self.vals[:,q] = fdata['fdata']['vals'][0,:,qcheck]
            self.lut = get_lut(self.qv)  # Rebuild the lookup table since qv has changed
            
    def read_dimensions(self,the_file,closefile=False):
        "Reads the quantity-code count, niter, version, and Endian information"
        self.fd = open(the_file,'rb')        
        specs = np.fromfile(self.fd,dtype='int32',count=4)
        bcheck = specs[0]       # If not 314, we need to swap the bytes
        self.byteswap = False
        if (bcheck != 314):
            specs.byteswap()
            self.byteswap = True
            
        self.version = specs[1]
        self.niter = specs[2]
        self.nq = specs[3]
        if (closefile):
            self.fd.close()
   
        
#dt = np.dtype([('time', [('min', np.int64), ('sec', np.int64)]),
#
#               ('temp', float)])
#
#x = np.zeros((1,), dtype=dt)
#
#x['time']['min'] = 10; x['temp'] = 98.25
#
#x
#array([((10, 0), 98.25)],
#      dtype=[('time', [('min', '<i8'), ('sec', '<i8')]), ('temp', '<f8')])
