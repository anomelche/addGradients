import os, sys, numpy, shutil
from ADFR.utils.maps import MapsFile
from Volume.IO.AutoGridWriter import WriteAutoGrid
from Volume.Grid3D import Grid3DF
from ADFRcc.adfrcc import AddGradients 

class addGradients:

    def __init__(self):
        self.mf = None
        self.spacing = None
        self._seen = {}
        self._mapdata = None
        self._writer = WriteAutoGrid()
        self._gradLow = 0.0
        
    def setTarget(self, target):
        if isinstance(target, str):
            self.mf = MapsFile(target)
        else:
            assert isinstance(target, MapsFile)
            self.mf = target
        self.spacing = self.mf.getBoxSpacing()
        self.mf.unzipMaps()
        self.mapTypes = self.mf.getMapTypes()
        map_ = self.mf.getMap('e')
        self._emapdata = map_.getGridDataPy()
        self._emin = min(self._emapdata.flatten())
        self._emax = max(self._emapdata.flatten())
        self._gradLow = int(-self._emin)+1
        print 'ELEC', self._emin, self._emax, self._gradLow
        
    def writeADmap(self, filename, ar, origin, spacing):
        grid = Grid3DF(ar, origin, (spacing, spacing, spacing),
                       {'GRID_PARAMETER_FILE':'None',
                        'GRID_DATA_FILE':'None',
                        'MACROMOLECULE':'None'})
        self._writer.write(grid, filename)

    
    def processMaps(self, outTargetFile, errorCut=0.01):
        print 'processing', 
        #typelist = self.mapTypes[:]
        #typelist.remove('HD')
        #typelist = ['HD'] + typelist
        #for atype in typelist:
        for atype in self.mapTypes:
            if atype in ['e', 'd', 'sd']: continue
            print atype,
            sys.stdout.flush()
            self._seen = {}
            map_ = self.mf.getMap(atype)
            origin = map_.getOriginPy().astype('f')
            self._mapdata = map_.getGridDataPy()

            ## loop to zero out small positive values due to desolvation
            nx, ny, nz = self._mapdata.shape
            for a in range(nx):
                for b in range(ny):
                    for c in range(nz):
                        if self._mapdata[a,b,c]>0.0 and self._mapdata[a,b,c]<errorCut:
                            self._mapdata[a,b,c] = 0.0

            ## find all grid points with negative values
            i,j,k = numpy.where(self._mapdata<=self._gradLow)
            indices = numpy.zeros((len(i),3),'i')
            indices[:,0] =  i
            indices[:,1] =  j
            indices[:,2] =  k
            neighborPts = 1
            cutOffValue = 16
            errorCut=0.01
            # C++ class
            addgrad = AddGradients(self._gradLow)
            addgrad.processMap(self._mapdata, indices, (self.spacing,self.spacing,self.spacing), neighborPts, cutOffValue, errorCut)
            # the values in self._mapdata are changed by the c++ processMap method.
            # which also computes the clusters
            
            # this is how to get the clusters
            clusters = addgrad.clusters
            clen = [len(x) for x in clusters]
            print clen
            ## for i, cl in enumerate(clusters):
            ##     if len(cl) <=100:
            ##         print "cluster ", i , " " , cl  
            #import pdb; pdb.set_trace()
            self.writeADmap(os.path.join(self.mf.folder, 'rigidReceptor.%s.map'%atype),
                            self._mapdata.astype('f'), origin, self.spacing)
        name = outTargetFile
        print 'creating target file ...', './%sGrad.trg'%name
        shutil.move(self.mf.folder, self.mf.folder+'Grad')
        filename = make_archive(self.mf.folder+'Grad', './%sGrad.zip'%name)
        
        #filename = shutil.make_archive(
        #    name+'Grad', 'zip', os.path.split(self.mf.folder)[0], './%sGrad'%name)
        shutil.move('./%sGrad.zip'%name, './%sGrad.trg'%name)
        print 'Done'

def make_archive(source, destination):
    # make_archive('/path/to/folder', 'path/to/folder.zip')
    base = os.path.basename(destination)
    name, _format = base.split('.')
    archive_from = os.path.dirname(source)
    archive_to = os.path.basename(source.strip(os.sep))
    filename = shutil.make_archive(name, _format, archive_from, archive_to)
    shutil.move('%s.%s'%(name, _format), destination)
    return filename

if __name__=='__main__':
    import sys
    filename = '7cpa_rec.trg'
    #filename = sys.argv[1]
    processor = addGradients()
    processor.setTarget(filename)
    #name = os.path.splitext(os.path.basename(filename))[0]
    name = os.path.splitext(filename)[0]
    processor.processMaps(name)
    #writeADmap('cgrad.Cg.map', cdata.astype('f'), (ox, oy, oz), spacing)
