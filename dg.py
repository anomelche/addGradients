import os, sys, numpy, shutil
from ADFR.utils.maps import MapsFile
from Volume.IO.AutoGridWriter import WriteAutoGrid
from Volume.Grid3D import Grid3DF

class DensityClustering:

    def __init__(self, spacing, neighborPts=14):
        self._indices = None
        self.neighborpts = neighborPts
        self.spacing = spacing
        self._clusters = None
        self._clen = []
        
    def findClustersD(self, ptsIndices, cVolcut=20):
        """returns list of contiguous clusters of filter points with every point in a
        cluster having the specified number of neighbors. """
	#print "wathch me", self.spacing
        # get number of neighbors needed for density clustering
        nneighbor = self.neighborpts
        sx,sy,sz = self.spacing
        si = sj = sk = 1.0
        clusters = [] # list of clusters identified
        cluster = [] # current cluster to build

        # keycache
        keycache = {}
        indcache = {}
        n = 0
        for i,j,k in ptsIndices:
            key = '%d,%d,%d'%(i,j,k)
            keycache[(i,j,k)] = key
            indcache[(i,j,k)] = n
            n += 1
            #if i%si == 0 and j%sj==0 and k%sk==0:
            #    newptsIndices.append((i,j,k))
            
        
        indices = {} # indices of points not clustered yet
        #n = 0
        for ii in range(len(ptsIndices)):
            i,j,k = ptsIndices[ii]
            if i%si == 0 and j%sj==0 and k%sk==0:
                indices[ii] = ii
                #n += 1
        #import pdb
        #pdb.set_trace()
        ci = 0 # index of point in cluster
        index = indices.values()[0]
        i, j, k = ptsIndices[index] # first point
        cluster.append( index )
        indices.pop(index)
        #indices.pop(indices.keys()[0])
        incluster = {index:True}
        while True:
            neighbors = [] # list of neighbors of current point

            # build list of grid points neighboring index and in ptsIndices
            if keycache.has_key( (i-si,j,  k  ) ): neighbors.append( indcache[(i-si,j,  k  )] )
            if keycache.has_key( (i+si,j,  k  ) ): neighbors.append( indcache[(i+si,j,  k  )] )
            if keycache.has_key( (i,  j-sj,k  ) ): neighbors.append( indcache[(i,  j-sj,k  )] )
            if keycache.has_key( (i,  j+sj,k  ) ): neighbors.append( indcache[(i,  j+sj,k  )] )
            if keycache.has_key( (i,  j,  k-sk) ): neighbors.append( indcache[(i,  j  ,k-sk)] )
            if keycache.has_key( (i,  j,  k+sk) ): neighbors.append( indcache[(i,  j  ,k+sk)] )

            if keycache.has_key( (i-si,j+sj,  k  ) ): neighbors.append( indcache[(i-si,j+sj,  k  )] )
            if keycache.has_key( (i+si,j+sj,  k  ) ): neighbors.append( indcache[(i+si,j+sj,  k  )] )
            if keycache.has_key( (i-si,j-sj,  k  ) ): neighbors.append( indcache[(i-si,j-sj,  k  )] )
            if keycache.has_key( (i+si,j-sj,  k  ) ): neighbors.append( indcache[(i+si,j-sj,  k  )] )
            if keycache.has_key( (i-si,  j ,k+sk ) ): neighbors.append( indcache[(i-si,  j ,k+sk )] )
            if keycache.has_key( (i+si,  j ,k+sk ) ): neighbors.append( indcache[(i+si,  j ,k+sk )] )
            if keycache.has_key( (i-si,  j ,k-sk ) ): neighbors.append( indcache[(i-si,  j ,k-sk )] )
            if keycache.has_key( (i+si,  j ,k-sk ) ): neighbors.append( indcache[(i+si,  j ,k-sk )] )
            if keycache.has_key( (i   ,j-sj,k+sk ) ): neighbors.append( indcache[(i   ,j-sj,k+sk )] )
            if keycache.has_key( (i   ,j+sj,k+sk ) ): neighbors.append( indcache[(i   ,j+sj,k+sk )] )
            if keycache.has_key( (i   ,j-sj,k-sk ) ): neighbors.append( indcache[(i   ,j-sj,k-sk )] )
            if keycache.has_key( (i   ,j+sj,k-sk ) ): neighbors.append( indcache[(i   ,j+sj,k-sk )] )

            if keycache.has_key( (i+si,j+sj,k+sk ) ): neighbors.append( indcache[(i+si,j+sj,k+sk )] )
            if keycache.has_key( (i-si,j+sj,k+sk ) ): neighbors.append( indcache[(i-si,j+sj,k+sk )] )
            if keycache.has_key( (i+si,j-sj,k+sk ) ): neighbors.append( indcache[(i+si,j-sj,k+sk )] )
            if keycache.has_key( (i+si,j+sj,k-sk ) ): neighbors.append( indcache[(i+si,j+sj,k-sk )] )
            if keycache.has_key( (i-si,j-sj,k+sk ) ): neighbors.append( indcache[(i-si,j-sj,k+sk )] )
            if keycache.has_key( (i-si,j+sj,k-sk ) ): neighbors.append( indcache[(i-si,j+sj,k-sk )] )
            if keycache.has_key( (i+si,j-sj,k-sk ) ): neighbors.append( indcache[(i+si,j-sj,k-sk )] )
            if keycache.has_key( (i-si,j-sj,k-sk ) ): neighbors.append( indcache[(i-si,j-sj,k-sk )] )
            
                
            if len(neighbors)>=nneighbor: # if point index (i,j,k) has enough neighbors
                neigh = [x for x in neighbors if not incluster.has_key(x)]
                cluster.extend(neigh)
                for n in neigh:
                    incluster[n] = True
                    indices.pop(n)             

            ci += 1
            if ci<len(cluster): # more points in this cluster can be checked for neighbors
                i,j,k = ptsIndices[cluster[ci]]
            else: # all points in cluster have been checked for neighbors
                #print 'CLUSTER', len(cluster)
                clusters.append(cluster)
                if len(indices)==0:
                     break
                cluster = []
                ci = 0
                index = indices.values()[0]
                i,j,k = ptsIndices[index]
                cluster.append( index )
                indices.pop(index)
                #indices.pop(indices.keys()[0])
                incluster[index] = True
                
        clen = [len(x) for x in clusters]  
        cl = []
        for i in numpy.argsort(clen)[::-1]:
            cl.append(clusters[i])#numpy.array(clusters[i], 'i'))
        clv = [x for x in cl if len(x) >= cVolcut]
        cl = clv
        #print 'STEP1', (clen)
        #import pdb;pdb.set_trace()
        #self._clusters = cl # clusters sorted by size
        #self._clen = [len(x) for x in cl]   # sorted clusters lengths
        self._clusters = cl # clusters sorted by size
        self._clen = [len(x) for x in cl]   # sorted clusters lengths

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

    def onepass(self, cdata, indices, generation):
        edge = []
        seen = self._seen
        for a,b,c in indices:
            for a1,b1,c1 in [[a-1,b,c],[a+1,b,c],[a,b-1,c],[a,b+1,c],[a,b,c-1],[a,b,c+1]]:
                if seen.get((a1,b1,c1), True):
                    continue
                edge.append((a1, b1, c1))
                seen[(a1, b1, c1)] = True
                cdata[(a1, b1, c1)] = generation
        return edge
    
    def processMaps1(self, outTargetFile, errorCut=0.01):
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
            print "pts array length: ",  len(indices)
            cc = 0
            for ii in indices:
                if cc < 10:
                    print "point ", cc,  "is ", ii[0], " ", ii[1], " ", ii[2]
                cc = cc+1
        return

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
            ## cluster grid poitns with negative values
            sc = DensityClustering((self.spacing,self.spacing,self.spacing), neighborPts=1)
            sc.findClustersD(indices, cVolcut=1)
            assert len(indices)==numpy.sum(sc._clen)
        
            for a in range(nx):
                for b in range(ny):
                    for c in range(nz):
                        self._seen[(a,b,c)] = False

            ## keep largest cluster of negative values as the "outside"
            for x in indices[sc._clusters[0]]:
                self._seen[tuple(x)] = True

            ## if atype=='HD': # 0 out electrostatic values inside the receptor
            ##     edata = self._emapdata
            ##     for x in indices[sc._clusters[0]]:
            ##         edata[tuple(x)] = 0.0

            ##     self._gradLow = int(-min(self._emapdata.flatten()))+1
            ##     self.writeADmap(os.path.join(self.mf.folder, 'rigidReceptor.e.map'),
            ##                     self._emapdata.astype('f'), origin, self.spacing)
            ##     print 'ELEC2', min(self._emapdata.flatten()), max(self._emapdata.flatten()), self._gradLow
            #import pdb; pdb.set_trace()
            generation = self._gradLow
            #generation = 1
            ind = indices[sc._clusters[0]]
            while(len(ind)):
                #print generation, len(ind)
                ind = self.onepass(self._mapdata, ind, generation)
                generation += 1

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
    #filename = '7cpa_rec.trg'
    filename = sys.argv[1]
    processor = addGradients()
    processor.setTarget(filename)
    #name = os.path.splitext(os.path.basename(filename))[0]
    name = os.path.splitext(filename)[0]
    processor.processMaps(name)
    #writeADmap('cgrad.Cg.map', cdata.astype('f'), (ox, oy, oz), spacing)
