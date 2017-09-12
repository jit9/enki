import pickle
from numpy import *
from searchNED import *
from getFlux import *

class pixel:
    """ object stores pixel numbers and value at coords"""

    def __init__ (self, i, j, val, s_n, weight=1.):
        self.i = i
        self.j = j
        self.s_n = s_n
        self.val = val
        self.weight = weight

    def __str__ (self):
        
        return "%d %d %f" %(self.i, self.j, self.s_n)
        

class cluster:
    """ creates a cluster object"""

    def __init__(self):

        self.centx = 0.0      # x coord of centroid
        self.centy = 0.0      # y coord of centroid
        self.num = 0.0        # num of pixels
        self.s_nMax = 0.0     # greatest S/N in cluster
        self.valMax = 0.0     # greatest value in cluster
        self.pixels = []      # list of pixels in cluster
        self.ra = 0.0         # centroid ra
        self.dec = 0.0        # centroid dec
        self.weight = 0.

    def add (self, pixel):

        # add pixel to list
        self.pixels.append(pixel)

        # find new centroid
        numerX = 0.0
        numerY = 0.0
        denom = 0.0
        for pix in self.pixels:
            numerX += (pix.j * pix.s_n)
            denom += pix.s_n
            numerY += (pix.i * pix.s_n)

        self.centx = numerX / denom
        self.centy = numerY / denom

        self.weight = (self.weight*self.num + pixel.weight)/(self.num+1)
        # increment num of pix
        self.num += 1
        
        if abs(pixel.val)>abs(self.valMax):
            self.valMax = pixel.val
        if abs(pixel.s_n)>abs(self.s_nMax):
            self.s_nMax = pixel.s_n

    def addCluster (self, cluster):
        """ combines two cluster objects"""
        
        # combine pixels list
        self.pixels += cluster.pixels

        self.weight = (self.weight*self.num + cluster.weight*cluster.num)/(self.num+cluster.num)
        # calc new centroid
        numerX = 0.0
        numerY = 0.0
        denom = 0.0
        for pix in self.pixels:
            numerX += (pix.j * pix.s_n)
            denom += pix.s_n
            numerY += (pix.i * pix.s_n)
            if abs(pix.val)>abs(self.valMax):
                self.valMax = pix.val
            if abs(pix.s_n)>abs(self.s_nMax):
                self.s_nMax = pix.s_n

        self.centx = numerX / denom
        self.centy = numerY / denom

        # combine num of pix
        self.num = self.num + cluster.num
   
    
    def inContact (self, pixel):
        """ returns true if cluster is next to input pixel"""
        
        for k in range(0, len(self.pixels)): 
            diff_i = self.pixels[k].i - pixel.i
            diff_j = self.pixels[k].j - pixel.j

            if (diff_i**2 + diff_j**2) <= 3:
                return True
            
        return False
        

    def centroidPix (self):
        """ returns cluster centroid in pixels"""
        return self.centx, self.centy

    def __str__ (self):
                
        return "%d %f" %(self.num, self.s_nMax)

class clusterList( list ):

    def save( self, filename ):
        f = file(filename, 'w')
        pickle.dump( self, f )
        f.close()

    def read( self, filename ):
        f = file(filename)
        cl = pickle.load( f )
        self += cl
        f.close()

    def inPaintOverClusters( self, mp, supermap = None,  backgroundAnnulus0 = 5, backgroundAnnulus1 = 10 ):

        smhw = 0.1
        for cl in self:
            cosdec = numpy.cos(cl.dec*numpy.pi/180.)
            if supermap != None:
                smp = supermap.selectSubMap( cl.ra-smhw/cosdec, cl.ra+smhw/cosdec, cl.dec-smhw, cl.dec+smhw )
            else:
                smp = mp.selectSubMap( cl.ra-smhw/cosdec, cl.ra+smhw/cosdec, cl.dec-smhw, cl.dec+smhw )
            x = range(0, smp.Nx)
            y = range(0, smp.Ny)
            xx, yy = numpy.meshgrid(x,y)
            centX, centY = smp.skyToPix(cl.ra, cl.dec)
            Dist = numpy.sqrt((xx-centX)**2 + (yy-centY)**2)
            a, b = numpy.where((Dist > backgroundAnnulus0)*(Dist < backgroundAnnulus1))
            avg = numpy.mean(smp.data[a,b])
            for pix in cl.pixels:
                mp.data[pix.i,pix.j] = avg




class clusterFinder:
    """finds clusters in a liteMap object and creates catalog
    inputs: filteredMap, type = 'sources' or 'clusters', threshold = S/N max
    returns: list of cluster objects"""

    def __init__(self, filteredMap, type, threshold, minPixPerCluster = 4, rms = None, weight = None):
        
        if type == "sources":
            self.sign = 1
        else:
            self.sign = -1

        self.thresh = threshold
        self.wcs = filteredMap.wcs
        self.clusters = clusterList()
        self.data = filteredMap.data
        if weight != None:
            print weight.shape, self.data.shape
            self.weightedData = self.data*(weight/weight.max())
        else:
            self.weightedData = self.data
        self.currCluster = 0
        self.map = filteredMap 
        self.minPixPerCluster = minPixPerCluster
        self.rms = rms

    def getSigmaPixels(self):
        """ return pixel indices of the given threshold"""
       
        data = self.weightedData
        
        # mask clusters/sources to get new rms value
    
        if self.rms == None:
            rms = std(data)
            i1, j1 = where(data < - self.thresh * rms)
            i2, j2 = where(data > self.thresh * rms)
     
            pixelDist = 8
            
            noise = copy(data)
            noise[i1, j1] = data[i1-pixelDist, j1-pixelDist]
            noise[i2, j2] = data[i2-pixelDist, j2-pixelDist]
 
            rms = std(noise)
        else:
            rms = self.rms

        print "rms = ", rms
            # return either bright or dark pixels
        if self.sign == -1:
            i, j = where(data < - self.thresh * rms)
        else:
            i, j = where(data >  self.thresh * rms)

        # put data map in units of S/N
        self.snr = self.weightedData / (self.sign*rms)

        return i, j 
        

    def findClusters(self):
        """ find cluster objects, prints output to catalog file and creates regions file"""

        # get pixels of significance
        i, j = self.getSigmaPixels()

        # create corresponding pixel objects
        pix = []    
        for k in range(0, len(i)):
            pix.append(pixel(i[k], j[k], self.data[i[k], j[k]], self.snr[i[k], j[k]]))
        
        # create initial clusters list  
        for k in range(0, len(i)):

            # create first cluster object
            if k == 0:
                clust = cluster()
                self.clusters.append(clust)
                clust.add(pix[k])

            # check if pix is in contact with
            else:
                 
                if (self.clusters[self.currCluster].inContact(pix[k])):
                    self.clusters[self.currCluster].add(pix[k])

                else:
                    
                    self.currCluster = self.currCluster + 1
                    clust = cluster()                            
                    clust.add(pix[k])
                    self.clusters.append(clust)
                    
        # continue to loop thru list until no touching clusters left
        oldClusterList = self.clusters
        if len(oldClusterList ) == 0:
            redundantClustersExist = False
            newClusterList = oldClusterList
        else:
            redundantClustersExist = True
         

        while(redundantClustersExist):
            redundantClustersExist = False
            newClusterList = clusterList( [oldClusterList[0]] )

            for k in xrange(1, len(oldClusterList)):
              
                redundantCluster = False

                for cl in newClusterList:
                    for m in xrange(0, oldClusterList[k].num):
                        
                        if cl.inContact(oldClusterList[k].pixels[m]):
                            cl.addCluster(oldClusterList[k])
                            redundantCluster = True
                            redundantClustersExist = True
                            break
                    
                if not redundantCluster:
                    newClusterList.append(oldClusterList[k])
            oldClusterList = newClusterList
           
        # get nonredundant clust list
        self.clusters = newClusterList

        # eliminate small clusters and find S/N max
        clustersNew = clusterList()
        for cl in self.clusters:
            if cl.num >= self.minPixPerCluster:
                clustersNew.append(cl)         

        self.clusters = clustersNew

        # return list of clusters
        return self.clusters


    def setRADecOfClusters( self ):
        for cl in self.clusters: 
            cent = cl.centroidPix() 
            cl.ra, cl.dec =  self.map.pixToSky(cent[0], cent[1])


    def removeClustersAroundGiants( self, sig = 50, rad = 25 ):
        """
        Remove clusters in rad arcminutes around clusters of sig significance
        """
        superclusters = []
        if self.clusters[0].ra == 0.0 and self.clusters[0].dec == 0.0:
            self.setRADecOfClusters()
        for cl in self.clusters:
            if cl.s_nMax > sig:
                superclusters.append(cl)

        for scl in superclusters:
            cosdec = numpy.cos(scl.dec*numpy.pi/180.)
            toRemove=[]
            for cl in self.clusters:
                dra = (cl.ra-scl.ra)*cosdec
                ddec = (cl.dec - scl.dec)
                dist = numpy.sqrt(dra**2 + ddec**2)
                if dist*60 < rad and cl != scl:
                    toRemove.append(cl)
            for cl in toRemove:
                self.clusters.remove(cl)
