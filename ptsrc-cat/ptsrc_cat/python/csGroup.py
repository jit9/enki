import pickle
import numpy
from searchNED import *
# from getFlux import *


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
        

class group:
    """ creates a group object"""

    def __init__(self):

        self.centx = 0.0      # x coord of centroid
        self.centy = 0.0      # y coord of centroid
        self.num = 0.0        # num of pixels
        self.s_nMax = 0.0     # greatest S/N in group
        self.valMax = 0.0     # greatest value in group
        self.pixels = []      # list of pixels in group
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

    def addGroup (self, group):
        """ combines two group objects"""
        
        # combine pixels list
        self.pixels += group.pixels

        self.weight = (self.weight*self.num + group.weight*group.num)/(self.num+group.num)
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
        self.num = self.num + group.num
   
    
    def inContact (self, pixel):
        """ returns true if group is next to input pixel"""
        
        for k in range(0, len(self.pixels)): 
            diff_i = self.pixels[k].i - pixel.i
            diff_j = self.pixels[k].j - pixel.j

            if (diff_i**2 + diff_j**2) <= 3:
                return True
            
        return False
        

    def centroidPix (self):
        """ returns group centroid in pixels"""
        return self.centx, self.centy

    def __str__ (self):
                
        return "%d %f" %(self.num, self.s_nMax)

class groupList( list ):

    def save( self, filename ):
        f = file(filename, 'w')
        pickle.dump( self, f )
        f.close()

    def read( self, filename ):
        f = file(filename)
        gr = pickle.load( f )
        self += gr
        f.close()

    def inPaintOverGroups( self, mp, supermap = None,  backgroundAnnulus0 = 5, backgroundAnnulus1 = 10 ):

        smhw = 0.1
        for gr in self:
            cosdec = numpy.cos(gr.dec*numpy.pi/180.)
            if supermap != None:
                smp = supermap.selectSubMap( gr.ra-smhw/cosdec, gr.ra+smhw/cosdec, gr.dec-smhw, gr.dec+smhw, safe = True )
            else:
                smp = mp.selectSubMap( gr.ra-smhw/cosdec, gr.ra+smhw/cosdec, gr.dec-smhw, gr.dec+smhw, safe = True )
            x = range(0, smp.Nx)
            y = range(0, smp.Ny)
            xx, yy = numpy.meshgrid(x,y)
            centX, centY = smp.skyToPix(gr.ra, gr.dec)
            Dist = numpy.sqrt((xx-centX)**2 + (yy-centY)**2)
            a, b = numpy.where((Dist > backgroundAnnulus0)*(Dist < backgroundAnnulus1))
#             a = a[numpy.where((a < smp.Ny)*(a >= 0))]
#             b = b[numpy.where((b < smp.Ny)*(b >= 0))]
            avg = numpy.mean(smp.data[a,b])
            for pix in gr.pixels:
                mp.data[pix.i,pix.j] = avg




class groupFinder:
    """finds groups in a liteMap object and creates catalog
    inputs: filteredMap, type = 'sources' or 'groups', threshold = S/N max
    returns: list of group objects"""

    def __init__(self, filteredMap, type, threshold, minPixPerGroup = 4, rms = None, weight = None ):
        
        if type == "sources":
            self.sign = 1
        else:
            self.sign = -1

        self.thresh = threshold
        self.wcs = filteredMap.wcs
        self.groups = groupList()
        self.data = filteredMap.data
        if weight != None:
            self.weightedData = self.data*(weight/weight.max())
        else:
            self.weightedData = self.data
        self.currGroup = 0
        self.map = filteredMap 
        self.minPixPerGroup = minPixPerGroup
        self.rms = rms

    def getSigmaPixels(self):
        """ return pixel indices of the given threshold"""
       
        data = self.weightedData
        
        # mask groups/sources to get new rms value
    
        if self.rms == None:
            rms = numpy.std(data[numpy.where(data != 0.)])
            i1, j1 = numpy.where(data < - self.thresh * rms)
            i2, j2 = numpy.where(data > self.thresh * rms)
            print "raw rms: ", rms 
            pixelDist = 8
            
            noise = numpy.copy(data)
            noise[i1, j1] = data[i1-pixelDist, j1-pixelDist]
            noise[i2, j2] = data[i2-pixelDist, j2-pixelDist]
 
            rms = numpy.std(noise[numpy.where(data != 0.)]) # right now we use zeros as an effective source mask
            self.rms = rms
        else:
            rms = self.rms

        print "rms = ", rms
            # return either bright or dark pixels
        if self.sign == -1:
            i, j = numpy.where(data < - self.thresh * rms)
        else:
            i, j = numpy.where(data >  self.thresh * rms)

        # put data map in units of S/N
        self.snr = self.weightedData / (self.sign*rms)

        return i, j 
        

    def findGroups(self):
        """ find group objects, prints output to catalog file and creates regions file"""

        # get pixels of significance
        i, j = self.getSigmaPixels()

        # create corresponding pixel objects
        pix = []    
        for k in range(0, len(i)):
            pix.append(pixel(i[k], j[k], self.data[i[k], j[k]], self.snr[i[k], j[k]]))
        
        # create initial groups list  
        for k in range(0, len(i)):

            # create first group object
            if k == 0:
                grust = group()
                self.groups.append(grust)
                grust.add(pix[k])

            # check if pix is in contact with
            else:
                 
                if (self.groups[self.currGroup].inContact(pix[k])):
                    self.groups[self.currGroup].add(pix[k])

                else:
                    
                    self.currGroup = self.currGroup + 1
                    grust = group()                            
                    grust.add(pix[k])
                    self.groups.append(grust)
                    
        # continue to loop thru list until no touching groups left
        oldGroupList = self.groups
        if len(oldGroupList ) == 0:
            redundantGroupsExist = False
            newGroupList = oldGroupList
        else:
            redundantGroupsExist = True
         

        while(redundantGroupsExist):
            redundantGroupsExist = False
            newGroupList = groupList( [oldGroupList[0]] )

            for k in xrange(1, len(oldGroupList)):
              
                redundantGroup = False

                for gr in newGroupList:
                    for m in xrange(0, int(oldGroupList[k].num)):
                        
                        if gr.inContact(oldGroupList[k].pixels[m]):
                            gr.addGroup(oldGroupList[k])
                            redundantGroup = True
                            redundantGroupsExist = True
                            break
                    
                if not redundantGroup:
                    newGroupList.append(oldGroupList[k])
            oldGroupList = newGroupList
           
        # get nonredundant grust list
        self.groups = newGroupList

        # eliminate small groups and find S/N max
        groupsNew = groupList()
        print "Found %d total groups" % len(self.groups)
        for gr in self.groups:
            if gr.num >= self.minPixPerGroup:
                groupsNew.append(gr)         
        print "Kept %d groups after number-of-pixels cuts" % len(groupsNew)
        self.groups = groupsNew

        # return list of groups
        return self.groups


    def setRADecOfGroups( self ):
        for gr in self.groups: 
            cent = gr.centroidPix() 
            gr.ra, gr.dec =  self.map.pixToSky(cent[0], cent[1])


    def removeGroupsAroundGiants( self, sig = 50, rad = 25 ):
        """
        Remove groups in rad arcminutes around groups of sig significance
        """
        supergroups = []
        if self.groups[0].ra == 0.0 and self.groups[0].dec == 0.0:
            self.setRADecOfGroups()
        for gr in self.groups:
            if gr.s_nMax > sig:
                supergroups.append(gr)

        for sgr in supergroups:
            cosdec = numpy.cos(sgr.dec*numpy.pi/180.)
            toRemove=[]
            for gr in self.groups:
                dra = (gr.ra-sgr.ra)*cosdec
                ddec = (gr.dec - sgr.dec)
                dist = numpy.sqrt(dra**2 + ddec**2)
                if dist*60 < rad and gr != sgr:
                    toRemove.append(gr)
            for gr in toRemove:
                self.groups.remove(gr)
