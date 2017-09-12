#
# A module for handling catalogs of sources
# 
# Notes to self: should we make "row" its own class?
#


import numpy
import pickle, copy
from flipper import *
import pyfits

defaultCols = { \
        'ra'   : { 'desc':'Right Ascension of the object'   , 'order': 0, 'type': float, 'fmt' : '%10.5f', 'default': -999.}, \
        'dec'  : { 'desc':'Declination of the object'       , 'order': 1, 'type': float, 'fmt' : '%10.5f', 'default': -999.}, \
        's/n'  : { 'desc':'signal to noise of the detection', 'order': 2, 'type': float, 'fmt' : '%10.5f', 'default': -999.}
        }

def deltaT2JyPerSr(freqGHz,T0 = 2.726):
    """
    @return the multiplicative conversion factor to go from deltaT_cmb to JyPerSr
    """
    kB = 1.380658e-16
    h = 6.6260755e-27
    c = 29979245800.
    nu = freqGHz*1.e9
    x = h*nu/(kB*T0)
    cNu = 2*(kB*T0)**3/(h**2*c**2)*x**4/(4*(numpy.sinh(x/2.))**2)
    cNu *= 1e23/T0
    return cNu

def deltaT2y(freqGHz, T0 = 2.726):
    """
    @return the multiplicative conversion factor to go from deltaT_cmb to compton y 
    """
    kB = 1.380658e-16
    h = 6.6260755e-27
    c = 29979245800.
    nu = freqGHz*1.e9
    x = h*nu/(kB*T0)
    f_nu = x*(numpy.exp(x)+1)/(numpy.exp(x)-1) - 4
    return 1./f_nu/T0

def compare_by (fieldname): 
    def compare_two_dicts (a, b): 
        return cmp(a[fieldname], b[fieldname]) 
    return compare_two_dicts 

def convertRADecDegreesToSexagesimal(ra, dec):
    """
    Accepts float values for ra and dec
    returns ra, dec in the form (hh,mm,ss), abs(deg,mm,ss), sign
    where the absolute value of the declination is returned and the 
    sign is the sign of the declination ('+' or '-')
    """
    if ra < 0:
        ra += 360.
    ra_hr = int(ra/15.)
    minsec = numpy.mod(ra,15.)/15.*60
    ra_min = int(minsec)
    ra_sec = numpy.mod(minsec,1.)*60

    if dec >= 0.:
        sign = '+'
    else:
        sign = '-'
    dec_deg = int(abs(dec))
    minsec = numpy.mod(abs(dec),1.)*60
    dec_min = int(minsec)
    dec_sec = numpy.mod(minsec,1.)*60
    return (ra_hr, ra_min, ra_sec), (dec_deg, dec_min, dec_sec), sign

def makeIAUName( prefix, ra, dec):
    ra_s, dec_s, sign = convertRADecDegreesToSexagesimal(ra, dec)
    return prefix +  ' J%02d%02d%04.1f%s%02d%02d%02d' % (ra_s[0], ra_s[1], ra_s[2], sign ,dec_s[0], dec_s[1], numpy.round(dec_s[2]))

def convertSexagesimalToRaDecDegrees(ra, dec):
    """
    Assumes ra, dec arguments are of the form "hh:mm:ss", "deg:mm:ss"
    Returns float values for ra and dec
    """
    ra_hr, ra_min, ra_sec= ra.split(':')

    sign = ra_hr[0]
    if sign == '-':
        sign_ra = -1.
        ra_hr = ra_hr[1:]
    else:
        sign_ra = 1.

    ra = sign_ra * 15. * ( float(ra_hr) + float(ra_min)/60. + float(ra_sec)/3600. )

    dec_deg, dec_min, dec_sec= dec.split(':')

    sign = dec_deg[0]
    if sign == '-':
        sign_dec = -1.
        dec_deg  = dec_deg[1:]
    else:
        sign_dec = 1.

    dec =  sign_dec * (float(dec_deg) + float(dec_min)/60. + float(dec_sec)/3600.)

    return ra, dec

class catalog( list ):

    def __init__( self, cols = None):
        if cols == None:
            self.cols = defaultCols
        else:
            self.cols = cols
        self.ncol = len(self.cols.keys())
        self.sep = ' '

    def select( self, colName, func ):
        sub = catalog(cols = self.cols)
        for row in self:
            if colName != None:
                if func(row[colName]):
                    sub.append(copy.copy(row))
            else:
                if func(row):
                    sub.append(copy.copy(row))
        return sub

    def addRow( self, row ):
        newRow = {}
        inputKeys = row.keys()
        for k in self.cols.keys():
            if k in inputKeys:
                newRow[k]=row[k]
            else:
                newRow[k]=self.cols[k]['default']
        self.append(newRow)

    def addCol( self, name, desc, order, dtype, fmt, default, defaultFmt = None ):
        """
        Add a column to the catalog
        """
        cnames = self.cols.keys()
        if name in cnames:
            raise ValueError("Name %s already in column list" % name)
        for n in cnames:
            if self.cols[n]['order'] >= order:
                self.cols[n]['order'] += 1
        if defaultFmt:
            self.cols[name] = {'desc': desc   , 'order': order, 'type': dtype, 'fmt' : fmt, 'default': default, 'defaultFmt': defaultFmt}
        else:
            self.cols[name] = {'desc': desc   , 'order': order, 'type': dtype, 'fmt' : fmt, 'default': default}
        self.ncol = len(self.cols.keys())


    def sortBy( self, col ):
        self.sort(compare_by(col))
        

    def __str__( self ):
        string = ""
        colNames = self.cols.keys()
        for colName in self.cols.keys():
            colNames[self.cols[colName]['order']] = colName
        for row in self:
            for k in colNames:
                if self.cols[k]['default'] == row[k] and 'defaultFmt' in self.cols[k].keys():
                    string += self.cols[k]['defaultFmt'] % row[k]
                    string += self.sep
                    continue
                try:
                    string += self.cols[k]['fmt'] % row[k]
                except:
                    print "Failed to format col %s with this value %s." % (k, str(row[self.cols[k]['order']])), self.cols[k]['default'] == row[k], self.cols[k]['default']
                    raise
                string += self.sep
            string += '\n'
        return string

    def colNamesSorted( self ):
        colNames = self.cols.keys()
        for colName in self.cols.keys():
            colNames[self.cols[colName]['order']] = colName
        return colNames

    def colsString( self, comment = False ):
        colNames = self.colNamesSorted()
        cs = ""
        for i in range(self.ncol):
            if comment:
                cs += "# %d : %s - %s\n" % (i, colNames[i], self.cols[colNames[i]]['desc'])
            else:
                cs += "%d : %s - %s\n" % (i, colNames[i], self.cols[colNames[i]]['desc'])
        return cs

    def writeRegionFile( self, filename, radius = .1, raCol = 'ra', decCol = 'dec', color = 'black', select = None ):
        string = ""
        f = file( filename, 'w' )
        header = "global color=%s font='helvetica 10 normal' select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n" % color 
        f.write(header)
        f.write("wcs\n")
        for row in self:
            if select == None:
                if ('s/n' in row.keys() and row['s/n'] >= 5) or ('snr' in row.keys() and row['snr'] >= 5):
                    f.write("circle( %f, %f, .1 )#color=red\n" % (row[raCol], row[decCol]))
                else:
                    f.write("circle( %f, %f, .1 )\n" % (row[raCol], row[decCol]))
            else:
                selected = False
                for sel in select:
                    if sel[0](row):
                        f.write("circle( %f, %f, .1 )#color=%s\n" % (row[raCol], row[decCol], sel[1]))
                        selected = True
                if not selected:
                    f.write("circle( %f, %f, .1 )\n" % (row[raCol], row[decCol]))
        f.close()

    def read( self, filename, binary = False ):
        if binary:
            pass

        f = file(filename)
        lines = f.readlines()
        for line in lines:
            if line[0] == '#' or len(line) <= 1:
                continue
            fields = line.split()
            if fields[0][0] == '#' or len(fields) == 0:
                continue
            self.append({})
            for colName in self.cols.keys():
                col = self.cols[colName]
                if col['order'] < 0: 
                    continue
                if fields[col['order']] == col['default']:
                    self[-1][colName] = col['default']
                else:
                    try:
                        self[-1][colName] = col['type'](fields[col['order']])
                    except:
                        try:
                            self[-1][colName] = eval(col['type'])(fields[col['order']])
                        except:
                            print "Could not read column %d: %s with type %s and value %s" % (col['order'], colName, col['type'], fields[col['order']])
        f.close()

    def write( self, filename ):
        f = file(filename, 'w')
        pickle.dump( self, f )
        f.close()

    def writeASCII( self, filename, header=True, mode = 'w' ):
        f = file( filename, mode ) 
        if header:
            f.write(self.colsString(comment=True))
        f.write(self.__str__())
        fd = flipperDict.flipperDict()
        fd['cols'] = self.cols
        fd.writeToFile("%s.dict" % filename)

    def arrayFromCol(self, colName):
        c = []
        for row in self:
            c.append(row[colName])
        return numpy.array(c)

    #XXX under construction
    def plotCol1VsCol2(self, col1, col2, err1 = None, err2 = None ):
        cat = catalog.readFromPickleFile("merged.cat")
        rf = cat.arrayFromCol("Recovered_val")
        cf = cat.arrayFromCol("InputRadio_S_148")
        cf /=1000.
        re = cat.arrayFromCol("Recovered_err")
        # pylab.errorbar(cf[numpy.where(cf!=-.999)], rf[numpy.where(cf!=-.999)], yerr=re[numpy.where(cf!=-.999)], fmt='.')
        pylab.errorbar(cf[numpy.where(cf!=-.999)], rf[numpy.where(cf!=-.999)], yerr=None, fmt='.')
        pylab.gca().set_yscale('log')
        pylab.gca().set_xscale('log')
        pylab.plot([.0001,10],[.0001,10])
        utils.saveAndShow()

    def findNearestSource( self, row, rad = None ):
        """
        Given a row/source, find the closest row/source in this catalog
        @return distance (deg), closest row
        """
        closestRow = self[0]
        distToClosestRow = distanceBetweenSources(closestRow, row)
        rowsWithinRadius = []
        if rad != None and distToClosestRow < rad:
            rowsWithinRadius.append([distToClosestRow, closestRow])
        for myrow in self[1:]:
            d = distanceBetweenSources(row, myrow)
            if d < distToClosestRow:
                distToClosestRow = d
                closestRow = myrow
            if rad != None and d < rad:
                rowsWithinRadius.append([d,myrow])
        return distToClosestRow, closestRow, rowsWithinRadius

def readFromASCIIFile( filename, colFile= None, binary = False ):
    if binary:
        pass

    if colFile == None:
        colFile = "%s.dict" % filename
    colDict = flipperDict.flipperDict()
    colDict.readFromFile(colFile)
    cat = catalog(cols = colDict['cols'])
    cat.read(filename)
    return cat

def readFromFITS( filename, hdu = 1, rows = None, selection = None ):
    """
    hdu -- which header/data unit to read in
    rows -- if None read in all rows, if [start, stop] read in rows between start stop
    """
    hdulist = pyfits.open(filename)
    if selection != None:
        inds = selection(hdulist[hdu].data)
        data = hdulist[hdu].data[inds]
    else:
        data = hdulist[hdu].data
    cols = hdulist[hdu].columns
    _cat = catalog(cols={})
    for i in xrange(len(cols.names)):
        format = cols.formats[i]
        if format[-1] in 'EFD':
            tp = float
            fmt = "%20.6e"
            default = -999.
        if format[-1] in 'IL':
            tp = int 
            fmt = "%20d"
            default = -999
        if format[-1] == 'A':
            tp = str
            fmt = "%20s"
            default = '--'
        _cat.addCol(cols.names[i], cols.units[i], len(_cat.cols.keys()), tp, \
                fmt, default, defaultFmt = None )
    print "Rows to read: ", len (data)
    if rows == None:
        rowstoread = data
    else:
        rowstoread = data[rows[0]:rows[1]]
    for row in rowstoread:
        r = {}
        for c in _cat.cols.keys():
            r[c] = row[_cat.cols[c]['order']]
        _cat.append(r)
    hdulist.close()
    return _cat

def readFromPickleFile( filename ):
    f = file(filename)
    cat = pickle.load(f)
    f.close()
    return cat

fromClusterCols = { \
        'ra'   : { 'desc':'Right Ascension of the object'   , 'order': 0, 'type': float, 'fmt' : '%10.5f', 'default': -999.}, \
        'dec'  : { 'desc':'Declination of the object'       , 'order': 1, 'type': float, 'fmt' : '%10.5f', 'default': -999.}, \
        'val'  : { 'desc':'val of the detection',           'order': 2, 'type': float, 'fmt' : '%10.3e', 'default': -999.}, \
        'err'  : { 'desc':'1-sigma error of the detection', 'order': 3, 'type': float, 'fmt' : '%10.3e', 'default': -999.}, \
        's/n'  : { 'desc':'signal to noise of the detection', 'order': 4, 'type': float, 'fmt' : '%10.3f', 'default': -999.}, \
        'npix' : { 'desc':'Number of pixels in cluster'     , 'order': 5, 'type': int,   'fmt' : '%10d  ', 'default': -999 }  \
        }


def catalogFromGroupList( clusterList ):
    cat = catalog(cols = fromClusterCols)
    for cl in clusterList:
        err = cl.valMax*(1./cl.s_nMax)
        cat.addRow( {'ra':cl.ra, 'dec':cl.dec, 's/n':cl.s_nMax, 'err':err, 'val':cl.valMax, 'npix':len(cl.pixels)} )
    return cat

def mergeCatalogs( cat1, name1, cat2, name2, dist=1.5, cat1Only = False ):
    """
    dist = distance between objects at which point you declare a match in arcmin
    """
    newCols = { }
    colNames1 = cat1.cols.keys()
    colNames2 = cat2.cols.keys()
    for colName in colNames1:
        if name1 != None:
            newCols["%s_%s" % (name1, colName)] = copy.deepcopy(cat1.cols[colName] )
        else:
            newCols[colName] = cat1.cols[colName] 
    for colName in colNames2:
        newCols["%s_%s" % (name2, colName)] = copy.deepcopy(cat2.cols[colName])
        newCols["%s_%s" % (name2, colName)]['order'] += len(colNames1)
    newCols['matched'] = {'desc':'match between catalogs','order': len(colNames1)+len(colNames2), 'type': float, 'fmt' : ' %s ', 'default': "N"}
    cat = catalog(cols = newCols)
    double2 = []
    for row1 in cat1:
        row = {}
        for colName in colNames1:
            if name1 != None:
                row["%s_%s" % (name1, colName)] = row1[colName]
            else:
                row[colName] = row1[colName]
        matches = False
        for row2 in cat2:
            cosdec = numpy.cos((row1['dec'] + row2['dec'])/2 * numpy.pi/180)
            xdiff = (row1['ra']-row2['ra'])*cosdec
            ydiff = row1['dec'] - row2['dec']
            d = numpy.sqrt(xdiff**2 + ydiff**2)
            if d < dist/60.:
                for colName in colNames2:
                    row["%s_%s" % (name2, colName)] = row2[colName]
                double2.append(row2)
                matches = True
                row['matched'] = 'Y'
        if not matches:
            for colName in colNames2:
                row["%s_%s" % (name2, colName)] = cat2.cols[colName]['default']
            row['matched'] = 'N'
        cat.addRow(row)
    if not cat1Only:
        for row2 in cat2:
            if row2 in double2:
                continue
            row = {}
            for colName in colNames1:
                row["%s_%s" % (name1, colName)] = cat1.cols[colName]['default']
            for colName in colNames2:
                row["%s_%s" % (name2, colName)] = row2[colName]
            row['matched'] = 'N'
            cat.addRow(row)
    return cat

def getDecimalRADecFromRow( row ):
    if 'ra' in row.keys():
        return row['ra'], row['dec']
    elif 'ra_s' in row.keys():
        return convertSexagesimalToRaDecDegrees( row['ra_s'], row['dec_s'] )
    else:
        raise ValueError('No key ra or ra_s in row: %s' % str(row))

def distanceBetweenSources(row1, row2):
    """
    Compute distance in degrees between to catalog entries
    """
    ra1, dec1 = getDecimalRADecFromRow( row1 )
    ra2, dec2 = getDecimalRADecFromRow( row2 )
    cosdec = numpy.cos(numpy.pi/180. * (dec1+dec2)/2)
    dra = (ra1-ra2)*cosdec
    ddec = dec1-dec2
    return (dra**2+ddec**2)**0.5


def read(catname):
    try:
        cat = readFromPickleFile( catname )
    except:
        cat = readFromASCIIFile( catname )
    return cat
