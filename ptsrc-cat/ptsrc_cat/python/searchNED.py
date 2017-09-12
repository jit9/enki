import urllib2
import re

def searchNED (ra, dec, radius = 0.5):

    searchString = 'http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?search_type=Near+Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon=%.4fd&lat=%.4fd&radius=%.3f&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=xml_names&zv_breaker=30000.0&list_limit=5&img_stamp=YES' %(ra, dec, radius)

    site = urllib2.urlopen(searchString)

    lines = site.readlines()
    for line in lines:
        print line
    site.close()
    radio = '--'
    for line in lines:
        if 'PMN ' in line and '<TD>' in line:
            radio = line[4:-6]
    if radio == '--':
        for line in lines:
            if 'SUMSS ' in line and '<TD>' in line:
                radio = line[4:-6]

    ir = '--'
    for line in lines:
        if 'IRAS  ' in line and '<TD>' in line and not 'ID' in line:
            ir = line[4:-6]
            print ir
    if ir == '--':
        for line in lines:
            if 'IRAS ' in line and '<TD>' in line:
                ir = line[4:-6]
    if ir == '--':
        for line in lines:
            if '2MASX ' in line and '<TD>' in line:
                ir = line[4:-6]

    return radio, ir


mainDict = {0:['Number',int], 1:['Name',str], 2:['RA',float], 3:['Dec',float], 4:['ID',str], 5:['Velocity',float], 6:['Redshift',float]}

def searchNEDMain(ra, dec, radius = 0.5):

    searchString = 'http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?search_type=Near+Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon=%.4fd&lat=%.4fd&radius=%.3f&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=xml_main&zv_breaker=30000.0&list_limit=5&img_stamp=YES' %(ra, dec, radius)

    site = urllib2.urlopen(searchString)

    lines = site.readlines()
#     for line in lines:
#         print line
    site.close()
    IDs = []
    for line in lines:
        if '<TR>' in line:
            IDs.append({})
            i = 0
            continue
        if '<TD>' in line and i < len(mainDict.keys()):
            IDs[-1][mainDict[i][0]] = mainDict[i][1](line[4:-6])
            i+=1
    return IDs


def searchNEDForFluxes( name ):
#     searchString = 'http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?search_type=Near+Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon=%.4fd&lat=%.4fd&radius=%.3f&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=xml_basic&zv_breaker=30000.0&list_limit=5&img_stamp=YES' %(ra, dec, radius)
    name=name.replace(" ", "%20")
    print name
    searchString="http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?extend=no&of=xml_all&objname=%s" % name

    site = urllib2.urlopen(searchString)
    lines = site.readlines()
    site.close()
    col = 0
    for line in lines:
        if "<TD>" in line:
            col+=1
            print col, line
        else:
            print line
#         if col == 8:
#             flux = 10**float(line.strip(" </TD>"))
#         if col == 9:
#             freq = float(line.strip(" </TD>"))
#     return freq, flux


def searchNEDForClusters (ra, dec, radius = 2.0):
    
    searchString = 'http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?search_type=Near+Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon=%.4fd&lat=%.4fd&radius=%.3f&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=xml_names&zv_breaker=30000.0&list_limit=5&img_stamp=YES&in_objtypes1=GClusters' %(ra, dec, radius)

    site = urllib2.urlopen(searchString)

    lines = site.readlines()
    site.close()
    tabledata=False
    for line in lines:
        if '<TABLEDATA>' in line:
            tabledata = True
        if not tabledata:
            continue
        if 'GClstr' in line:
            return lastline[4:-6]
        lastline = line
    return '--'
