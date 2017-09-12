from flipper import *

def dndsVec( cat, col, scale, bins, area ):
    colData = cat.arrayFromCol(col)*scale
    bin0 = bins[0] 
    dnds = []
    sdnds = []
    svec = []
    for bin1 in bins[1:]:
        inds = numpy.where( (colData>bin0) * (colData<bin1) )
        dataInBin = colData[inds]
        dnds.append(len(dataInBin)/area/(bin1-bin0))
        svec.append(numpy.mean(dataInBin))#XXX IS THERE A BETTER WAY TO PICK BIN CENTERS?
        bin0 = bin1
    return numpy.array(dnds), numpy.array(svec)

def dndsVecs(cats, cols, scales, bins, area):
    dndss = []
    svecs = []
    for i in range(len(cats)):
        dnds, svec = dndsVec(cats[i], cols[i], scales[i], bins, area)
        dndss.append(dnds)
        svecs.append(svec)
#     print "svecs", svecs
    svec = numpy.mean(svecs, axis=0)
    dndsMean = numpy.mean(dndss, axis=0)
    dndsRMS  = numpy.std(dndss, axis=0)
    return numpy.array(svec), dndsMean, dndsRMS



def integralPlot(cat , col, xlim, ylim, xlabel, ylabel, filename, scale = 1.0):
    cat.sortBy(col)
    cat.reverse()
    detected = 0.
    detectedVec = []
    total = 0.
    totalVec = []
    for row in cat:
        total += 1.
        totalVec.append(total)
        if row['matched'] == 'Y':
            detected += 1.
        detectedVec.append(detected)
    colData = cat.arrayFromCol(col)*scale
    frac = numpy.array(detectedVec)/numpy.array(totalVec)
    pylab.plot(colData, frac)
    pylab.xlim(xlim)
    pylab.ylim(ylim)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)
    pylab.savefig(filename)
    pylab.clf()

