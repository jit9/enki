import numpy

def linearFit( x, y, err, fixZero=False ):
    err = numpy.array(err)
    x = numpy.array(x)
    y = numpy.array(y)
    
    if fixZero:
        M = numpy.zeros([len(x),1])
        for i in xrange(len(x)):
            M[i] = x[i]
    else:
        M = numpy.ones([len(x),2])
        for i in xrange(len(x)):
            M[i][0] = x[i]
    Ninv = numpy.diag(1/err**2)
    cov = numpy.linalg.inv(numpy.dot(M.transpose(), numpy.dot(Ninv, M)))
    ans = numpy.dot(cov, numpy.dot(M.transpose(), numpy.dot(Ninv, y)))
    return ans, cov

