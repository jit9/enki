from flipper import *
from optimalFilter import *
import cluster
import catalog

paramFile = sys.argv[1]
params = flipperDict.flipperDict(paramFile)
params.readFromFile(paramFile)

def printTime( then ):
    now = time.time()
    print "Time Elapsed = %s sec." % (now - then)
    return now

print "Initializing..."
then = time.time()

data = "%s/%s" % (params['inputDir'], params['data'])
sub1 = "%s/%s" % (params['inputDir'], params['sub1'])
sub2 = "%s/%s" % (params['inputDir'], params['sub2'])
noise = "%s/%s" % (params['inputDir'], params['noise'])
weight = "%s/%s" % (params['inputDir'], params['weight'])
dataMap = liteMap.liteMapFromFits(data)

if not os.path.exists(params['outputDir']):
    os.makedirs(params['outputDir'])


def filterAndExtract( pickleName, filteredMapName, sourceType, threshold, then, fwhm = 1.4, filterMask = None, detectionMask=None, force = False):
    if not os.path.exists(pickleName) and not force: 
        if not os.path.exists( filteredMapName ): 
            print "Filtering to find %s" % sourceType
            filteredMap = optimalFilter(data, sub1, sub2, noise, weight, fwhm = fwhm, masking = filterMask)
            filteredMap.writeFits( filteredMapName, overWrite=True)
        else:
            print "Found filtered %s map at %s" % (sourceType, filteredMapName)
            filteredMap = liteMap.liteMapFromFits( filteredMapName )
        then = printTime(then)
        print "Extracting %s" % sourceType
        if detectionMask != None:
            for m in detectionMask:
                x1, y0 = filteredMap.skyToPix(m[0], m[2])
                x0, y1 = filteredMap.skyToPix(m[1], m[3])
                filteredMap.data[y0:y1, x0:x1] = 0.
        cf = cluster.clusterFinder(filteredMap, sourceType, threshold, )
        cf.findClusters()
        cf.setRADecOfClusters()
        cf.clusters.save(pickleName)
        clusters = cf.clusters
    else:
        print "Found previous %s pickle" % sourceType
        clusters = cluster.clusterList()
        clusters.read(pickleName)
    return clusters , then

then = printTime(then)
sourcePickle = "%s/sources.pickle" % params['outputDir']
fm1 = "%s/%s" % (params["outputDir"], params["filteredMap1"])
sources, then = filterAndExtract( sourcePickle, fm1, "sources", params['sourceThreshold'], 
        then, force = params['forceSourceExtraction'] )

if params['saveSourceCatalog']:
#     cf.catalog(clusters, dataMap, 3, outputDir= params["outputDir"], getAltName = params["getAltName"])
    cat = catalog.catalogFromClusterList(sources)
    cat.sortBy('s/n')
    cat.reverse()
    cat.write("%s/%s" % (params["outputDir"], params["sourceCatalog"]))

then = printTime(then)
# clusterPickle = "%s/clusters.pickle" % params['outputDir']
# fm2 = "%s/%s" % (params["outputDir"], params["filteredMap2"])
# clusters, then = filterAndExtract( clusterPickle, fm2, "clusters", params['clusterThreshold'], then, 
#         detectionMask = params['zeroMask'], force = params['forceClusterExtraction'], filterMask = sources, fwhm=params['clusterFWHM'] )
# 
# then = printTime(then)
# if params['saveClusterCatalog']:
#     cat = catalog.catalogFromClusterList(clusters)
#     cat.sortBy('s/n')
#     cat.reverse()
#     cat.write("%s/%s" % (params["outputDir"], params["clusterCatalog"]))
#     cat.writeRegionFile("%s/%s" % (params["outputDir"], params["clusterRegionFile"]))
