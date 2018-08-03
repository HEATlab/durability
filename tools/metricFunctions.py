import randomWalk
import test
import stnConverter
import sys

def sphericalMetric(S, dim, met):
	if met == "radRatio": return radRatio(S)
	elif met == "wilson" : return test.wilsonFlex(S)
	elif met == "spherical": return test.spherical(S, dim)
	elif met == "burgers": return test.hunsbergerFlex(S)
	else:
		return test.naiveFlex(S)

# Failed Metrics
def radRatio(inequalDict):
    chebRad = randomWalk.chebyshev(inequalDict)[1]
    cheb = randomWalk.chebyshev(inequalDict)[0]

    maxRad = -sys.maxint
    for l in inequalDict.values():
        for inequal in l:
            r = inequal[2] - cheb[inequal[1]] + cheb[inequal[0]]
            maxRad = max(r, maxRad)

	return chebRad / maxRad

def avgRadRatio(inequalDict):
    chebRad = randomWalk.chebyshev(inequalDict)[1]
    cheb = randomWalk.chebyshev(inequalDict)[0]

    radSum = 0.0
    count = 0
    for l in inequalDict.values():
        for inequal in l:
            r = inequal[2] - cheb[inequal[1]] + cheb[inequal[0]]
            radSum += r
            count += 1

    avgRad = radSum / count

    return chebRad / avgRad

def avgMinRadRatio(inequalDict):
    chebRad = randomWalk.chebyshev(inequalDict)[1]
    cheb = randomWalk.chebyshev(inequalDict)[0]

    radSum = 0.0
    for l in inequalDict.values():
    	currRad = sys.maxint
        for inequal in l:
            r = inequal[2] - cheb[inequal[1]] + cheb[inequal[0]]
            currRad = min(r, currRad)
        radSum += currRad

    avgRad = radSum / len(inequalDict.values())

    return chebRad/ avgRad

def sphereVol(S, dim):
    chebsph = randomWalk.chebyshev(S)
    return test.sphereVolume(float(chebsph[1]),dim)/(test.naiveVolFlex(stnConverter.dictToMatrix(S))+0.00000001)

