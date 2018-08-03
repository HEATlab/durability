import json
import os
import glob



def getIneq4TU(f):
	ineqDict = {}
	
	for line in f:
		if not line.find('specifies each arc that should be considered a refinement constraint', 0, len(line)) > 0:
			if line.find('<refine>', 0, len(line)) > 0:
				values = line.split()
				start = int(values[2]) - 1
				end = int(values[3]) - 1
				constraint = int(values[4])
				
				ineq = [start, end, constraint]

				if start not in ineqDict:
					ineqDict[start] = []
				if end not in ineqDict:
					ineqDict[end] = []
				
				ineqDict[start].append(ineq)
				ineqDict[end].append(ineq)
				
	return ineqDict

the_path = './Data/Small STNs'
the_stn_list = []
file_list = glob.glob(os.path.abspath(the_path + '/*.txt'))

for file in file_list:
    f = open(file, 'r')
    the_stn_list.append(getIneq4TU(f))

with open(os.path.abspath('./Data/dataset.json'), 'w') as fp:
    json.dump(the_stn_list, fp)

# ineqDict4TU = getIneq4TU(f)

# for i in range(41):
# 	print str(i) + ': ' 
# 	print ineqDict4TU[i]
# 	print '\n'

# with open(os.path.abspath('2_10_50_50_0.json'), 'w') as fp:
#     json.dump(ineqDict4TU, fp)
