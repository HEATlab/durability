# This file contains functions that convert robot brunch STNs to our format
# for different STNs
# Authors: Viva Ojha and Joon Lee
# Contact Info: vmojha@g.hmc.edu and joolee@g.hmc.edu

import os
import glob
import json
from stntools import stn, loadSTNfromJSONfile
from randomWalk import getInequal


the_path = '../stngen/old_stn_gen/'


for i in range(1,5):
    the_stn_list = []
    curr_path = os.path.abspath(the_path + 'normal_run' + str(i) + '/*/*.json')
    file_list = glob.glob(curr_path)
    for file in file_list:
        #print str(folder)
        #print glob.glob(os.path.abspath(str(folder) + '*.json'))
        #for file in glob.glob(folder):
        ineqDict = getInequal(loadSTNfromJSONfile(file,0)['stn'])
        # newIneqDict = {}
        # for key in ineqDict.keys():
        #     newIneqDict[int(key)] = []
        #     val = ineqDict[key]
        #     for ineq in val:
        #         if ineq[2] != float('inf'):
        #             newIneq = []
        #             newIneq.append(ineq[0])
        #             newIneq.append(ineq[1])
        #             newIneq.append(ineq[2]/1000.0)
        #             newIneqDict[int(key)].append(newIneq)
        #print newIneqDict
        the_stn_list.append(ineqDict)
    with open(os.path.abspath('./Data/' + 'the_many_many_stns_' + str(i) +
              '.json'), 'w') as fp:
        json.dump(the_stn_list, fp)


# '''
# ineqDict = getInequal(loadSTNfromJSONfile(the_path +
#                       'normal_run4/STN_a2_i4_s1_t2000/minimal_3.json')['stn'])
# newIneqDict = {}
# for key in ineqDict.keys():
#     newIneqDict[key] = []
#     val = ineqDict[key]
#     for ineq in val:
#         newIneq = []
#         newIneq.append(ineq[0])
#         newIneq.append(ineq[1])
#         newIneq.append(ineq[2]/1000.0)
#         newIneqDict[key].append(newIneq)
# print newIneqDict
# the_stn_list.append(newIneqDict)
# '''

# file_path = os.path.abspath('../handmade_stns/dc_test.json')

# ineqDict = getInequal(loadSTNfromJSONfile(file_path, 0)['stn'])

# '''
# newIneqDict = {}
# for key in ineqDict.keys():
#     newIneqDict[int(key)] = []
#     val = ineqDict[key]
#     for ineq in val:
#         if ineq[2] != float('inf'):
#             newIneq = []
#             newIneq.append(ineq[0])
#             newIneq.append(ineq[1])
#             newIneq.append(ineq[2]/1000.0)
#             newIneqDict[int(key)].append(newIneq)
# '''
# # print newIneqDict
# # the_stn_list.append(newIneqDict)
# with open(os.path.abspath('./Data/dc_test_ineq.json'), 'w') as fp:
#     json.dump(ineqDict, fp)
