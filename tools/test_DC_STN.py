from dc_stn import *
import sys

S = loadSTNfromJSONfile("json/out.json")['stn']
print S
S_copy = S.copy()
if not S_copy.floydWarshall():
  print "Input STN is not consistent"
  sys.exit(0)
elif not S.isStronglyControllable(debug=True):
    print "Input STN is not strongly controllable"
    sys.exit(0)
else:
  print "Input STN is consistent\n"

S_DC = STNtoDCSTN(S)

print S_DC
dc = S_DC.is_DC(debug_flag=True)
print "\nOutput of is_DC: {}\n".format(dc)
print S_DC
