from tools.dataparselib import DataParser

obj = DataParser()
obj.load_matrix("/home/robotbrunch/Desktop/20170717_all_robustness.dat")
early=(obj
  .get_zero_filtered()
  .only_columns(['std','early'])
  .sort_on_key('std')
)
stat=(obj
  .get_zero_filtered()
  .only_columns(['std','stat'])
  .sort_on_key('std')
)
n_stat=(obj
  .get_zero_filtered()
  .only_columns(['std','n_stat'])
  .sort_on_key('std')
)
dyn=(obj
  .get_zero_filtered()
  .only_columns(['std', 'dyn'])
  .sort_on_key('std')
)
dc=(obj
  .get_zero_filtered()
  .only_columns(['std','dc'])
  .sort_on_key('std')
)

early_groupings = early.get_groupings(['std'])
stat_groupings = stat.get_groupings(['std'])
n_stat_groupings = n_stat.get_groupings(['std'])
dyn_groupings = dyn.get_groupings(['std'])
dc_groupings = dc.get_groupings(['std'])

print("early")
for k in early_groupings:
  print(early_groupings[k].get_means().to_csv(include_headers=False)
    +","+str(early_groupings[k].get_90_confint('early')))
print("stat")
for k in stat_groupings:
  print(stat_groupings[k].get_means().to_csv(include_headers=False)
    +","+str(stat_groupings[k].get_90_confint('stat')))
print("n_stat")
for k in n_stat_groupings:
  print(n_stat_groupings[k].get_means().to_csv(include_headers=False)
    +","+str(n_stat_groupings[k].get_90_confint('n_stat')))
print("dyn")
for k in dyn_groupings:
  print(dyn_groupings[k].get_means().to_csv(include_headers=False)
    +","+str(dyn_groupings[k].get_90_confint('dyn')))
print("dc")
for k in dc_groupings:
  print(dc_groupings[k].get_means().to_csv(include_headers=False)
    +","+str(dc_groupings[k].get_90_confint('dc')))
