from tools.dataparselib import DataParser

obj = DataParser()
#obj.get_matrix("/home/robotbrunch/Desktop/test.dat")
obj.load_matrix("/home/robotbrunch/Desktop/20170717_all_robustness.dat")
early=(obj
  .get_zero_filtered()
  .only_columns(['tght','early'])
  .sort_on_key('tght')
)
stat=(obj
  .get_zero_filtered()
  .only_columns(['tght','stat'])
  .sort_on_key('tght')
)
n_stat=(obj
  .get_zero_filtered()
  .only_columns(['tght','n_stat'])
  .sort_on_key('tght')
)
dyn=(obj
  .get_zero_filtered()
  .only_columns(['tght', 'dyn'])
  .sort_on_key('tght')
)
dc=(obj
  .get_zero_filtered()
  .only_columns(['tght','dc'])
  .sort_on_key('tght')
)

early_groupings = early.get_groupings(['tght'])
stat_groupings = stat.get_groupings(['tght'])
n_stat_groupings = n_stat.get_groupings(['tght'])
dyn_groupings = dyn.get_groupings(['tght'])
dc_groupings = dc.get_groupings(['tght'])

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
