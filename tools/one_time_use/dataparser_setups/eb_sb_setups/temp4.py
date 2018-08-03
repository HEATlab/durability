from tools.dataparselib import DataParser

obj = DataParser()
obj.load_matrix("/home/robotbrunch/Desktop/20170717_all_robustness.dat")

X="ia"

srea_better=(obj
  .get_zero_filtered()
  .get_filtered_less('early','stat')
)

early_better=(obj
  .get_zero_filtered()
  .get_filtered_less('stat','early')
)

# Split into the proper columns

sb_early=(srea_better
  .only_columns([X,'early'])
  .sort_on_key(X)
)
eb_early=(early_better
  .only_columns([X,'early'])
  .sort_on_key(X)
)

sb_nstat=(srea_better
  .only_columns([X,'n_stat'])
  .sort_on_key(X)
)
eb_nstat=(early_better
  .only_columns([X,'n_stat'])
  .sort_on_key(X)
)

sb_dyn=(srea_better
  .only_columns([X,'dyn'])
  .sort_on_key(X)
)
eb_dyn=(early_better
  .only_columns([X,'dyn'])
  .sort_on_key(X)
)

sb_dc=(srea_better
  .only_columns([X,'dc'])
  .sort_on_key(X)
)
eb_dc=(early_better
  .only_columns([X,'dc'])
  .sort_on_key(X)
)

# Get groupings

sb_early_gs = sb_early.get_groupings([X])
eb_early_gs = eb_early.get_groupings([X])

sb_nstat_gs = sb_nstat.get_groupings([X])
eb_nstat_gs = eb_nstat.get_groupings([X])

sb_dyn_gs = sb_dyn.get_groupings([X])
eb_dyn_gs = eb_dyn.get_groupings([X])

sb_dc_gs = sb_dc.get_groupings([X])
eb_dc_gs = eb_dc.get_groupings([X])


# Print out

print("sb_early")
for k in sb_early_gs:
  grp = sb_early_gs[k]
  print(grp.get_means().to_csv(include_headers=False)
    +","+str(grp.get_90_confint('early')))
print("eb_early")
for k in eb_early_gs:
  grp = eb_early_gs[k]
  print(grp.get_means().to_csv(include_headers=False)
    +","+str(grp.get_90_confint('early')))
print("-----------------------")
print("sb_nstat")
for k in sb_early_gs:
  grp = sb_nstat_gs[k]
  print(grp.get_means().to_csv(include_headers=False)
    +","+str(grp.get_90_confint('n_stat')))
print("eb_nstat")
for k in eb_nstat_gs:
  grp = eb_nstat_gs[k]
  print(grp.get_means().to_csv(include_headers=False)
    +","+str(grp.get_90_confint('n_stat')))
print("-----------------------")
print("sb_dyn")
for k in sb_dyn_gs:
  grp = sb_dyn_gs[k]
  print(grp.get_means().to_csv(include_headers=False)
    +","+str(grp.get_90_confint('dyn')))
print("eb_dyn")
for k in eb_dyn_gs:
  grp = eb_dyn_gs[k]
  print(grp.get_means().to_csv(include_headers=False)
    +","+str(grp.get_90_confint('dyn')))
print("-----------------------")
print("sb_dc")
for k in sb_dc_gs:
  grp = sb_dc_gs[k]
  print(grp.get_means().to_csv(include_headers=False)
    +","+str(grp.get_90_confint('dc')))
print("eb_dc")
for k in eb_dc_gs:
  grp = eb_dc_gs[k]
  print(grp.get_means().to_csv(include_headers=False)
    +","+str(grp.get_90_confint('dc')))
