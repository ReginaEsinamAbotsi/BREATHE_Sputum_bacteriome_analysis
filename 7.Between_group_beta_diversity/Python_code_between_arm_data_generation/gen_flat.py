# %%
import os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import glob
import numpy as np
from itertools import combinations

meta = pd.read_csv('./meta_data_final31stAugust.csv', index_col=0)
meta = meta[['visit','trial_arm']]

PATHS_1 = glob.glob('Fla*/*/*bray1*')
PATHS_2 = glob.glob('Fla*/*/*CLR*')
PATHS_3 = glob.glob('new_data/*')

def gen_plots(PATHS, target='visit'):
	output_df = pd.DataFrame()

	for i in PATHS:
		df = pd.read_csv(i, index_col=0)
		name = i.split('/')[-1]
		meta_sel = meta.loc[df.index]
		print(meta_sel.columns.tolist())
		groups = meta_sel[target].unique()
		g1 = groups[0]
		g2 = groups[1]
	
		g1_idx = meta_sel[meta_sel[target]==g1].index.tolist()
		g2_idx = meta_sel[meta_sel[target]==g2].index.tolist()
		g3_idx = meta_sel.index.tolist()	
	
		df1 = df.loc[g1_idx, g1_idx]
		df2 = df.loc[g2_idx, g2_idx]
		df3 = df.loc[g3_idx, g3_idx]

		idx_3_dic = {k:v for k,v in zip(meta_sel.index.tolist(),meta_sel[target].tolist())}
		print(df3.sum().sum())
		for n, idx in enumerate(combinations(g3_idx, 2)):
			if idx_3_dic[idx[0]] == idx_3_dic[idx[1]]:
				df3.loc[idx[0],idx[1]] = 0
				df3.loc[idx[1],idx[0]] = 0
		print(df3.sum().sum())
	 
		ser1 = np.tril(df1).flatten()
		ser1 = ser1[ser1!=0]
		ser2 = np.tril(df2).flatten()
		ser2 = ser2[ser2!=0]
		ser3 = np.tril(df3).flatten()
		print(ser3.shape)
		ser3 = ser3[ser3!=0]
		print(ser3.shape)
		
		name1_lst = [g1 for x in range(len(ser1))]
		name2_lst = [g2 for x in range(len(ser2))]
		name3_lst = ['Bewtween_group_matrix' for x in range(len(ser3))]
	
		ser_total = list(ser1) + list(ser2) + list(ser3)
		name_total = name1_lst + name2_lst + name3_lst
	
		if name.startswith('azm'):
			arm_lst = ['AZM' for x in range(len(name_total))]
		elif name.startswith('Placebo'):
			arm_lst = ['Placebo' for x in range(len(name_total))]
		else:
			arm_lst = []
		#output_df = pd.DataFrame({'Group':name_total, 'Values':ser_total, 'Arm':arm_lst})
		output_df = pd.DataFrame({'Group':name_total, 'Values':ser_total})
	   # output_df = pd.concat([output_df, add_plot_df], axis=0)	
	
		output_df.sort_values('Group',ascending=True,inplace=True)
		sns.violinplot(x='Group', y='Values', data=output_df)
		plt.tight_layout()
		plt.savefig(name+'.png')
		plt.close()
		output_df.to_csv(name + '.csv')
	
		#azm_only_df = output_df[output_df.Arm=='AZM']
		#sns.violinplot(x='Visit', y='Values', data=azm_only_df)
		#plt.tight_layout()
		#plt.savefig(name+'_azm_only.png')
		#plt.close()
		#azm_only_df.to_csv(name+'_azm_only.csv')
	
	#	placebo_only_df = output_df[output_df.Arm=='Placebo']
		#sns.violinplot(x='Visit', y='Values', data=placebo_only_df)
		#plt.tight_layout()
		#plt.savefig(name+'_placebo_only.png')
		#plt.close()
		#placebo_only_df.to_csv(name+'_placebo_only.csv')

	return

#gen_plots(PATHS_1)
#gen_plots(PATHS_2)
gen_plots(PATHS_3, 'trial_arm')