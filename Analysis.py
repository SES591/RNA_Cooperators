import numpy as np
import matplotlib.pylab as plt
import random
import Parameters
import pandas as pd
#########################################################################################################
def plot_all_assembilies(num_runs):
	
	for exp in range(num_runs):
		fig  = plt.figure()
		for i in range(3):
			print i
			ax = fig.add_subplot(3,2,i+1)
			fname = Parameters.dirname + '/%i_assembly_%i.dat' % (exp, i)
			t, pop = np.loadtxt(fname, unpack = True)
			ax.plot(t[:-5],pop[:-5])
		plt.savefig(Parameters.dirname + '/%i_all_assembilies.png' %exp)
		plt.close()
#########################################################################################################
def plot_selfish_cooperative(num_runs):
	import seaborn as sns
	for exp in range(num_runs):
		fig  = plt.figure()
		cooperators = [3,4,5]
		selfish = [0,1,2]

		
		fname = Parameters.dirname + '/%i_assembly_%i.dat' % (exp, 0)
		t, pop = np.loadtxt(fname, unpack = True)

		cooperative_pop = [0.0]*len(pop[:-5])
		selfish_pop = [0.0]*len(pop[:-5])

		for c in cooperators:
			fname = Parameters.dirname + '/%i_assembly_%i.dat' % (exp, c)
			t, pop = np.loadtxt(fname, unpack = True)
			cooperative_pop = np.add(cooperative_pop,pop[:-5] )

		for s in selfish:
			fname = Parameters.dirname + '/%i_assembly_%i.dat' % (exp, s)
			t, pop = np.loadtxt(fname, unpack = True)
			selfish_pop = np.add(selfish_pop,pop[:-5] )

		ax = fig.add_subplot(1,1,1)
		ax.plot(t[:-5], selfish_pop, label = 'selfish', color = 'r')
		ax.plot(t[:-5], cooperative_pop, label = 'cooperative', color = 'g')
		ax.legend(loc = 'upper left')
		plt.xlabel('System Time')
		plt.ylabel('Total Abundance')

		#plt.show()
		plt.savefig(Parameters.dirname + '/%i_cooperative_vs_selfish.png' % exp)
		plt.close()
#########################################################################################################
def pickle_time_series_data(num_runs):
	import pickle

	seq_names =['AA', 'CU', 'UA', 'AU', 'CG', 'UG']
	timeseries = {'UA': {}, 'AU': {}, 'CG': {}, 'UG' : {},'AA': {},'CU': {} }
	for exp in range(0,num_runs):

		filename = ('%s/%i_rank_coarse_assemblies.csv' % (Parameters.dirname, exp))
		coarse_df = pd.read_csv(filename, header = 0, index_col= False)
		for seq in seq_names:
			time = coarse_df[seq].tolist()
			timeseries[seq][exp] = time
	fname = 'initfrags_%i_ksexp_%i.pickle' % (Parameters.init_frags, Parameters.ksexp)

	with open(fname, 'wb') as handle:
  		pickle.dump(timeseries, handle)
		
#########################################################################################################
def gen_ka_df():
	from Parameters import ka_dict
	import pandas as pd
	import matplotlib.colors as colors
	import seaborn as sns
	seq_names = ['UA','AU','CG','UG','AA','CU']
	igss = ['U', 'A', 'C', 'U']
	tagss = ['A', 'U', 'G']

	a = np.zeros(shape=(len(seq_names),len(seq_names)))
	ka_df =pd.DataFrame(a, index=seq_names, columns=seq_names)

	for seq1 in seq_names:
		for seq2 in seq_names:
			igs = seq1[0]
			tag = seq2[1]
			ka_df[seq1].loc[seq2] = ka_dict[(igs, tag)]

	sns.heatmap(ka_df, annot = True) #, norm = colors.SymLogNorm(linthresh=0.001, linscale=0.05)
	plt.show()
	ka_df.to_pickle('ka_df2.pickle')

#########################################################################################################	
#plot_all_assembilies(10)
plot_selfish_cooperative(10)
#pickle_time_series_data(10)
#gen_ka_df()
