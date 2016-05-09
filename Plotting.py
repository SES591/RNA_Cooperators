import numpy as np
import matplotlib.pylab as plt
import Parameters
######################################################################################################################################
def output_plots(exp, sequences, catalysts, substrates):
	"""Use this function to control which plots are generated at the end of the experimental run """
	plot_monomers(exp)
	plot_survivors(exp)
	plot_all(exp)
	if Parameters.hypercycle == True:
		plot_network_polymers(exp, sequences, catalysts, substrates)

######################################################################################################################################
def plot_monomers(exp):
	""" Monomer populations as time series, coloring is based on rgba, maps to network visualization """
	t, Ws = np.loadtxt(Parameters.dirname+'/%i_monomer_W.dat' % exp, unpack = True)
	t, Xs = np.loadtxt(Parameters.dirname+'/%i_monomer_X.dat' % exp, unpack = True)
	t, Ys = np.loadtxt(Parameters.dirname+'/%i_monomer_Y.dat' % exp, unpack = True)
	t, Zs = np.loadtxt(Parameters.dirname+'/%i_monomer_Z.dat' % exp, unpack = True)

	plt.plot(t,Ws, label = 'W')
	plt.plot(t,Xs, label = 'X')
	plt.plot(t,Ys, label = 'Y')
	plt.plot(t,Zs, label = 'Z')
	plt.legend()
	plt.savefig(Parameters.dirname+'/%i_monomers.png'%exp)
	plt.close()
######################################################################################################################################
def plot_network_polymers(exp, sequences, catalysts, substrates):

	for i in range(len(Parameters.hypercycles)):
		Polymers = Parameters.hypercycles[i]
		for j in range(len(Polymers)):
			ID = Polymers[j]
			label = sequences[ID].sequence
			t,population = np.loadtxt(Parameters.dirname+'/%i_sequence_%i.dat' % (exp,ID), unpack = True)
			plt.plot(t, population, label= label)

		plt.legend()
		plt.savefig(Parameters.dirname+'/%i_hypercycle_population_%i.png' % (exp,i))
		plt.close()
######################################################################################################################################
def plot_survivors(exp):

	survivors, IDs = np.loadtxt('%s/%i_surviving_species.dat' % (Parameters.dirname, exp), dtype = str, unpack = True)
	w = 0.0
	x = 0.0
	y = 0.0
	z = 0.0

	for i in range(len(survivors)):
		ID = int(IDs[i])+1
		label = survivors[i]
		t,population = np.loadtxt(Parameters.dirname+'/%i_sequence_pops.dat' % exp, unpack = True, usecols = (0, ID))
		plt.plot(t, population, label= label)
		final_pop = population[-1]
		for m in range(len(label)):
			if label[m] == "W":
				w += final_pop
			elif label[m] == "X":
				x += final_pop
			elif label[m] == "Y":
				y += final_pop
			elif label[m] == "Z":
				z += final_pop
	title = "Total Composition: W = %i, X=%i, Y=%i, Z=%i" % (w,x,y,z)
	plt.title(title)
	plt.legend()
	plt.savefig(Parameters.dirname+'/%i_survivors.png' % (exp))
	plt.close()
######################################################################################################################################
def plot_all(exp):

	IDs, seqs = np.loadtxt('%s/%i_seq_dictionary.dat' % (Parameters.dirname, exp), dtype = str, unpack = True)
	w = 0.0
	x = 0.0
	y = 0.0
	z = 0.0

	for i in range(1,len(IDs)):
		ID = int(IDs[i])
		label = seqs[i]
		t,population = np.loadtxt(Parameters.dirname+'/%i_sequence_pops.dat' % exp, unpack = True, usecols = (0, ID))
		final_pop = population[-1]
		mean_pop = np.mean(population)
		if mean_pop != 0.0:
			plt.plot(t, population, label= label)


		for m in range(len(label)):
			if label[m] == "W":
				w += final_pop
			elif label[m] == "X":
				x += final_pop
			elif label[m] == "Y":
				y += final_pop
			elif label[m] == "Z":
				z += final_pop
	title = "Total Composition: W = %i, X=%i, Y=%i, Z=%i" % (w,x,y,z)
	plt.title(title)
	plt.legend()
	plt.savefig(Parameters.dirname+'/%i_all_seqs.png' % (exp))
	plt.close()
######################################################################################################################################
def plot_network(exp, sequences, catalysts, substrates):
	import networkx as nx
	import random
	network = nx.DiGraph()
	edges = []
	node_colors = []
	edge_colors = []
	labels = {}

	for i in range(len(catalysts)):
		cat_ID = catalysts[i]
		sub_ID = substrates[i]

		cat_label = sequences[cat_ID].sequence
		sub_label = sequences[sub_ID].sequence

		cat_color = sequences[cat_ID].seq_list
		cat_color = [float(x)/4.0 for x in cat_color]
		cat_color = tuple(cat_color)
		sub_color = sequences[sub_ID].seq_list
		sub_color = [float(x)/4.0 for x in sub_color]
		sub_color = tuple(sub_color)

		node_colors.append(cat_color)
		node_colors.append(sub_color)

		if cat_ID in Parameters.reps:
			edge_colors.append('g')
		elif cat_ID in Parameters.recs:
			edge_colors.append('r')

		network.add_node(cat_ID, label = cat_label, node_color = cat_color)
		labels[cat_ID] = cat_label
		network.add_node(sub_ID, label = sub_label, node_color = sub_color)
		labels[sub_ID] = sub_label

		edge = (catalysts[i], substrates[i])
		edges.append(edge)



	network.add_edges_from(edges)

	nx.draw(network, node_color = node_colors, edge_color = edge_colors, with_labels = True, labels = labels, alpha = 0.6, node_size = 1500)

	plt.savefig(Parameters.dirname+'/%i_Network_visualization.png' %exp)
	plt.close()
