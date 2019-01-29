"""This program creates the Propagated SPOKE Entry Vectors
Usage:
  make_psevs_by_node_type.py [-a=N] [-b=N] [-w=N] [-t=NAME] [-i=NAME] [-d=NAME] [-n=NAME] [-m=NAME] [-p=NAME] [-e=NAME] [-s=NAME]
  make_psevs_by_node_type.py -h | --help

Options:
  -h --help     Show this screen.
  -a=N     Acceptable difference threshold [default: 0.001].
  -b=N     Probability of Random Jump [default: 0.1].
  -w=N     Number of workers for parallel processing [default: 6].
  -t=NAME     SPOKE node type [default: Gene].
  -i=NAME     Input directory [default: build/].
  -d=NAME     Filename containing direct hit list [default: all_patient_cohort_direct_hit_list.npy].
  -n=NAME     Filename containing SPOKE node list [default: all_patient_cohort_spoke_node_list.npy].
  -m=NAME     Filename containing patient direct hit matrix [default: all_patient_cohort_binary_direct_hit_matrix].
  -p=NAME     Filename transition probability matrix [default: transition_probability_matrix.npy].
  -e=NAME     Filename spoke edges [default: neo4j_edges.tsv].
  -s=NAME     SPOKE directory [default: spoke_v_1/].
"""
from docopt import docopt

import pickle as pic
import operator
import numpy as np
import re
import os
import multiprocessing as mp

arguments = docopt(__doc__)

# Load input arguments
acceptable_diff=float(arguments['-a'])
probability_random_jump=float(arguments['-b'])
num_workers=int(arguments['-w'])
direct_hit_list_filename = arguments['-d']
node_list_filename = arguments['-n']
patient_matrix_filename = arguments['-m']
transition_prob_matrix_filename = arguments['-p']
spoke_edge_filename = arguments['-e']
input_directory = arguments['-i']
spoke_node_type = arguments['-t']
spoke_directory = arguments['-s']

# Load ordered list of direct hit nodes
direct_hit_list = np.load(input_directory+direct_hit_list_filename, allow_pickle=False)
# Of direct hits
number_direct_hits = len(direct_hit_list)

# Load ordered list of all SPOKE
node_list = np.load(input_directory+node_list_filename, allow_pickle=False)
# Total number SPOKE nodes
number_nodes = len(node_list)
# All SPOKE nodes id to index dict
node_to_index_dict = dict(zip(node_list, range(number_nodes)))

# Load node type list 
node_type_list = np.load(spoke_directory+'node_type_list.npy', allow_pickle=False)
# Make direct hit node type list
direct_hit_node_type_list = np.array([node_type_list[node_to_index_dict[direct_hit]] for direct_hit in direct_hit_list])
# Make index of direct hits in entire spoke node list
direct_hit_index_list = np.array([node_to_index_dict[direct_hit] for direct_hit in direct_hit_list])

def load_neo4j_connect_matrix(node_to_index_dict):
	'''
	Loads a symmetric connectivity matrix
	If there is an edge between two SPOKE nodes then the cell will be 1 else 0
	'''
	connectivity_matrix = np.zeros((number_nodes, number_nodes))
	with open(spoke_directory+spoke_edge_filename) as edge_file:
		for lines in edge_file:
			node_1, relationship, node_2 = lines.strip('\n').split('\t')
			connectivity_matrix[node_to_index_dict[node_1]][node_to_index_dict[node_2]] = 1
	return connectivity_matrix

def make_all_direct_hit_arrays():
	'''
	Make transition probability matrix using all direct hits 
	'''
	matrix = np.zeros((number_direct_hits,number_direct_hits))
	with open(input_directory+patient_matrix_filename) as sig_file:
		for lines in sig_file:
			direct_hit = np.array(lines.strip('\n').split(), dtype=float) == 1
			add_val = np.nan_to_num(direct_hit.astype(float)/np.sum(direct_hit))
			matrix[direct_hit ,] = matrix[direct_hit ,] + add_val
	matrix = np.transpose(matrix)
	matrix = np.transpose(np.nan_to_num(matrix/np.sum(matrix,axis=0)))
	return matrix

# Make transition probability matrix
direct_hit_array_list = np.zeros((2,2))
if os.path.exists(input_directory+transition_prob_matrix_filename):
	direct_hit_array_list = np.load(input_directory+transition_prob_matrix_filename, allow_pickle=False)
else:
	direct_hit_array_list = make_all_direct_hit_arrays()
	np.save(input_directory+transition_prob_matrix_filename, direct_hit_array_list, allow_pickle=False)


# Find index for desired direct hits
desired_direct_hit_index = direct_hit_node_type_list == spoke_node_type
# Narrow down rows to only desired direct hit
direct_hit_array_list = direct_hit_array_list[np.arange(number_direct_hits)[desired_direct_hit_index]]
# Make direct hit list of only desired type
direct_hit_list = direct_hit_list[desired_direct_hit_index]
# Update number of direct hits
number_direct_hits = len(direct_hit_list)

# Load SPOKE connectivity matrix
connectivity_matrix = load_neo4j_connect_matrix(node_to_index_dict)

# Make connectivity transision probability
connectivity_matrix = np.transpose(np.nan_to_num(connectivity_matrix/np.sum(connectivity_matrix, axis=0)))

def get_rank_vector(index_and_probability_random_jump):
	index, probability_random_jump = index_and_probability_random_jump
	full_direct_hit_array = np.zeros(number_nodes)
	full_direct_hit_array[direct_hit_index_list] = direct_hit_array_list[index]
	connect_matrix = np.transpose((full_direct_hit_array*probability_random_jump)+(connectivity_matrix*(1 - probability_random_jump)))
	rank_vector = np.ones(number_nodes)/float(number_nodes)
	rank_diff, i = 1, 0
	while (rank_diff > acceptable_diff) and (i < 40):
		new_rank_vector = np.dot(connect_matrix, rank_vector)
		rank_diff = np.sum(np.abs(new_rank_vector - rank_vector))
		rank_vector = new_rank_vector
		i+= 1
	return rank_vector/np.sum(rank_vector)

def make_all_direct_hit_page_rank_parallel(probability_random_jump):
	p = mp.Pool(num_workers)
	index_and_probability_random_jump = zip(range(number_direct_hits), np.full(number_direct_hits, probability_random_jump))
	rank_vector_list_of_lists = np.array(p.map(get_rank_vector, index_and_probability_random_jump))
	p.close()
	p.join()
	return rank_vector_list_of_lists


rank_vector_list_of_lists = make_all_direct_hit_page_rank_parallel(probability_random_jump)
filename = spoke_node_type+'_direct_hit_page_rank_E_' + '_'.join(str(probability_random_jump).split('.')) + '_A_' + '_'.join(str(acceptable_diff).split('.'))
np.save(input_directory+filename, rank_vector_list_of_lists, allow_pickle=False)

