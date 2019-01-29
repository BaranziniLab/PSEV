"""Connects all overlapping concepts between EHRs and SPOKE
Usage:
  connect_all_patients.py [-i=NAME] [-s=NAME] [-o=NAME]
  connect_all_patients.py -h | --help

Options:
  -h --help     Show this screen.
  -i=NAME     Input directory [default: build_4/].
  -s=NAME     SPOKE directory [default: spoke_v_1/].
  -o=NAME     Output filename [default: all_patient_cohort].
"""
from docopt import docopt

import pickle as pic
import operator
import numpy as np
from collections import Counter
import time
import os
import pandas as pd

arguments = docopt(__doc__)

directory = arguments['-i']
spoke_directory = arguments['-s']
save_filename = arguments['-o']

node_list = np.load(spoke_directory+'spoke_node_list.npy', allow_pickle=False)
node_type_list = np.load(spoke_directory+'node_type_list.npy', allow_pickle=False)
#
doid_to_icd_filename = 'direct_doid_to_icd_9_10_group_filled_in_final_used_for_query.tsv'
icd_9_cui_filename = 'ehr_cui_icd9.tsv'
icd_10_cui_filename = 'ehr_cui_icd10.tsv'
symptom_cui_and_mesh_filename = 'symptom_cui_and_mesh.tsv'
db_id_to_clarity_file = 'db_id_to_clarity.pkl'
loinc_filename='distinct_cui_loinc.tsv'
path_filename='cui_to_loinc_path.tsv'
cui_cpt_filename='distinct_cui_cpt.tsv'
#
patient_filename = directory+'patient_cohort.tsv'
diagnoses_table_filename = directory+'distinct_patient_icd_9_10_by_year.tsv'
medication_table_filename = directory+'distinct_patient_meds_by_year.tsv'
lab_table_filename = directory+'distinct_patient_lab_by_year.tsv'

def get_icd_to_hetio_disease(filename):
	disease_query_table = []
	with open(directory+'constants/'+filename, 'rU') as doid_file:
		header = doid_file.readline()
		for lines in doid_file:
			doid, disease_name, icd_9_list, icd_10_list = lines.strip('\n').split('\t')
			for icd in icd_9_list.split()+icd_10_list.split():
				disease_query_table.append([doid, icd])
	disease_query_table = pd.DataFrame(data=np.array(disease_query_table), columns=['SPOKE', 'EHR'])
	return disease_query_table

def get_icd_cui_mesh_tables():
	disease_query_table = get_icd_to_hetio_disease(doid_to_icd_filename)
	icd9_table = pd.read_table(directory+'constants/'+icd_9_cui_filename, sep='\t', header=0)
	icd10_table = pd.read_table(directory+'constants/'+icd_10_cui_filename, sep='\t', header=0)
	mesh_table = pd.read_table(directory+'constants/'+symptom_cui_and_mesh_filename, sep='\t', header=0)
	cui_icd_table = pd.concat((icd9_table, icd10_table), ignore_index=True).rename(index=str, columns={"CODE": "EHR", "CUI":"SPOKE"}).drop_duplicates()
	mesh_cui_icd_table = pd.merge(mesh_table, cui_icd_table, left_on='CUI', right_on='SPOKE').drop(['CUI', 'SPOKE'], axis=1).rename(index=str, columns={"CODE":"SPOKE"}).drop_duplicates()
	icd_to_spoke = pd.concat((disease_query_table, cui_icd_table, mesh_cui_icd_table), ignore_index=True)
	icd_to_spoke = icd_to_spoke[icd_to_spoke['SPOKE'].isin(np.array(node_list, dtype=str))]
	return icd_to_spoke

def get_clarity_to_dbid(db_to_clarity_id_dict):
	clarity_to_spoke = []
	for dbid, clarity_list in db_to_clarity_id_dict.items():
		for clarity_id in clarity_list:
			clarity_to_spoke.append([dbid, clarity_id])
	clarity_to_spoke = pd.DataFrame(data=np.array(clarity_to_spoke), columns=['SPOKE', 'EHR'])
	clarity_to_spoke['EHR'] = np.array(clarity_to_spoke['EHR'], dtype=str)
	return clarity_to_spoke

def get_loinc_map_stats(path_df):
	id_and_type = path_df.drop(['EHR'], axis=1).drop_duplicates()
	print '# of LOINC mapped to SPOKE:', path_df.EHR.unique().shape[0]
	print 'Mapped SPOKE breakdown:', Counter(id_and_type.Node_Type.values)

def filter_transform_labs_to_spoke(path_df, path_len_filter=2):
	cui_to_all_spoke = pd.read_table(spoke_directory+'All_SPOKE_CUI.tsv', sep='\t', header=0)
	cui_cpt_df = pd.read_table(directory+'constants/'+cui_cpt_filename, sep='\t', header=0)
	cui_loinc_df = pd.read_table(directory+'constants/'+loinc_filename, sep='\t', header=0)
	path_len_list, first_list, second_list = [], [], []
	for cui, path in path_df.values:
		path_list = path.split('->')
		path_len_list.append(len(path_list))
		first_list.append(path_list[0])
		second_list.append(path_list[1])
	path_df['Path_Len'] = np.array(path_len_list)
	path_df['First'] = np.array(first_list)
	path_df['Second'] = np.array(second_list)
	path_df = path_df.rename(index=str, columns={"CUI":"CUI_SPOKE"}).drop(['Path'], axis=1)
	path_df = path_df[path_df['Path_Len'] <= (path_len_filter+1)]
	cpt_path = path_df[path_df['Second'].isin(cui_cpt_df.CUI.unique())]
	path_df = pd.concat((path_df[path_df['Path_Len'] <= path_len_filter], cpt_path), ignore_index=True)
	path_df = pd.merge(path_df, cui_loinc_df, left_on='First', right_on='CUI').drop(['First', 'Second', 'Path_Len', 'CUI'], axis=1)
	path_df = pd.merge(cui_to_all_spoke, path_df, left_on='CUI', right_on='CUI_SPOKE').drop(['CUI', 'CUI_SPOKE'], axis=1).rename(index=str, columns={'Node_ID':'SPOKE', 'CODE':'EHR'}).drop_duplicates()
	get_loinc_map_stats(path_df)
	path_df = path_df.drop(['Node_Type'], axis=1).drop_duplicates()
	return path_df

def get_patient_to_spoke_df(icd_to_spoke, clarity_to_spoke, path_df):
	patient_diag_ehr = pd.read_table(diagnoses_table_filename, sep='\t', header=0)
	patient_diag_ehr = pd.concat((pd.merge(patient_diag_ehr.drop(['ICD10_Code'], axis=1), icd_to_spoke, left_on='ICD9_Code', right_on='EHR').drop(['EHR', 'ICD9_Code'], axis=1), pd.merge(patient_diag_ehr.drop(['ICD9_Code'], axis=1), icd_to_spoke, left_on='ICD10_Code', right_on='EHR').drop(['EHR', 'ICD10_Code'], axis=1)), ignore_index=True).drop_duplicates()
	#
	patient_med_ehr  = pd.read_table(medication_table_filename, sep='\t', header=0)
	patient_med_ehr['Medication_ID'] = np.array(patient_med_ehr['Medication_ID'], dtype=str)
	patient_med_ehr = pd.merge(patient_med_ehr, clarity_to_spoke, left_on='Medication_ID', right_on='EHR').drop(['EHR', 'Medication_ID'], axis=1).drop_duplicates()
	#
	patient_lab_ehr = pd.merge(pd.read_table(lab_table_filename, sep='\t', header=0), path_df, left_on='LOINC_Code', right_on='EHR').drop(['EHR', 'LOINC_Code'], axis=1).drop_duplicates()
	patient_lab_ehr = patient_lab_ehr[patient_lab_ehr['Lab_Result_Abnormal']=='Yes']
	patient_lab_ehr = patient_lab_ehr.drop(['Lab_Result_Abnormal'], axis=1)
	return pd.concat((patient_diag_ehr, patient_med_ehr, patient_lab_ehr), ignore_index=True)

def print_binary_matrix(patient_to_hetio_connect):
	# print patient matrix 
	output = open(directory+save_filename+'_binary_direct_hit_matrix', 'w')
	for row in patient_to_hetio_connect:
		output.write(' '.join(np.array(row, dtype=str))+'\n')
	output.close()


patient_list = np.array(pd.read_table(patient_filename, sep='\t', header=0)['Patient_ID'].values, dtype=str)
icd_to_spoke = get_icd_cui_mesh_tables()
clarity_to_spoke= get_clarity_to_dbid(pic.load(open('db_id_to_clarity.pkl', 'r')))
#
path_df = pd.read_table(directory+'constants/'+path_filename, sep='\t', header=0)
path_df = filter_transform_labs_to_spoke(path_df, path_len_filter=2)
#
patient_to_spoke = get_patient_to_spoke_df(icd_to_spoke, clarity_to_spoke, path_df)
patient_to_spoke['Patient_ID'] = np.array(patient_to_spoke['Patient_ID'].values, dtype=str)
patient_to_spoke['SPOKE'] = np.array(patient_to_spoke['SPOKE'].values, dtype=str)
patient_to_spoke = patient_to_spoke[patient_to_spoke['Patient_ID'].isin(patient_list)]
#
direct_hit_list = patient_to_spoke.SPOKE.unique()
bin_direct_hit = np.array([node in direct_hit_list for node in node_list])
direct_hit_list = np.array(node_list[bin_direct_hit], dtype=str)
print 'SEP STATS:', direct_hit_list.shape, Counter(node_type_list[bin_direct_hit])
np.save(directory+save_filename+'_direct_hit_list', direct_hit_list, allow_pickle=False)
np.save(directory+save_filename+'_spoke_node_list', node_list, allow_pickle=False)
#
dh_index_df = pd.DataFrame(data=np.array([direct_hit_list, np.arange(len(direct_hit_list))]).T, columns=['SPOKE', 'SEP_Index'])
patient_index_df = pd.DataFrame(data=np.array([patient_list, np.arange(len(patient_list))]).T, columns=['Patient_ID', 'Patient_Index'])
patient_to_spoke = pd.merge(pd.merge(patient_index_df, patient_to_spoke, on='Patient_ID'), dh_index_df, on='SPOKE')
#
patient_to_hetio_connect = np.zeros((len(patient_list), len(direct_hit_list)), dtype=int)
patient_to_hetio_connect[np.array(patient_to_spoke.Patient_Index.values, dtype=int) , np.array(patient_to_spoke.SEP_Index.values, dtype=int)] = 1
print_binary_matrix(patient_to_hetio_connect)
#
