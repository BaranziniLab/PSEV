# PSEV

The data provided here is not the same as that used in the paper. 
Patient IDs have been scrambled.
This is only a fraction of the patient data used.

###### Clone repository
git clone https://github.com/baranzini-lab/PSEV.git
###### Go to constants directory
cd PSEV/constants/
###### Unzip lab path file
tar -xvzf cui_to_loinc_path.tsv.tar.gz
###### Go to build directory
cd ../build/
###### Unzip EHR files
tar -xvzf EHR_DIAGNOSIS.tsv.tar.gz 

tar -xvzf EHR_LABS.tsv.tar.gz

tar -xvzf EHR_MEDICATION.tsv.tar.gz 
###### Go to SPOKE directory
cd ../spoke_v_1/
###### Unzip SPOKE edge file
tar -xvzf neo4j_edges.tsv.tar.gz 
###### Go to main PSEV directory
cd ..
###### Create matrix that contains patient connections to SPOKE nodes
python connect_all_patients.py 
###### stdout (using default settings):
>\# of LOINC mapped to SPOKE: 557
>
>Mapped SPOKE breakdown from LABS: Counter({'SideEffect': 226, 'Compound': 136, 'Gene': 49, 'Symptom': 12})
>
>SEP STATS: (3183,) Counter({'SideEffect': 1851, 'Compound': 922, 'Symptom': 257, 'Disease': 137, 'Gene': 16})
###### Create PSEV matrix for specific SEPs (by node type: -t Gene) 
python make_psevs_by_node_type.py 
>**NOTE**
>
>**Default creates PSEVs for Gene SEPs.**
>_To change default settings (including number of cores to use) please see python make_psevs_by_node_type.py --help_
>
>**Beware this program will use a lot of memory and take a lot of time. Run on cluster if possible.**
###### stdout (using default settings):
>make_psevs_by_node_type.py:83: RuntimeWarning: invalid value encountered in divide
>  add_val = np.nan_to_num(direct_hit.astype(float)/np.sum(direct_hit))
>make_psevs_by_node_type.py:111: RuntimeWarning: invalid value encountered in divide
>  connectivity_matrix = np.transpose(np.nan_to_num(connectivity_matrix/np.sum(connectivity_matrix, axis=0)))
