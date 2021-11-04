import sys
import re
import string
from collections import OrderedDict
from functools import lru_cache

import atom3d.datasets as da
import pandas as pd
import pprint

path = '/rigel/free/users/gt2453/atom3d_data/'

########### PSR Dataset #########
def has_number(str):
    return bool(re.search(r'\d', str))

def casp_dict():
    casp_dataset = pd.read_csv('/rigel/free/users/gt2453/casp_data/casp13_targetlist.csv', sep=';')
    description = casp_dataset['Description'].str.split(' ')
    target = casp_dataset['Target'].str.upper()
    new_dataset = pd.concat([target, description], axis=1)
    
    dict = OrderedDict()
    for line in new_dataset.itertuples():
        for w in line[2]:
            if len(w) == 4 and has_number(w) and w.islower():
                cur_target = line[1]
                dict[cur_target] = w
    
    # manually add a few PDB IDs from CASP13 that were available on the website
    # but not the file that was downloadable
    dict['H0980'] = '6gnx'
    dict['T0980S1'] = '6gnx'
    dict['T0980S2'] = '6gnx'
    return dict

def generate_psr():
    atom3d_dataset = da.load_dataset(path + 'psr/raw/casp5_to_13/data', 'lmdb')

    counter = 0
    target_list = []
    for data in atom3d_dataset:
        casp_id = data['id']
        target_id = casp_id.split(' ')[0].translate(str.maketrans('', '', string.punctuation))
        target_list.append(target_id)
        
        if counter % 200 == 0:
            print(counter, ':\t', target_id)
            target_list = list(set(target_list)) # remove duplicates
        counter += 1

    target_list = list(set(target_list)) # remove duplicates
    target_list = sorted(target_list) # sort
    target_str = ' \n'.join(target_list)

    with open('/rigel/free/users/gt2453/casp_data/casp_target_ids.txt', 'w') as f: # open and write to file    
        f.write(target_str) # write to file

########## PPI Dataset ########
@lru_cache()
def generate_ppi():
    atom3d_dataset = da.load_dataset(path + 'ppi/raw/DIPS/data', 'lmdb')

    pdb_id_list = [] # store in list
    counter = 0
    print('Len:\t', len(atom3d_dataset))
    for data in atom3d_dataset:
        id_list = data['id'].split('_') # parse id from dataset, refer to atom3D documentation
        id1 = id_list[0].split('.')[0] # parse for pdb id 1
        id2 = id_list[3].split('.')[0] # parse for pdb id 2

        pdb_id_list.append(id1) # add pdb id 1 to list
        
        # check if id1 and id2 are different
        if id1 != id2:
            pdb_id_list.append(id2) # add pdb id 2 to list

        # print counter updates
        if counter % 5000 == 0:
            print(counter)
            pdb_id_list = list(set(pdb_id_list))
        counter += 1
    
    pdb_ids = list(set(pdb_id_list))
    num_pdb_ids = len(pdb_ids)
    print('Num PDB:\t', num_pdb_ids)

    index = 0
    while index < num_pdb_ids:
        try:
            pdb_id_tmp = ' '.join(pdb_ids[index:index+1000])
        except:
            pdb_id_temp = ' '.join(pdb_ids[index:])
        
        num_file = int(index / 1000)
        with open(f'/rigel/free/users/gt2453/pdb_ids/ppi/ppi_ids_{num_file}.txt', 'w') as f: # open and write to file    
            f.write(pdb_id_tmp) # write to file
        index += 1000

def generate_ppi2():
    atom3d_dataset = da.load_dataset(path + 'ppi/raw/DIPS/data', 'lmdb')

    pdb_id_list = [] # store in list
    counter = 0
    print('Len:\t', len(atom3d_dataset))
    for data in atom3d_dataset:
        id_list = data['id'].split('_') # parse id from dataset, refer to atom3D documentation
        id1 = id_list[0].split('.')[0] # parse for pdb id 1
        id2 = id_list[3].split('.')[0] # parse for pdb id 2

        pdb_id_list.append(id1) # add pdb id 1 to list
        
        # check if id1 and id2 are different
        if id1 != id2:
            pdb_id_list.append(id2) # add pdb id 2 to list

        # print counter updates
        if counter % 5000 == 0:
            print(counter)
            pdb_id_list = list(set(pdb_id_list))
        counter += 1
    
    pdb_ids = list(set(pdb_id_list))
    num_pdb_ids = len(pdb_ids)
    print('Num PDB:\t', num_pdb_ids)

    # write to file
    pdb_id_tmp = '\n'.join(pdb_ids)
    with open(f'/rigel/free/users/gt2453/pdb_ids/ppi/ppi_ids.txt', 'w') as f: # open file    
            f.write(pdb_id_tmp) # write to file

#generate_ppi2()
######### MSP Dataset #######
def generate_msp():
    dataset = da.load_dataset(path + 'msp/raw/MSP/data', 'lmdb')

    # get all PDB IDs from dataset
    id_list = [] # store pdb IDs in list
    counter = 0
    print('Len of Dataset:\t', len(dataset))
    for data in dataset:
        counter += 1

        id = data['id'].split('_')[0]
        id_list.append(id)
        
        if counter % 100 == 0:
            print(counter)
    
    id_list = list(set(id_list)) # ensure no duplicates
    
    # write PDB IDs to file. Max 1000 IDs in one file because of Protein Data Bank downloader tool restrictions 
    index = 0
    while index < len(id_list):
        try:
            tmp_ids = ' '.join(id_list[index:index+1000])
        except:
            tmp_ids = ' '.join(id_list[index:])
        
        num_file = int(index / 1000)
        with open(f'/rigel/free/users/gt2453/pdb_ids/msp/msp_ids_{num_file}.txt', 'w') as f: # open and write to file    
            f.write(tmp_ids) # write to file
        index += 1000
    
def generate_msp2():
    dataset = da.load_dataset(path + 'msp/raw/MSP/data', 'lmdb')

    # get all PDB IDs from dataset
    id_list = [] # store pdb IDs in list
    counter = 0
    print('Len of Dataset:\t', len(dataset))
    for data in dataset:
        counter += 1

        id = data['id'].split('_')[0]
        id_list.append(id)
        
        if counter % 100 == 0:
            print(counter)
    
    id_list = list(set(id_list)) # ensure no duplicates
    
    # write PDB IDs to file. Max 1000 IDs in one file because of Protein Data Bank downloader tool restrictions 
    index = 0
    while index < len(id_list):
        try:
            tmp_ids = '\n'.join(id_list[index:index+1000])
        except:
            tmp_ids = '\n'.join(id_list[index:])
        
        num_file = int(index / 1000)
        with open(f'/rigel/free/users/gt2453/pdb_ids/msp/msp_ids.txt', 'w') as f: # open and write to file    
            f.write(tmp_ids) # write to file
        index += 1000

######### RES Dataset #########

def generate_res():
    atom3d_dataset = da.load_dataset(path + 'res/raw/RES/data', 'lmdb')

    pdb_id_list = [] # store in list
    counter = 0
    print('Len:\t', len(atom3d_dataset))
    for data in atom3d_dataset:
        if data == None:
            continue

        id = data['id']# parse id from dataset, refer to atom3D documentation
        pdb_id_list.append(id) # add pdb id to list
        
        # print counter updates
        if counter % 5000 == 0:
            print(counter)
            pdb_id_list = list(set(pdb_id_list))
        counter += 1
    
    pdb_ids = list(set(pdb_id_list))
    num_pdb_ids = len(pdb_ids)
    print('Num PDB:\t', num_pdb_ids)

    # write to file
    pdb_id_tmp = '\n'.join(pdb_ids)
    with open(f'/rigel/free/users/gt2453/pdb_ids/res/res_ids.txt', 'w') as f: # open file    
            f.write(pdb_id_tmp) # write to file

generate_res()