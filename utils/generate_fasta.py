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

@lru_cache()
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

atom3d_dataset = da.load_dataset(path + 'psr/raw/casp5_to_13/data', 'lmdb')
#dict= casp_dict()
#pprint.pprint(dict)

counter = 0
id_list = []
target_list = []
for data in atom3d_dataset:
    casp_id = data['id']
    target_id = casp_id.split(' ')[0].translate(str.maketrans('', '', string.punctuation))

    target_list.append(target_id)
    if counter % 200 == 0:
        print(counter, ':\t', target_id)
    counter += 1

id_list = list(set(id_list)) # remove duplicates
target_list = list(set(target_list))
print(len(target_list))

target_str = ' '.join(target_list)
with open('/rigel/free/users/gt2453/casp_data/casp_target_ids.txt', 'w') as f: # open and write to file    
    f.write(target_str) # write to file


    