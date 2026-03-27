# Randomization-based GO term enrichment analysis between two EggNogMapper annotations
#
# Usage: python go_enrichment.py [annotation_file] [background_annotation_file] n
#
# Example: python go_enrichment.py UFM1.correlated.go_annotation dump.refined_mcl_clusters.filt_expanded.go_annotation 100000
#

import sys
from random import shuffle
from collections import Counter
import random
import requests


def random_shuffle(d):
    keys = list(d.keys())
    shuffle(keys)
    return dict(list(zip(keys, list(d.values()))))

def reverse_sp_dict(d):
    reversed_dict = {}
    for k, v in d.items():
        reversed_dict.setdefault(v, [])
        reversed_dict[v].append(k.split(".")[0])
    return reversed_dict

def file_import(fname):
    infile = open(fname,'r').readlines()
    goids = []
    for line in infile:
        gos = []
        if len(line.split('\t')[1]) > 1:
            gos = line.split('\t')[1].split(',')
        if len(gos) > 0:
            for ID in gos:
                goids.append(ID.strip())
    return goids

def number_of_annotated_seqs(fname):
    infile = open(fname,'r').readlines()
    m = 0
    for line in infile:
        if len(line.split('\t')[1].strip()) > 1:
            m += 1
    return m

def get_go_term_description(go_id):
    url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_id}"
    headers = {"Accept": "application/json"}

    response = requests.get(url, headers=headers)
    
    if response.status_code == 200:
        data = response.json()
        try:
            return data['results'][0]['name']
        except (KeyError, IndexError):
            return "Description not found."
    else:
        return f"Error: {response.status_code}"

# import the data
exp = file_import(sys.argv[1])

go2count_obs = Counter(exp)

# define number of replicates
rounds = int(sys.argv[3])

annotated_seqs = number_of_annotated_seqs(sys.argv[1])

bkgd = open(sys.argv[2],'r').readlines()

# load in the background file and import annotated sequences
seq2go_bkgd = {}
for line in bkgd:
    if len(line.split('\t')[1].strip()) > 1:
        seq2go_bkgd[line.split('\t')[0]] = ([i.strip() for i in line.strip().split('\t')[1].split(',')])

# make a record dictionary for the permutation test
go2count_check = {}
for key in list(go2count_obs.keys()):
    go2count_check[key] = 0

for n in range(rounds):
    # Just the progress bar, if the dataset is big
    sys.stdout.write('\r')
    i = int(((n+1)/float(rounds))*100)
    sys.stdout.write("[%-100s] %d%%" % ('='*i, i))
    sys.stdout.flush()
    ##############################################
    # randomly select same number of annotated seqs as in test set
    random_gos = []
    random_gos = [go for k in random.sample(list(seq2go_bkgd.keys()),annotated_seqs) for go in seq2go_bkgd[k]]
    go2count_random = Counter(random_gos)
    for key in list(go2count_obs.keys()):
        if key in list(go2count_random.keys()):
            if abs(go2count_random[key]) >= abs(go2count_obs[key]):
                go2count_check[key] += 1
    
print()
print()

go2pval = {}
for key in list(go2count_check.keys()):
    go2pval[key]= go2count_check[key]/float(rounds)

m = 0
n = 0
for go in sorted(list(go2pval.keys()), key=lambda l: go2pval[l], reverse=True):
    if go2pval[go] <= 0.05:
        print(go + '\t' + str(go2pval[go]))
        n += 1
        m += 1
    else:
        n += 1

print(str(m) + ' enriched GO ids found out of ' + str(n) + ' (' + str(float(m)/n) + ')')

bkgd_seqs = number_of_annotated_seqs(sys.argv[2])
bkgd = file_import(sys.argv[2])
go2count_bkgd = Counter(bkgd)

# get descriptions
go_to_desc = {}
for go in go2pval:
    try:
        description = get_go_term_description(go)
        go_to_desc[go] = description
    except:
        go_to_desc[go] = 'NA'

out = open(sys.argv[1]+'.go_enrichment.tab','w')
out.write('goID\tdesc\tPvalue\tFoldEnrch\ttest_count\tbkgd_count\texpected\n')
for go in list(go2pval.keys()):
    if go2pval[go] <= 1:
        # calculate fold enrichment
        expect = float(annotated_seqs) * (go2count_bkgd[go]/float(bkgd_seqs))
        if expect != 0:
            fe = go2count_obs[go]/expect
        else:
            fe = 'INF'
        out.write(go + '\t' + go_to_desc[go] + '\t' + str(go2pval[go]) + '\t' + str(fe) + '\t' + str(go2count_obs[go]) + '\t' + str(go2count_bkgd[go])+'\t'+str(expect)+'\n')
out.close()
    

            
    


    
