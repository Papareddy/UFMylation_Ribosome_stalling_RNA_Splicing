# Randomization-based cluster enrichment analysis between 2 phenotypic groups
#
# Usage: python profiling.py [orthogroup_file] [clade_file] [n] [min taxa] [binary (True/False)] [threads]
# Example: python profiling.py dump.refined_mcl_clusters.filt_expanded UFM1.clade_file.txt 100000 False 32
#
# Notes:
# 
# 1. The orthogroup file should be a tab delimited file where each line is an orthogroup and each field is a sequence. Eg. dump file from mcl. 
# 2. Sequences should be formatted as Genus_species.seqID (or at least taxonomy.ID) - you can use different delimiters if you alter the script
# 3. The clade file should comprise two columns. 
#    Column one is the phenotype class (eg. present and absent). 
#    Column two is a tab delimited list of species.
#    Each row is a different clade.

import sys
from random import shuffle
from collections import Counter
from statistics import mean
import pandas as pd
from scipy.stats import pointbiserialr
import numpy as np
import multiprocessing

pd.set_option("display.max_columns", None)

num_threads = int(sys.argv[6])
binary = bool(sys.argv[5])

def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (sum(s[n//2-1:n//2+1])/2.0, s[n//2])[n % 2] if n else None

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

def get_metric_dict(cl2seq,x,binary):
    cl2metric = {}
    cl2m1 = {}
    cl2m2 = {}
    pbc_d = {}
    if x == 1:
        out = open(sys.argv[2].split('.')[0]+'.metric_per_cluster.tab','w')
        out.write('OG\t'+classes[0]+'\t'+classes[1]+'\tdiff\tpbc\n')
    for cl in cl2seq:
        counted = Counter(cl2seq[cl])
        if binary == True:
            for t in counted:
                if counted[t] > 0:
                    counted[t] = 1
        value_d = {}
        for group in list(group_d.keys()):
            value_d['v'+str(group.strip('group'))] = median([counted[sp] for sp in group_d[group]])
        class1_values = []
        class2_values = []
        for v in list(value_d.keys()):
            if 'group' + v.strip('v') in class1:
                class1_values.append(value_d[v])
            else:
                class2_values.append(value_d[v])
        m1 = mean(class1_values)
        m2 = mean(class2_values)
        cl2metric[cl] = m1 - m2
        class1_consistency = 1-(class1_values.count(0)/len(class1_values))
        class2_consistency = 1-(class2_values.count(0)/len(class2_values))
        cl2m1[cl] = [median(class1_values),mean(class1_values),class1_consistency]
        cl2m2[cl] = [median(class2_values),mean(class2_values),class2_consistency]
        # calculate correlation
        y = 'corr'
        if y == 'corr':
            # calculate PBC (point biserial correlation)
            corr_d = {'phenotype':[], 'og':[]}
            for i in class1:
                corr_d['phenotype'].append(1)
                corr_d['og'].append(value_d['v'+i.strip('group')])
            for i in class2:
                corr_d['phenotype'].append(0)
                corr_d['og'].append(value_d['v'+i.strip('group')]) 
            df = pd.DataFrame(corr_d)
            if len(set(corr_d['og'])) > 1:
                r_pb, p_value = pointbiserialr(corr_d['phenotype'],corr_d['og'])
            else:
                r_pb = 0
            pbc_d[cl] = round(r_pb,4)
        # set metric as correlation
        cl2metric[cl] = pbc_d[cl]
        pbc = round(r_pb,3)
        # output
        if x == 1:
            out.write(cl + '\t' + str(float(m1)) + '\t' + str(float(m2)) + '\t' + str(m1-m2) + '\t' + str(pbc) +'\n')
    if x == 1:    
        out.close()
    return cl2metric, cl2m1, cl2m2, pbc_d

# Randomization rounds
rounds = int(sys.argv[3])

# Group definition
group_file = open(sys.argv[2],'r').readlines()
classes = list(set([c.split('\t')[0].strip() for c in group_file]))
classes.sort()
print(classes)

group_d = {}
class1 = []
class2 = []

n = 1
for g in group_file:
    group_d['group'+str(n)] = [t.strip() for t in g.split('\t')[1:]]
    if g.split('\t')[0] == classes[0]:
        class1.append('group'+str(n))
    else:
        class2.append('group'+str(n))
    n += 1

# Format the orthogroup file
ortho = open(sys.argv[1],'r').readlines()
out = open(sys.argv[2].split('.')[0]+'.orthogroup_assignment.tab','w')

# determine maximum size of orthogroups - very large ones become ambiguous
sizes = [og.count(',')+1 for og in ortho]
max_size = 200000

for orthogroup in ortho:
    og = orthogroup.split('\t')[0]
    sequences = [s.strip() for s in orthogroup.split('\t')[1].split(',') if s != '']
    # only consider orthogroups with at least n taxa
    if len(set([t.split('.')[0] for t in sequences])) >= int(sys.argv[4]):
        # exclude giant orthogroups
        if len(sequences) <= max_size:
            for seq in sequences:
                out.write(seq.strip()+'\t'+og+'\n')
out.close()

# The tab delimited sequence to cluster mapping as first argument
seq2cl = dict([line.split("\t")[0], line.strip().split("\t")[1]] for line in open(sys.argv[2].split('.')[0]+'.orthogroup_assignment.tab'))
cl2seq_observed = reverse_sp_dict(seq2cl)
# Choice of metric to use in function 
cl2metric_observed, cl2m1, cl2m2, pbc_observed = get_metric_dict(cl2seq_observed,1,binary)

cl2counts = {}

print('\nComparing ' + classes[0] + ' to ' + classes[1] + '\n')

n_d = {}
def randomizer(n):
    seq2cl_random = random_shuffle(seq2cl)
    cl2seq_random = reverse_sp_dict(seq2cl_random)
    cl2metric_random, r1, r2, pbc_random = get_metric_dict(cl2seq_random,2,binary)
    result = []
    for cl in cl2metric_observed:
        cl2counts.setdefault(cl, 0)
        # Two-tail test. over or under-representation
        if abs(cl2metric_random[cl]) >= abs(cl2metric_observed[cl]):
            result.append(1)
        else:
            result.append(0)
    return result

with multiprocessing.Pool(processes=num_threads) as pool:
    results = pool.map(randomizer, list(range(0,rounds)))

n = 0
for cl in cl2metric_observed:
    v = np.sum([r[n] for r in results])
    cl2counts[cl] = v
    n += 1

cl2pval = {cl: cl2counts[cl] / float(rounds) for cl in cl2counts.keys()}

e = 0
d = 0
t = 0

out = open(sys.argv[2].split('.')[0]+'.correlated_clusters.tab','w')
out.write('cluster\ttype\tpvalue\t'+classes[0]+'_median\t'+classes[0]+'_average\t'+classes[0]+'_consistency\t'+classes[1]+'_median\t'+classes[1]+'_average\t'+classes[1]+'_consistency\tpbc\n')
for cl in cl2pval:
    pval = cl2pval[cl]
    if pval <= 0.01:
        if cl2metric_observed[cl] > median([cl2metric_observed[k] for k in list(cl2metric_observed.keys())]): # check what the median difference is. If metric is above = enriched
            enr = classes[0]+'_enriched'
            if cl2m1[cl][2] >= 0.5:
                e += 1
        else:
            enr = classes[0]+'_depleted'
            if cl2m2[cl][2] >= 0.5:
                d += 1
        if ('_enriched' in enr) and (cl2m1[cl][2] >= 0.5):
            out.write(cl + '\t' + enr + '\t' + str(pval) + '\t' + str(cl2m1[cl][0]) + '\t' + str(cl2m1[cl][1]) + '\t' + str(cl2m1[cl][2]) + '\t' + str(cl2m2[cl][0]) + '\t' + str(cl2m2[cl][1]) + '\t' + str(cl2m2[cl][2]) + '\t' + str(round(pbc_observed[cl],4)) + '\n')
        elif ('_depleted' in enr) and (cl2m2[cl][2] >= 0.5):
            out.write(cl + '\t' + enr + '\t' + str(pval) + '\t' + str(cl2m1[cl][0]) + '\t' + str(cl2m1[cl][1]) + '\t' + str(cl2m1[cl][2]) + '\t' + str(cl2m2[cl][0]) + '\t' + str(cl2m2[cl][1]) + '\t' + str(cl2m2[cl][2]) + '\t' + str(round(pbc_observed[cl],4)) + '\n')
    else:
        enr = 'not significant'
        out.write(cl + '\t' + enr + '\t' + str(pval) + '\t' + str(cl2m1[cl][0]) + '\t' + str(cl2m1[cl][1]) + '\t' + str(cl2m1[cl][2]) + '\t' + str(cl2m2[cl][0]) + '\t' + str(cl2m2[cl][1]) + '\t' + str(cl2m2[cl][2]) + '\t' + str(round(pbc_observed[cl],4)) + '\n')
    t += 1

print('\n' + str(round(100*(float(e)/t),2)) + '% enriched in ' + classes[0])
print(str(round(100*(float(d)/t),2)) + '% depleted in ' + classes[0])
print('from a total of ' + str(t) + ' orthogroups\n')

out.close()

    

