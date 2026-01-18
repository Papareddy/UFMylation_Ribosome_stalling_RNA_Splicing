
import sys
import math
from collections import defaultdict
import os

# Paths to data files (relative to this script)
# We assume the maxentpy repo was cloned to src/maxentpy_repo
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "maxentpy_repo", "maxentpy", "data")

bgd_5 = {'A': 0.27, 'C': 0.23, 'G': 0.23, 'T': 0.27}
cons1_5 = {'A': 0.004, 'C': 0.0032, 'G': 0.9896, 'T': 0.0032}
cons2_5 = {'A': 0.0034, 'C': 0.0039, 'G': 0.0042, 'T': 0.9884}

bgd_3 = {'A': 0.27, 'C': 0.23, 'G': 0.23, 'T': 0.27}
cons1_3 = {'A': 0.9903, 'C': 0.0032, 'G': 0.0034, 'T': 0.0030}
cons2_3 = {'A': 0.0027, 'C': 0.0037, 'G': 0.9905, 'T': 0.0030}

try:
    from string import maketrans
except ImportError:
    maketrans = str.maketrans

def hashseq(fa):
    table = maketrans('ACGT', '0123')
    seq = fa.translate(table)
    return sum(int(j) * 4**(len(seq) - i - 1) for i, j in enumerate(seq))

def load_matrix5():
    matrix_f = os.path.join(DATA_DIR, 'score5_matrix.txt')
    matrix = {}
    if not os.path.exists(matrix_f):
        raise FileNotFoundError(f"MaxEnt5 matrix not found at {matrix_f}")
    with open(matrix_f, 'r') as f:
        for line in f:
            entry = line.split()
            matrix[entry[0]] = float(entry[1])
    return matrix

def load_matrix3():
    matrix_f = os.path.join(DATA_DIR, 'score3_matrix.txt')
    matrix = defaultdict(dict)
    if not os.path.exists(matrix_f):
        raise FileNotFoundError(f"MaxEnt3 matrix not found at {matrix_f}")
    with open(matrix_f, 'r') as f:
        for line in f:
            n, m, s = line.split()
            matrix[int(n)][int(m)] = float(s)
    return matrix

# Global caches
_matrix5 = None
_matrix3 = None

def score5(fa):
    '''Calculate 5' splice site strength (9-mer)'''
    global _matrix5
    fa = fa.upper()
    if len(fa) != 9: return 0.0 # Robustness: return 0 if wrong length instead of crash
    
    if not _matrix5: _matrix5 = load_matrix5()
    
    # check bases valid
    if any(c not in 'ACGT' for c in fa): return 0.0

    key = fa[3:5]
    score = cons1_5[key[0]] * cons2_5[key[1]] / (bgd_5[key[0]] * bgd_5[key[1]])
    rest = fa[:3] + fa[5:]
    rest_score = _matrix5[rest]
    return math.log(score * rest_score, 2)

def score3(fa):
    '''Calculate 3' splice site strength (23-mer)'''
    global _matrix3
    fa = fa.upper()
    if len(fa) != 23: return 0.0 # Robustness
    
    if not _matrix3: _matrix3 = load_matrix3()

    if any(c not in 'ACGT' for c in fa): return 0.0

    key = fa[18:20]
    score = cons1_3[key[0]] * cons2_3[key[1]] / (bgd_3[key[0]] * bgd_3[key[1]])
    
    rest = fa[:18] + fa[20:]
    rest_score = 1
    rest_score *= _matrix3[0][hashseq(rest[:7])]
    rest_score *= _matrix3[1][hashseq(rest[7:14])]
    rest_score *= _matrix3[2][hashseq(rest[14:])]
    rest_score *= _matrix3[3][hashseq(rest[4:11])]
    rest_score *= _matrix3[4][hashseq(rest[11:18])]
    rest_score /= _matrix3[5][hashseq(rest[4:7])]
    rest_score /= _matrix3[6][hashseq(rest[7:11])]
    rest_score /= _matrix3[7][hashseq(rest[11:14])]
    rest_score /= _matrix3[8][hashseq(rest[14:18])]
    return math.log(score * rest_score, 2)

if __name__ == "__main__":
    # Test
    print(f"Testing MaxEnt Scorer...")
    try:
        s5 = score5("CAGGTAAGT") # Should be ~10.86
        s3 = score3("TTCCAAACGAACTTTTGTAGGGA") # Should be ~2.89
        print(f"5'SS (CAGGTAAGT): {s5:.2f}")
        print(f"3'SS (TTCC...): {s3:.2f}")
    except Exception as e:
        print(f"Error: {e}")
