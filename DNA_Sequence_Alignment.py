"""
Author: James Roberts
Last Updated: 10/20/2022

This program performs:
global alignment
Semi-global alignment (prefers gaps at ends)
Local alignment (No gaps)

using the Needleman-Wunsch algorithm.

You can change the sequences I provide in code (seq1 or seq2)
or you can read a FASTA sequence in to use with code. 
"""
# Constants
base_idx = {'A': 0, 'G': 1, 'C': 2, 'T': 3}
PTR_NONE, PTR_GAP1, PTR_GAP2, PTR_BASE = 0, 1, 2, 3
ALIGN_GLOBAL, ALIGN_SEMIGLOBAL, ALIGN_LOCAL = 0, 1, 2
S = [
    # A   G   C   T
    [ 3, -1, -2, -2],  # A
    [-1,  3, -2, -2],  # G
    [-2, -2,  3, -1],  # C
    [-2, -2, -1,  3]   # T
]
gap_penalty = 4
seq1 = 'TACGCAG'
seq2 = 'ACGTG'

# Functions
def seqalignDP(seq1, seq2, subst_matrix, gap_penalty, align_type):
    """
    Return the score of the optimal Needdleman-Wunsch alignment for seq1
    and seq2 along with traceback and Score Matrix.
    Note: gap_penalty should be positive (it is subtracted)
    """
    F = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    TB = [[PTR_NONE for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]


    # Initializations based on Align type
    if align_type == ALIGN_GLOBAL or align_type == ALIGN_SEMIGLOBAL:
        count = 1
        for i in range(len(seq1)):
            TB[i + 1][0] = PTR_GAP2
            count += 1

        count = 1
        for j in range(len(seq2)):
            TB[0][j + 1] = PTR_GAP1
            count += 1

    # Global
    if align_type == ALIGN_GLOBAL:
        count = 1
        for i in range(len(seq1)):
            F[i + 1][0] = -1 * gap_penalty * count
            count += 1

        count = 1
        for j in range(len(seq2)):
            F[0][j + 1] = -1 * gap_penalty * count
            count += 1

    # Semi Local
    if align_type == ALIGN_SEMIGLOBAL:
        if len(seq2) >= len(seq1):
            count = 1
            for i in range(len(seq1)):
                F[i + 1][0] = -1 * gap_penalty * count
                count += 1

        else:
            count = 1
            for j in range(len(seq2)):
                F[0][j + 1] = -1 * gap_penalty * count
                count += 1

    # Local (row 1 and col 1 are all zeros)
    if align_type == ALIGN_LOCAL:
        pass

    # Filling in the Score Matrix (F)
    S = subst_matrix
    if align_type == ALIGN_SEMIGLOBAL or align_type == ALIGN_GLOBAL:
        for i in range(len(seq1)):
            for j in range(len(seq2)):

                # Check match or mismatch
                if seq1[i] == seq2[j]:
                    match = S[0][0] #Score of a match
                elif seq1[i] == 'A' and seq2[j] == 'G':
                    match = S[1][0]
                elif seq1[i] == 'A' and (seq2[j] == 'T' or seq2[j] == 'C'):
                    match = S[2][0]
                elif seq1[i] == 'G' and seq2[j] == 'A':
                    match = S[1][0]
                elif seq1[i] == 'G' and (seq2[j] == 'T' or seq2[j] == 'C'):
                    match = S[2][0]
                elif seq1[i] == 'C' and seq2[j] == 'T':
                    match = S[1][0]
                elif seq1[i] == 'C' and (seq2[j] == 'A' or seq2[j] == 'G'):
                    match = S[2][0]
                elif seq1[i] == 'T' and seq2[j] == 'C':
                    match = S[1][0]
                elif seq1[i] == 'T' and (seq2[j] == 'A' or seq2[j] == 'G'):
                    match = S[2][0]

                m = PTR_BASE
                gap_y = PTR_GAP1
                gap_x = PTR_GAP2

                score_m = F[i][j] + match
                score_gap_y = F[i + 1][j] - gap_penalty
                score_gap_x = F[i][j + 1] - gap_penalty

                score = max([
                    score_m,
                    score_gap_y,
                    score_gap_x])

                if score_m == score:
                    comefrom = m
                if score_gap_x == score:
                    comefrom = gap_x
                if score_gap_y == score:
                    comefrom = gap_y

                F[i+1][j+1] = score
                TB[i+1][j+1] = comefrom

    if align_type == ALIGN_LOCAL:
        for i in range(len(seq1)):
            for j in range(len(seq2)):

                # Check match or mismatch
                if seq1[i] == seq2[j]:
                    match = S[0][0]  # Score of a match
                elif seq1[i] == 'A' and seq2[j] == 'G':
                    match = S[1][0]
                elif seq1[i] == 'A' and (seq2[j] == 'T' or seq2[j] == 'C'):
                    match = S[2][0]
                elif seq1[i] == 'G' and seq2[j] == 'A':
                    match = S[1][0]
                elif seq1[i] == 'G' and (seq2[j] == 'T' or seq2[j] == 'C'):
                    match = S[2][0]
                elif seq1[i] == 'C' and seq2[j] == 'T':
                    match = S[1][0]
                elif seq1[i] == 'C' and (seq2[j] == 'A' or seq2[j] == 'G'):
                    match = S[2][0]
                elif seq1[i] == 'T' and seq2[j] == 'C':
                    match = S[1][0]
                elif seq1[i] == 'T' and (seq2[j] == 'A' or seq2[j] == 'G'):
                    match = S[2][0]

                m = PTR_BASE
                gap_y = PTR_GAP1
                gap_x = PTR_GAP2

                score_m = F[i][j] + match
                score_gap_y = F[i + 1][j] - gap_penalty
                score_gap_x = F[i][j + 1] - gap_penalty

                score = max([
                    score_m,
                    score_gap_y,
                    score_gap_x])

                if score < 0:
                    score = 0

                if score_m == score:
                    comefrom = m
                elif score_gap_x == score:
                    comefrom = gap_x
                elif score_gap_y == score:
                    comefrom = gap_y
                else:
                    comefrom = PTR_NONE

                F[i+1][j+1] = score
                TB[i+1][j+1] = comefrom

    # Max score for each align_type
    if align_type == ALIGN_GLOBAL:
        F_max = F[len(seq1)][len(seq2)]
        index = [len(seq2), len(seq1)]
    elif align_type == ALIGN_SEMIGLOBAL:

        if len(seq2) >= len(seq1):
            end = []
            for i in range(len(seq2) + 1):
                end.append(F[-1][i])
            a = max(end)
            index = [end.index(a), len(seq1)]

            F_max = a
        if len(seq1) >= len(seq2):
            end = []
            for i in range(len(seq1)+1):
                end.append(F[i][-1])
            a = max(end)
            index = [len(seq2), end.index(a)]
            for i in range(len(seq1) - index[1]):
                TB[len(seq1)][-1 * (i+1)] = PTR_GAP2

            F_max = a
    elif align_type == ALIGN_LOCAL:
        index = []
        value = 0
        for i in range(len(seq1)+1):
            for j in range(len(seq2)+1):
                if F[i][j] > value:
                    index.append([i,j])
                    value = F[i][j]
        F_max = value
        index = index[-1]
    else:
        assert False

    return  F, F_max ,TB, index

def traceback(seq1, seq2, TB, align_type, index=0):
    """
    This function returns optimal alignment and traceback.
    """
    s1 = ""
    s2 = ""

    i = len(seq1)
    j = len(seq2)

    if align_type == ALIGN_GLOBAL:
        while TB[i][j] != PTR_NONE:
            if TB[i][j] == PTR_BASE:
                s1 = seq1[i - 1] + s1
                s2 = seq2[j - 1] + s2
                TB[i][j] = 'p'
                i = i - 1
                j = j - 1
            elif TB[i][j] == PTR_GAP1:
                s1 = '-' + s1
                s2 = seq2[j - 1] + s2
                TB[i][j] = 'p'
                j = j - 1
            elif TB[i][j] == PTR_GAP2:
                s1 = seq1[i - 1] + s1
                s2 = '-' + s2
                TB[i][j] = 'p'
                i = i - 1
            else:
                assert False

    if align_type == ALIGN_SEMIGLOBAL:
        while TB[i][j] != PTR_NONE:
            if TB[i][j] == PTR_BASE:
                s1 = seq1[i - 1] + s1
                s2 = seq2[j - 1] + s2
                TB[i][j] = 'p'
                i = i - 1
                j = j - 1
            elif TB[i][j] == PTR_GAP1:
                s1 = '-' + s1
                s2 = seq2[j - 1] + s2
                TB[i][j] = 'p'
                j = j - 1
            elif TB[i][j] == PTR_GAP2:
                s1 = seq1[i - 1] + s1
                s2 = '-' + s2
                TB[i][j] = 'p'
                i = i - 1
            else:
                assert False

    if align_type == ALIGN_LOCAL:
        i = index[0]
        j = index[1]
        while TB[i][j] != PTR_NONE:
            if TB[i][j] == PTR_BASE:
                s1 = seq1[i - 1] + s1
                s2 = seq2[j - 1] + s2
                TB[i][j] = 'p'
                i = i - 1
                j = j - 1
            elif TB[i][j] == PTR_GAP1:
                s1 = '-' + s1
                s2 = seq2[j - 1] + s2
                TB[i][j] = 'p'
                j = j - 1
            elif TB[i][j] == PTR_GAP2:
                s1 = seq1[i - 1] + s1
                s2 = '-' + s2
                TB[i][j] = 'p'
                i = i - 1
            else:
                assert False

    if align_type == ALIGN_GLOBAL or align_type == ALIGN_SEMIGLOBAL:
        TB[0][0] = 'p'

    # Distance Metric
    distance = len(s2)
    for i in range(len(s2)):
        if s2[i] == s1[i]:
            distance += -1

    return s1, s2, TB, distance

def readSeq(filename):
    """Reads in a FASTA sequence. Assumes one sequence in the file"""
    seq = []

    with open(filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.rstrip().upper())

    return "".join(seq)


print('---------------------------------------------------------------------')
F_Score_TB_index = seqalignDP(seq1=seq1,
                              seq2=seq2,
                              subst_matrix=S,
                              gap_penalty=gap_penalty,
                              align_type=ALIGN_GLOBAL)
sequences_distance = traceback(seq1 = seq1,
                               seq2 = seq2,
                               TB = F_Score_TB_index[2],
                               align_type=ALIGN_GLOBAL)
print('An optimal global alignment for sequences_distance ACGTG and TACGCAG:')
print(sequences_distance[0])
print(sequences_distance[1])

print('\nScore of alignment: ' + str(F_Score_TB_index[1]) + '\n')

print('Score Matrix F: ')
for b in range(len(sequences_distance[1])+1):
    print(F_Score_TB_index[0][b])

print()

print('Optimal Path are the p: ')
for b in range(len(sequences_distance[1])+1):
    print(F_Score_TB_index[2][b])

print('------------------------------------------------------------')
F_Score_TB_index = seqalignDP(seq1=seq1,
                              seq2=seq2,
                              subst_matrix=S,
                              gap_penalty=gap_penalty,
                              align_type=ALIGN_SEMIGLOBAL)
sequences_distance = traceback(seq1 = seq1,
                               seq2 = seq2,
                               TB = F_Score_TB_index[2],
                               align_type=ALIGN_SEMIGLOBAL)
print('An optimal semi-global alignment for sequences_distance ACGTG and TACGCAG:')
print(sequences_distance[0])
print(sequences_distance[1])

print('\nScore of alignment: ' + str(F_Score_TB_index[1]) + '\n')

print('Score Matrix F: ')
for b in range(len(sequences_distance[1])+1):
    print(F_Score_TB_index[0][b])

print()

print('Optimal Path are the p: ')
for b in range(len(sequences_distance[1])+1):
    print(F_Score_TB_index[2][b])

print('------------------------------------------------------------')
F_Score_TB_index = seqalignDP(seq1=seq1,
                              seq2=seq2,
                              subst_matrix=S,
                              gap_penalty=gap_penalty,
                              align_type=ALIGN_LOCAL)
sequences_distance = traceback(seq1 = seq1,
                               seq2 = seq2,
                               TB = F_Score_TB_index[2],
                               align_type=ALIGN_LOCAL,
                               index=F_Score_TB_index[3])
print('An optimal Local alignment for sequences_distance ACGTG and TACGCAG:')
print(sequences_distance[0])
print(sequences_distance[1])

print('\nScore of alignment: ' + str(F_Score_TB_index[1]) + '\n')

print('Score Matrix F: ')
for b in range(len(seq1)+1):
    print(F_Score_TB_index[0][b])

print()

print('Optimal Path are the p: ')
for b in range(len(seq1)+1):
    print(F_Score_TB_index[2][b])
