#Find a k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern
import collections

def create_k_pos (k):
    nt = ["A","T","C","G"]
    possible_k = []
    if k == 1:
       return nt
    else:
        for nuc in nt:
           prev_comb = create_k_pos(k-1)
           i = 0
           while i < len(prev_comb):
                comb =  nuc + prev_comb[i]
                possible_k.append(comb)
                i+=1
    return possible_k


def det_seq(sequence,k):
    nt = ["A","T","G","C"]
    poss_arr = []
    for i in range(0,len(sequence)):
        possible_kmer = sequence[i:i+k]
        if len(possible_kmer) == k:
            possible_kmer_list = list(possible_kmer)
            temp_poss_arr = []
            for n in range(0,len(possible_kmer_list)):
                if n == 0:
                        store = possible_kmer_list[n]
                else:
                        possible_kmer_list[n-1] = store
                        store = possible_kmer_list[n]
                i = 0
                while i <= 3:
                        possible_kmer_list[n] = nt[i]
                        possible_kmer_restored = "".join(possible_kmer_list)
                        temp_poss_arr.append(possible_kmer_restored)
                        i+=1  
            rm_dup = []
            for kseq in temp_poss_arr:
                if not kseq in rm_dup:
                        rm_dup.append(kseq)
            poss_arr.extend(rm_dup)
    return poss_arr


def all_mismatches (sequence, k, mismatch):
    if mismatch == 0:
        return sequence
    elif mismatch == 1:   
        return (det_seq(sequence,k))
    else:
        overall_seq = []
        temp_dic = {}
        for i in range(0,len(sequence)):
            possible_kmer = sequence[i:i+k]
            if len(possible_kmer) == k:
                recur_mis = all_mismatches(possible_kmer, k, mismatch-1)
                for recur_seq in recur_mis:
                    temp_arr = det_seq(recur_seq,k)
                    for ele in temp_arr:
                        temp_dic[ele] = 0
            overall_seq.extend(temp_dic.keys())
            temp_dic = {}
        return overall_seq


def minimize_score (sequence, k):
    poss = create_k_pos(k)
    score_dic = {}
    for seq in poss:
        for indseq in sequence:
            i = 0
            flag = False
            while not flag:
                if i == k:
                   flag = True
                   score_dic[seq] += k
                   break
                if i == 0:
                    check_poss = [seq]
                else:
                    check_poss = all_mismatches(seq,k,i)
                for check in check_poss:
                    if check in indseq and len(check) == (k):   
                        if not seq in score_dic:
                            score_dic[seq] = 0
                        score_dic[seq] += i
                        flag = True
                        break
                i+=1
    for key in score_dic:
        if score_dic[key] == min(score_dic.values()):
           inter_sequence = key
        # print(score_dic)
    return inter_sequence

input = [
"TCTATAGGTTGTGATCTGGACGTTATAAGCAAGCCTTTGACC",
"AAGCCTGGTCTCGTGCGGCCGACCTTGTCCCCCGCGTACTAT",
"TCCCTGATGTAATTGGGTTATATGGACATGGGGATTAAGGCT",
"GGTAAAATGTACTGGTGCTCTCTGAAGGCTCCACAATTGCCC",
"AAGTCTAGCCACGATCTCAGTCAGGCGGAGCTTTCTTCAACC",
"GCTCAAAATTGACTGCTATACGCGAAGGCTGACAAAAAGCGC",
"TCGCTAAAGACTAGTAGACACGGATTGGTCATCTGACTGGGT",
"AAGTCTTTGATTATGCGTATTACCTTGTAGTGAGGCCAGTCA",
"TATTGTACTGTAATAGCCCAGCATGGCCGTCGCAATAAGACT",
"GGCTTTAAGACTGAGAAATGAGTTGCGCGCCTAAGGAGAAAT"]

print (minimize_score(input,6))