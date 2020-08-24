#most frequent k-mer with up to d mismatches -> inefficient but works 
import collections

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
    if mismatch == 1:   
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


def freq (arr):
    freq_arr = []
    freq_dic = collections.defaultdict(int)
    for sequence in arr:
        freq_dic[sequence] += 1
    for key in freq_dic:
         if freq_dic[key] == max(freq_dic.values()):
            freq_arr.append(key)
    return freq_arr

input = "AAGTTACTTACTCCGAAGATGCATGCTTACTCCGTTACGGATATGCGGATGGATATGCTCCGATGCATGCTTACTTACATGCTCCGGGATTTACTTACTTACGGATGGATTCCGTTACATGCATGCTCCGATGCTCCGATGCTCCGAAGTTACTTACTTACTCCGTTACTCCGTCCGTCCGGGATAAGTTACTCCGATGCTCCGTTACGGATATGCGGATTTACGGATTTACATGCAAGTTACAAGTTACGGATAAGGGATTTACATGCTTACAAGAAGATGCAAGAAGATGCTCCGAAGATGCAAGGGATATGCATGCGGATGGATAAGATGCAAGAAGGGATAAGATGCATGCATGCTCCGTTACTCCGAAGGGATTTAC"

arr1 = (all_mismatches(input,5,3))

print (freq(arr1))




