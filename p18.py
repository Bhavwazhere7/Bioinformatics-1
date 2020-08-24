#Gibbs Sampler
def all_poss_kmer (seq,k):
    kmer_arr = []
    for i in range(0,len(seq)):
        new_seq = seq[i:i+k]
        if len(new_seq) == k:
            kmer_arr.append(new_seq)
    return kmer_arr

def score (seq_arr,k):
    score_arr = []
    dic_nt = {"A":0, "C":0, "G":0, "T":0}
    for i in range(0,k):
        for kmer in seq_arr:
            dic_nt[kmer[i]] += 1
        score_value = abs(max(dic_nt.values()) - sum(dic_nt.values()))
        score_arr.append(score_value)
        dic_nt = {"A":0, "C":0, "G":0, "T":0}
    return sum(score_arr)

def matrix_creator(seq_arr,k):
    dic_nt = {"A":0, "C":0, "G":0, "T":0}
    k_dic = {"A":[],"C":[],"G":[],"T":[]}
    for i in range(0,k):
        for seq in seq_arr:
            dic_nt[seq[i]] += 1
        for nt in k_dic:
            k_dic[nt].append(float(dic_nt[nt])/sum(dic_nt.values()))
        dic_nt = {"A":0, "C":0, "G":0, "T":0}
    return k_dic



def gibbs_sampler