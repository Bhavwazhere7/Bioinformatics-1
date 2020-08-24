# Determine total number of occurrences of Pattern in Text with at most d mismatches

def humming_distance (seq1, seq2):
    count = 0
    i = 0
    while i < len(seq1):
        if not seq1[i] == seq2[i]:
            count += 1
        i += 1
    return count


def approx_occurances (seq,pattern,mismatch):
    count = 0
    pattern_length = len(pattern)
    for i in range(0,len(seq)):
        temp_approx = seq[i:i+pattern_length]
        if len(temp_approx) == pattern_length:
            act_mismatch = humming_distance(pattern,temp_approx)
        else:
            break
        if act_mismatch <= mismatch:
            count+=1
    return count

seq_input = "GCAGTATAATGAGTTCTAACTAGGTGATAAGAGAGTTAGCTTAATAGTCTCTCAGCGTGTTGACAAAATACCTCACGTTTCAATGGTGAACACGTCTATATCGCTGGCGAACCGTAAAGGGTAGCGTATCGCTTAAAGGCGTATTGGTGCCCGGTTTTGGAGGCTGCACGAAGAGTATTCAGTGAGAGTTAAAGACGATTCCCGATAACGTACCTGCTTCCCGCGGGAGGGATGTCCTGCGTGAGAATGGCTTTCTGCAGCGCGGATGGGGTCTGAAGGCCCTTGCTCGACCCTTCGGATCCCCATGGTCAATTGTGCGTCCTGCTAGCCTGAGCCGTCGTCACATTTGCCCTAGTAATGGCATAGATATTGTCTGTCACACATAACC"
pattern_input = "ATAACC"
print(aprox_occurances(seq_input,pattern_input,3))