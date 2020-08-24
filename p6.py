#Find a position in a genome where the skew diagram attains a minimum.

def min_skew(sequence):
    skew = {"A": 0, "T": 0, "C" : -1, "G" : 1}
    tot_value = 0
    collection_skew = {}
    min_pos = []
    for nt in range(0,len(sequence)):
        if not ((sequence[nt] in skew)):
            continue
        collection_skew[nt] = tot_value 
        tot_value += skew[sequence[nt]]
    min_skew_value = min(collection_skew.values())
    for key in collection_skew:
        if collection_skew[key] == min_skew_value:
           min_pos.append(key)
    return min_pos

my_file = open("dataset_7_6-2.txt", "r")
temp = my_file.read()
print (len(temp))


# input = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"     
print(min_skew(temp))
