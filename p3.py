#Reverse Compliment

def reverse_compliment(sequence):
    nucleotide_pairs = {"A":"T", "T":"A", "C":"G", "G":"C"}
    reversed = ""
    i = len(sequence) - 1
    while i >= 0: 
        reversed += nucleotide_pairs[sequence[i]]
        i -= 1
    return reversed

input = "ACTACAGGAAATGACTTGCATGAATTAGCCTGTCACCGTGAAATGAGAATGTGGGTACCGTAAAGCAAAAGGCACAGATGATGTATAACACTGGTAGCCATTCAGCGACGACACTTCACATTTGGCCTCTCAGAGATCCAACGATGGTACGGTTTAGGCGCCTAACTATTCGGTGTTGCGTACCGAGTGTGCCCGGACGAAGGATTTGGACAGGCTATTACGAAATGAATCGTTGGCCCAATGACAGAACTACTCAACTGCCGGGAAATGGTGGCGTCAATTGGATATTCATACCCTCTTTGGCGCTCACGGCCCAGGCATATCTCTCGTAGTACCGCTCAGATAATAGTTGTGGTCTTTCACTAGCAACGCAGAGTGAGGTGAGCACACGCCTAGCGTGTGGGGCTTGCGGTTGTATCAGGGCACCGTCGAACCTCCTCGACTGGCCAGGCAAGGACACTCTGTATTCCGGATGACTAAGATGGGGGCTTAGGCGCGCACGGACACGTGCGACCACTGGGGGTCAACCGTATAAGCCACTGTAGGCTCGTGTGTCTCTGGTAAAAAGACTCAGCGGATGGTCAGTTAGTAAACGACCCCGATCTAATGTATCCCCAGCTGAATGAAGCCAAATATGACGAGAGTGGTGAAGCGTTAAGTTGTAAGCACCTAAGCCAGATTCGGATCTTTGTACACTGTGAGGTCGTATTGGCACATGAAACTCCTGAGAAGTCCCGGTTCACCCTGCAGAAACTGTAGCCAAGCCGTGCCGTAGGCATCCTGCATCCATGCCCATAAATAGGATCTTTAGTTAAGCCAATCGACCCCCTAACATGTGACGACCTTAAGTGTACTCTTGAGCTACTGCTGGCGTCATCAACATCCAGACCCTGGAGAATGCTCCGTTTGTTGACACAGTCCTGCCATACCGGCGGGCAGTATCCCGCAAGAGTCCCATCTCCTCGCGTACGGTACGAACACTAGGTAGGAATCTAGTGGGCGTGGTCCTTCTAAGGACGCAAAACACCAATTGCGCGGCCTAGATCTGGCTTTACAAGCTATTATGTACCTCAAATAAAAAGAGGATAAGGGCCGGGCCAGCGGCTCGACAGACGCCAGGGGTAGGTGCCTGCGCCCTCCAAAACGTGTCTCGGGTATGTAATACTGCAAGCCGGACGCACCCTCGTCCAAACATCGACAATGGCGTAAGACCGCGTAACTCACAGTGAGTGGCCCTACAGGCGTAGGAACTCTTTCTCCACATGAGGAGTGTCCCGCGAGTAGCAACTTACGGTCGTTTCGTCCTCTATCCACTCACATGTGCTTAGCACCGATGCTGTAGGGAACACGAGCAATCCAAGAAGTGACAAGGTGGGTTTTGCCCACAAGTGGTTGGCACAGATGCCCAATTACGATGATGTCGAATTTTGACTCATCAACAATGGGCTTTGCCACTAGGTTGCTGCCGCTCGCGCAACCGACAGGAATAATGGTCGTCTCAGAAGGGTCTGAATCGCTCACCTCATGCAAGTTAGAGACAAGAAACGTCCTTCCCCGACGAGACGGACGCGTCCAGGTTGAAATATACTGGTTTCACGTACGCATTCAGGTACAGTTGGTTATCCTTACTCGGGTAGCACCCAATTTATGAGCCGGGTCGCATGAGAACACGGACCACTCAGTCACGCTCTACAATTAGTGCGGTAGTTAGGAGAACGCCTGTATTCACCGAAATTAACTTGATTGAACACAAACTGACGTAACTGCATCGCTGATATATCCATTGACGTATCGTCGCTAGCGATTTAATGCACTGTTCCCGAATGTAAGACGAAGTGCTGCTGCTTGCTAGCGAGGAAGCTTGAACAGTCCTGGCCGCCATTTCAGTCATTGCAACAAGGTGATAGCATCGTCATGCACGGAATGAAAGTTAATCCGATCGAATAAGGAACGAGTATGCTCGACGGGATATGATAGACCTTATTACCAGCTCGACCGCGTGGGTTGATGAACGTTGGGGTGATGAGTGAAGGGGATCCAATAACTATGATAGCCCAAGTCTCCGGTGACGGGTGCTGACTCGCACAGGCATGGCCTCTTCATTTCGCAAAGTTGTCTTGTACAGCGTCTTACTACGTTTTGCGAGGCAGATACCCGAGAATTGACTCTTCCTGCGGGAGCCTAAAATGGGCGCAAGGACCGTGTCCGAGTTCGAGAAAGCGGGACCACGCACGCAGGGCCATCAGGACTGAGTCCAGCTCGTTGACTTATAACCGAGGGAGCACTTGATTGCTATGTGCTCTGCCATAAATTCCAGTGAGGTCTTAACACACACCGAGTCCAGCCTGGAAACACAGCGGCTCGAACTGGAAGGTACTGGGGTTATACGTACGCTAAAAAGTAGTGTCGGATGGGTCAACACGAAAAAGGCCTGTTTCGCAGAATTCAGTGCGCAAGTCTTTCGGAGAGTAAAGTTCCTATTTGAAGTGGATCTGGAATAACTTTCTAGTCGCACCCCCGTTCGTACACATTACCAGTTACTCGTAAATGCGCGCATAGTAGGGCCAAATGCTTACTAAGCTTGTACGAGTTTTACAATAATTTAGCCCGCAGTTAGAGGGTGATGCCTTTACTTCTCCCCGTCGCTTTGTATAATGAAGCTCAACATCGCTTAAAGTTTGCAAGGACAGGATTGGATTTAGAGCTACCATGCTATAACCCAATAGTTACGCCATTTCTTCAACATATGTTACTCCCCATGTCCGCTTACTAGTCACTGTATATGAAAACTTGAATCTACTGAGTGCTGGCTCCGGGATCACTTACTTATTCCGGCCCTCCCCGTGTCAGTACACCTGATACCACTGAAGAAACATGTTTGCCAGAACGTATGCCTTCTAGGGACTAAAGGCGCAGCTCCTGGGTCAAAGTCTCCGTCCGAGCAATCACAGATAGATAGTTTGGATGGATCGAGTGTCCACATCTGGCCGACTACACTAACCCCTCTCCCTGGATACAATTACTCGACATAGGTGCTTAGGACGTGGGCTCTTCGTGTGTTAGACAAAAGTGCAAGCGACCATATAAGCCAGCTGAGATGACTCGATGGTTGACCAAACTCCGTACGGAGTATCCTTACCCGACTTTAAGCCTTAGCTTAGAGGACATAAATGCCTTCTAGCCGTTGCTCATATGGGTGGTCATTGCGGCAACAGGAGTCATGCGAGAATTCCAATAGATCCAACATAGCCTTTCCACCGGTGACGCATAAGAGTCGACAACTCCCCGGAGGGGGGAGTTAACTTCCCGACGACCGTACCACCGATATCGACACAGCCACCTCGGTGGATGGCCTCTAGTCACTTCGTTTTGCGTGGAGGCAATTACTTTCGTCCAAACATGGTGTCTATCCTAGCTGATTCGCTACTTGTTTGACGTGTAGCAACGGCAAGCATTCCGCTGTCTTAGCGAAGGGTCCGTTCAGAACAAACAAGTAGGAGCACCAATGCGTAAAAACCGGTGTTATAGAGCAACAATTAGGGTCCTATCGATGCTAATCGCCACGTAGGGATCCTGTTGAATTCACTCGATTGGGTGTGGTATCCGACCTGAGCGGCGTTATGCAAGACGTTGGCACTAAATCTGAGCTAGGATGGATACCACCGAGTAACTAGACTCAGGGTTATAAATAAGGCACCTGTCACAGGATACAGCATATGCTAGCTCCACCCCGGTGGGATCATAGTGACAACAGCGTGTTTTCATCCGCTACGCCAATATTGGCAAAGTGTGCCGTACGTTATTCCCACCTCGGTACACATAACTGAATTGCCAACCTTAATCCCATTTGAGCTCAAGATATGGTCTTAGGACTAGGAAGAACCGTACAAATTAATGACGAATCTACGGCCCGTCGTCGGATCATTGACCGTGCACCAATTTAAGGTAGTGTCATTCTGGTAGACGATAGGTGGCTCGATTCAACAAAAGCGCTGAATGTTTGTTAAACCCCCTATTCCATTCGTCCCGCCACTCTGCGTATCTCCCTGGTCGAGACACCACATAACAAAATCTGAAGAGTATGCGTGTTGAAACAGCCTTCTCAGAACGCCTTATAGATGAAAAGGAATGCCCGTAGGGGGATCTCACTCACGCTTTATAATGACGCGGTTCAACTGAGAAAAGCGGAACTTTGCAGTTGAGTGCCATAGACAGGGTATCTATTGATGCTAATGGCTACGGTGCTTCCGGGCCCGTCAAAGGCGGAGCTAACATTGCCTGGAAGCGTTTCATCTGATGAGGGCCCGAAACTTCCATGCTGGATAGTGTGACCCATGGCCCTCGGTACGGTAGACGATGCCATGGCACATGATGTAGGACTCATCAATGAGTAATTAAATGCGGAGTATTGCGGCCCGGGGTTCGGTACTCGAATTTCACAGGTAGTCCTACATTCAGGTGCGAAATCGCTTCATTCCCGTTCTTTCCGATAAATCTATCCTCGTTTCTATAGGGTCACCTCATGACTCTATAGCATACCCAATAAAGGATTTGCGACTCTATATCTGACCGCCAGACGGGCTCAGAATTTGCTAGAGTTTTATGCTCTCTGTGATCCGCTATCACAAGACACGGATCTGGCGAACAGCGAAATCTCAGAAGGTTCTTTTCCACGCTCTCATATAGAGGGTCTACTTTATAATACAAGGGATCGTCAGAAGCGTTGTGAGAGATTGCTGAAAAGCAGTATAGCAGAATGTCTATTGATTAGTTTTTACGGGTGTTGATATCCTACAAAAAAGAGTAAGCCGCCGTCGATCTAACGAGCCGGGGTATTGAGCCAGAATTTGTCGTAAGAGTGAACCAAGCCGCGACGTCCGACAAATTAAAAGTACACGCCTGTAGGGAGCCTCGCTAGATTCCTCTCGATGGGCTCGATCTCTAGAACATATTAGCTGTGTTCCTTCTCTCAAAAAAGAAAATCCAAGTGCGCCTATTCAGATTTAACAGGGTAACAGTATTGCAATGCGAGTTAAAGCGCGGGTTGGTACCTCACCTAAACGAGAACCAAGTGCTTTCGTAATCATCTAGCCTGCCGTTTATCTGTAACCATCGACAATGATATTCTGGTGTCCATAGGGGAGTGAATGACAATCACTCGTTTAGTCCGCTTCCGCAACGAGCCCTGCTGAAAATGGGGGAAATTGTCCTCCGCGCAGCATTGATGACCAGACCCGCATTGAGATGGGTAAGCTAGTGGCGCGTCATATATGGTCAACCTAGAGGGCTAGGGTCTCATTTGGCTGAACGTCCATGCCTGCTCGCTTAGCAAGATTTCATAATAGTTCTATACCGACAACAGTACACTCAGGCCACCAGCTTCCTCGCGTATGCTGCCCTGTTGGCAACACTTAGCTCTACTTGGTTGCTGAACTTGTGAGCAGGTCAGATGAGAATAATAAGGCTAGCGAAGAGATCACGTAATACCCATCCTCTCAAAGTTGTTAGGTGGGGCGCCACCTAACATAAGGAGGACTCCCCCATCGTTCCAGCATCAAAATGCGAGTCGTCTATGGCTGGGACATGAGGCAAGTCAATATCTAAGTGATTAGCTTGCAATAGCCAGAGGTACTTTAAAGTACCGGAATGTGTATTGTATTCTAGCAAGCAATTGTGAGAGCCAAGGAGGGGTTCTATTGGTTCTCAGTAGAACATTTTCGCATTGACGTGGTCTCCGTCCTTCGAGGTGGCTATCTCTTGCGGGACTTCCGGCAACATGAGTCTTCCCTGATAGCCAGCGCTATCATATCCTTATTGGACACTAACAGGACCTTGAACGTGGCCCTTCCAAGTTCCTGAGGTGGTTCTGTCGGTGGTGGGGGCCGACTTATGGTCCGGGGAGGATACCAACCGATTGCCTTGGAGAAGGGTTCTATATGGTCGTTCAGATTTCCTAGGTGCTAGGTCCAAACGACCACCCTCCGTAAGGTATCGTCACCGTGTGTCAATCGGCGGCCCTGTTATCGCCCGGAGACTTTCAAATCGGGGTGGATGCATACACCGACGCAAGCCGTATAGACTGTTTACTAGGGGACGATTGAAGCGTGCCGTAACGATTACGTGGATTCTCTAAACACCTGAGCATGCTCGTTGAAGCTGTTTGCAGGACGGGGCCATGATACTGTCGCTTCATGTGCACAGACCTCCATACGGTCTCACGTATATTGAAAGTCAGGGCTTAGTCTGGTACTTGGTGGGCGTCCCTGATTGCGGTAATGCGATTTAACTGCAATCCGCCCCATTCTAAACGTGATCGGGAATCGATGTGCCAGGGCGTGTATGCTCGGCTGCGTCGCACGAGACGCTGCCTTAGAAGTCCTCTAACGTCCGTGTAACGCTCTAGACTTATCACCCGTCAATAGGGTGTTCTATCCGTAAACTCACTGTACGGCATGTACCAACTTTTAGAAATGTAGCAGGCAGCATATTGAAGCCTAGCTAATTCTACAGAGGATTGTCGGGCCATCGCGTGTCTCACAGCACCGGACAGACGTGTAGACTCATCATGTCGAAGTCCTTCATGGGCGGTAGTTCAGAGCGGTTCCTTGACCCTCCGACGCCTGGTACCCACCTCAACTTTGTCCGTAAAATGGTTTACCGCTCAGGAGGATTGGTGCAACCGGTGCTCCGTAACTGTAAAGTAAGTTCCCACACGTCTCAGGCTGCAAAACTAAACGTAGGGACGGACGTTACCCAATTTTGGCCCGGTGTCTGGTAGTAGGACCGTGGCGTAACGATGGGGTTCACTTCATCGTCTGGAGCTTATCGCATCTGATTATTCAGTAGACCCCAGGTCTATTGTAACAACCGGAAGAACCAGGTTTTGCTTGCCCTTGGTCAGTTTTCTGAGATAGCAATCCGAAATGACAAATACCATTTCCAATGGCCAAGGAATACATGTTTGATGACCTGCTTAGTACGTCATTTACTCACCAGAGTTTCTAGCCTCCTCGTGAACGCGACGCAGTAGCGCCGTTTCTGCCCTCCGCAGTAGTCTGCGGCATCCATTCTGCCAAAGGATCGGAATAGATCAGGTCTTCAATACTCGGATTGCCAGGCTGGGGAGTTTTATCGGTGAACGTCCTCTGGGCATGCCAGCGTTGAAAAAGACAGAGCGAAAGCAATTATATTACAGTTTATAGACAGATGGTCCTTAAGTAGATTTGGCAAAAAAGAGAGTTGAAAGACGGGCCTTAGCGTGGTCAGGTAGGGCCCATGAAGCTACTTCTAATACAGCGGAGCAGGTTCCTAAGACGGAGTTCCTATCTACAGAGGCACTGCCATTATCTGAGCAAAGCGTCGAGGAGTGTCTGTTCCAGAGGGGAAACTAGACCGCCATATACATCTTCTAGGAATATGATGCTTAATCCGTTGTCTTATAGCGGTATATGCAGTCACGAGGCCGCTAAACTCACTAGACGCATCCACAACGGTGCGCCGCATGGAACGGACATACGGGGGTGTGGGAGTTAGGCCAGGAAAGCAGCACTAGTAGGAACTGATGCTGTGAATGATGAATGCCTGGAGCTGGACGTTAGGCTACATTGAGGCGCAGCAGGCGATCCAAGTAGATCTGCTCGGCTGTAACTAGTCACACCACGGCTTTGGATCTTTGGTGTCCGAGGTATATAGAGCGTGCTTTAGCTATTAACACAAGCTCTTGCCGTACACCACGATCTCTCTCACTCTCAAGATATTATGTAGTACCCATTGCGCTTTATTTAGCTTAGGATATTGAGATATGTCATCTCTAGAGCTTATAACGGTGAACGGTGATGCCCTCGGTATCACGCGGGGGTTCATTAGCGATTACATTGCACGAGTGCTTGCTATACTTAAGTCTGCTATGGCGCGACCTCGGGCGTGTACCATATCAATCCAATACGGGTGGTGTAATCCAGGGGGATGCTGCCGGTAAAAACTTATTTCCTGCTAAAAGTAGGATCCGCTCTATGGAGTACTCGGGCGCATGAATGAGACGTAGTTGAGATCTGCATGTGGTTAGGCGATTAGCTGCACCTATTTCCTAGGCCGAGTTTCTCTTCCCAGCTTCTAACGCTCGCCGCGCGCCCGCTATGCACTTGGATTTATGTTTGACCAGGAAAATCTAGAGCTCGACTATCACCTTGTGTCATGGCCCGAGCTCAAAATTCGTCCTCTCGCTGAGCAGCACGACACACCCTACAGCATCGTGTTGCAATAACTAAGGGGAGTATCTGGCGTCACGGATCATTTCTCGGTACGACGCAACCAAGATCCTGGATGGTCGTAGCTGTCCTTACATTCTCCATAGTCTTGCAACAGGTTATGTCCATTTCCACGGACCATTAACACTAAGGACTCGATGCACACGTGCATCGGCACAGGGCCCGTAACTCCGGTATAGGAGCTAACCCTAGACTTTCTTGTACAATGGTGGATGCGATAATTCCGTTCCGCCTCCTCTTCGGGATCATAATTGAACCCCAGTTTAGTGAGGCCCCCCATAGTAATCGTATGTGGCTACGGAGCACACATAGAATCTCCAGCTTGAGCCCTACCGGGCAATTTCCTTCGTAGGTTGGGCCATGCGAGGCCGTTGGAGAATAGGAACCTGACTACTTATGACTTCCCCATTTCTAGAACAGATCTGCTCAGACAAGGCTGTCACACGGGCGGTCCTTCCCTAGTCTTCTCTTGATAGTGGGCGGTAATCTTCTACGAAAGGGACCACACAACTAATGGAGGGGACAAACCCCTTCTCTACCCCCATTCGCTCTAATTAGCATGGAACCGTTCAACTGGGCGTATCGTCTGCATCTTATAAAAAAGGCAATTGTGCGACAGAGCAGGGCATCATCCCCTCGTGCTCCTCCATACCAATTGCGTGGAAGCGCCGCGAGGAAAATGCATGATAACTAGGGTCCCTACAATTTATTCGCGAGATCCAAAGTCGGGTCATAGTACATGATCCGAGGGTAGATTAGCGAATGAGAACCAATTGGCAAGATCAGAAGGGAACCAACGCAGTGCGGTCTACAACACCAGGAACTCTCAGCCAGTTTTTGCATCAGTCACAATGAAACGCCTAAAGCCGTGATTCGTGGTAAAGATTCTTGAGATGACCGTAGGGAGCATAGGAATACCGCATCCTTTTCCGTCAGGACAAATCACAAGGGGCACGGCAGATTACTAAGTTCGCAGTGAACTCTGTGCTAGGGAAACGGGGCATCGCAAAGTTGCCTCTCCCTCACCATGTGCTCCCTGGGCAGCAGGAAGTCCAATGAATATATGCCAGAGGGGGTTGCAAGGCTTCAGCACGTTATCATTCTTATATGGCCGGATCTGGTCGCACAGGTGAGGTGGTATGGCTACTAGCGTGTGATATGATAATGGACGGTAAGTATCCCTAAGATTGCTGCACCAGCTCAAGTCACTGGCTAGAATGGGGCCAACCTGTTCACGTTAGTTCGTTTGAACCAACTTTTACCACAGTAGACCAAGGGCTACAGCATGACAAATAGTAGTCTACCCCTGTTAAATAATACTATTATTATCCGCGCTATTGCTATAGATACATTTGCTGAAATTGAATAAGTGAGACGTGTCTGTTACGCTACACCTTGTCGCCCG"

print (reverse_compliment(input))