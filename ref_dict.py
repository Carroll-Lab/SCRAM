'''
Created on 22 Jan 2016

@author: steve
'''
"""
A module for loading reference sequences in fasta format
"""
import time

def load_ref_file(ref_file):
    """
    Load ref file in fasta format

    Product ref_dict --> header:sequence
    """
    start = time.clock()
    ref_dict = {}
    ref_count = 0
    loaded_ref = open(ref_file, 'rU')
    full_len_seq = ''
    key = ''
    first_header = True
    for line in loaded_ref:
        if line[0] == '>' and full_len_seq == '':
            key = line.strip()
            ref_count += 1
            if first_header is True:
                first_seq = key
                first_header = False
        elif line[0] == '>' and full_len_seq != '':
            ref_dict.update({key: full_len_seq})
            key = line.strip()
            full_len_seq = ''
            ref_count += 1
        elif line[0] == '' and full_len_seq != '':
            ref_dict.update({key: full_len_seq})
            key = line.strip()
            full_len_seq = ''
        elif line[0] == '':
            pass
        else:
            full_len_seq += line.strip().upper()

    ref_dict.update({key: full_len_seq})

    print '\n ---- {0} reference sequences \
loaded for alignment ----'.format(ref_count)
    print "\n{0} length = {1} bp".format(ref_file.split('/')[-1],
                                         len(full_len_seq))
    print "\nReference sequence loading time = "\
     + str((time.clock() - start)) + " seconds\n"
    print "-"*50
    return [ref_dict, first_seq]
