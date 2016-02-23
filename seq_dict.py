'''
Created on 22 Jan 2016

@author: steve
'''
"""
Module for loading sRNA seq files in fasta format
"""
import time

# a single seq dict


def load_av_seq_files(seq_file_1, seq_file_2, 
                      sRNA_max_len_cutoff, min_reads,sRNA_min_len_cutoff):
    """
    load 2 seq dicts in .fasta format.

    Calculate RPMR and apply in the function

    sRNA_len_cutoff --> only reads  of length <=
    sRNA_len_cutoff accepted (including for RPMR calculation)
    min_reads --> only reads >= min_reads accepted

    produce seq_dict --> sRNA:av_RPMR
    """
    start = time.clock()
    seq_dict_1 = {}
    seq_dict_2 = {}
    read_count_1 = 0
    read_count_2 = 0
    loaded_seq_1 = open(seq_file_1, 'rU')
    loaded_seq_2 = open(seq_file_2, 'rU')

    av_seq_dict = {}  # This is the combined seq dict

    # Complete dict for seq_1
    for line in loaded_seq_1:
        if line[0] == '>':
            count = int(line.strip().split('-')[1])
            next_line = True #ensure that the nextline has the correct seq.
        elif count >= min_reads and len(line.strip()) <= sRNA_max_len_cutoff \
        and len(line.strip()) >= sRNA_min_len_cutoff and next_line == True:
            seq_dict_1[line.strip()] = count
            read_count_1 += count
            next_line = False
        else:
            pass

    loaded_seq_1.close()

    # Complete dict for seq_2
    for line in loaded_seq_2:
        if line[0] == '>':
            count = int(line.strip().split('-')[1])
            next_line = True #ensure that the nextline has the correct seq.
        elif count >= min_reads and len(line.strip()) <= sRNA_max_len_cutoff \
        and len(line.strip()) >= sRNA_min_len_cutoff and next_line == True:
            seq_dict_2[line.strip()] = count
            read_count_2 += count
            next_line = False
        else:
            pass

    loaded_seq_2.close()

    # Av sRNA RPMR if sRNA present in both
    for sRNA, count in seq_dict_1.iteritems():
        if sRNA in seq_dict_2:
            av_seq_dict[sRNA] = ((count * (float(1000000) / read_count_1)
                                  ) + (seq_dict_2[sRNA] *
                                       (float(1000000) / read_count_2))) / 2

    print "\nSequence file loading time = "\
     + str((time.clock() - start)) + " seconds\n"
    print "{0} has {1} loaded reads\n".format(seq_file_1.split('/')[-1],
                                              read_count_1)
    print "{0} has {1} loaded reads\n".format(seq_file_2.split('/')[-1],
                                              read_count_2)
    print "-"*50
    
    return av_seq_dict


def load_triple_seq_files(seq_file_1, seq_file_2, seq_file_3, 
                          sRNA_max_len_cutoff, min_reads,
                          sRNA_min_len_cutoff):
    """
    load 3 seq dicts in .fasta format.

    Calculate RPMR and apply in the function

    sRNA_len_cutoff --> only reads  of length <=
    sRNA_len_cutoff accepted (including for RPMR calculation)
    min_reads --> only reads >= min_reads accepted

    produce seq_dict --> sRNA:av_RPMR
    """
    seq_dict_1 = {}
    seq_dict_2 = {}
    seq_dict_3 = {}
    read_count_1 = 0
    read_count_2 = 0
    read_count_3 = 0
    loaded_seq_1 = open(seq_file_1, 'rU')
    loaded_seq_2 = open(seq_file_2, 'rU')
    loaded_seq_3 = open(seq_file_3, 'rU')

    av_seq_dict = {}  # This is the combined seq dict

    # Complete dict for seq_1
    for line in loaded_seq_1:
        if line[0] == '>':
            count = int(line.strip().split('-')[1])
            next_line = True #ensure that the nextline has the correct seq.
        elif count >= min_reads and len(line.strip()) <= sRNA_max_len_cutoff \
        and len(line.strip()) >= sRNA_min_len_cutoff and next_line == True:
            seq_dict_1[line.strip()] = count
            read_count_1 += count
            next_line = False
        else:
            pass

    loaded_seq_1.close()

    # Complete dict for seq_2
    for line in loaded_seq_2:
        if line[0] == '>':
            count = int(line.strip().split('-')[1])
            next_line = True #ensure that the nextline has the correct seq.
        elif count >= min_reads and len(line.strip()) <= sRNA_max_len_cutoff \
        and len(line.strip()) >= sRNA_min_len_cutoff and next_line == True:
            seq_dict_2[line.strip()] = count
            read_count_2 += count
            next_line = False
        else:
            pass

    loaded_seq_2.close()

    # Complete dict for seq_2
    for line in loaded_seq_3:
        if line[0] == '>':
            count = int(line.strip().split('-')[1])
            next_line = True #ensure that the nextline has the correct seq.
        elif count >= min_reads and len(line.strip()) <= sRNA_max_len_cutoff \
        and len(line.strip()) >= sRNA_min_len_cutoff and next_line == True:
            seq_dict_3[line.strip()] = count
            read_count_3 += count
            next_line = False
        else:
            pass

    loaded_seq_3.close()

    # Av sRNA RPMR if sRNA present in both
    for sRNA, count in seq_dict_1.iteritems():
        if sRNA in seq_dict_2 and sRNA in seq_dict_3:
            av_seq_dict[sRNA] = ((count * (float(1000000) / read_count_1)
                        ) + (seq_dict_2[sRNA] *
                        (float(1000000) / read_count_2)) + (seq_dict_3[sRNA] *
                        (float(1000000) / read_count_3))) / 3

    return av_seq_dict


def load_seq_file(seq_file, sRNA_max_len_cutoff, min_reads, \
    sRNA_min_len_cutoff):
    """
    load 1 seq dict in .fasta format.

    Calculate RPMR and apply in the function

    sRNA_len_cutoff --> only reads  of length <= sRNA_len_cutoff
    accepted (including for RPMR calculation)
    min_reads --> only reads >= min_reads accepted

    produce seq_dict --> sRNA:av_RPMR
    """
    start = time.clock()
    seq_dict_1 = {}
    read_count_1 = 0
    loaded_seq_1 = open(seq_file, 'rU')
    final_seq_dict = {}  # This is the returned seq dict

    # Complete dict for seq_1
    for line in loaded_seq_1:
        if line[0] == '>':
            count = int(line.strip().split('-')[1])
            next_line = True #ensure that the nextline has the correct seq.
        elif count >= min_reads and len(line.strip()) <= sRNA_max_len_cutoff \
        and len(line.strip()) >= sRNA_min_len_cutoff and next_line == True:
            seq_dict_1[line.strip()] = count
            read_count_1 += count
            next_line = False
        else:
            pass

    loaded_seq_1.close()

    # final RPMR - could simplify in future
    for sRNA, count in seq_dict_1.iteritems():
        final_seq_dict[sRNA] = count * (float(1000000) / read_count_1)
    print "\nSequence file loading time = "\
     + str((time.clock() - start)) + " seconds\n"
    print "{0} has {1} loaded reads\n".format(seq_file.split('/')[-1],
                                              read_count_1)
    print "-"*50
    return final_seq_dict

def load_seq_list(seq_list):
    """
    Load seq files names/paths from file
    """
    seqs=[]
    loaded_list = open(seq_list, 'rU')
    for line in loaded_list:
        seqs.append(line.strip())
    loaded_list.close()
    return seqs

def load_av_seq_list(paired_seq_list):
    """
    2 Seq files on one line deimited by tabs are replicates
    """
    seqs=[]
    loaded_list = open(paired_seq_list, 'rU')
    for line in loaded_list:
        seqs.append((line.strip().split('\t')[0], line.strip().split('\t')[1]))
   
    loaded_list.close()

    return seqs


















