'''
Created on 22 Jan 2016

@author: steve
'''
"""
Align to seq/s to ref
"""
import operator


def align_reads_to_seq(seq_dict, ref, sRNA_length):
    """
    Return mapped reads for a single ref_seq
    pos is 5' end of read relative to 5' end of fwd strand

    alignment_dict --> read:[[pos,count]]
    """

    alignment_dict = {}

    count_start = 0
    # indv_seq_align_count = 0

    ref_complement = complement(ref)

    while count_start < (len(ref) - (sRNA_length - 1)):
        query_seq_fwd = ref[count_start:(count_start + sRNA_length)]
        query_seq_rvs = ref_complement[count_start:(count_start + sRNA_length)]
        if query_seq_fwd in seq_dict and query_seq_fwd not in alignment_dict:
            alignment_dict[query_seq_fwd] = [
                [count_start, seq_dict[query_seq_fwd]]]
        elif query_seq_fwd in seq_dict and query_seq_fwd in alignment_dict:
            alignment_dict[query_seq_fwd].append(
                [count_start, seq_dict[query_seq_fwd]])
        if query_seq_rvs in seq_dict and query_seq_rvs not in alignment_dict:
            alignment_dict[query_seq_rvs] = [
                [len(ref)-count_start-1, 0 - seq_dict[query_seq_rvs]]]
        elif query_seq_rvs in seq_dict and query_seq_rvs in alignment_dict:
            alignment_dict[query_seq_rvs].append(
                [len(ref)-count_start-1, 0 - seq_dict[query_seq_rvs]])
        count_start += 1

    return alignment_dict


def aln_by_ref_pos(alignment_dict):
    """
    Create 2 lists - fwd alignment and rvs alignment
    Each contains tuple(pos,count)
    Returned as ordered.
    """

    fwd_alignment={}
    rvs_alignment={}
    aln_count = 0
    for alignment in alignment_dict.itervalues():
        for i in alignment:
            if i[1] > 0:
                fwd_alignment[i[0]]=i[1] # pos,count
                aln_count += i[1]
            elif i[1] < 0: 
                rvs_alignment[i[0]]=i[1] #pos:count
                
                aln_count -= i[1]


    sorted_fwd_alignment = sorted(fwd_alignment.items(), 
                                  key=operator.itemgetter(0))


    sorted_rvs_alignment = sorted(rvs_alignment.items(), 
                                  key=operator.itemgetter(0))
#    detect_phase(sorted_fwd_alignment, sorted_rvs_alignment)
    print "\n{0} reads per million reads have aligned\n".format(aln_count)
    print "-"*50
    return [sorted_fwd_alignment, sorted_rvs_alignment, aln_count]


def count_align_reads_to_seq(seq_dict, ref, sRNA_length):
    """
    Return mapped reads for a single ref_seq
    pos is 5' end of read relative to 5' end of fwd strand

    returns an integer
    """
    #start = time.clock()
    
    aligned_count = 0 #number of reads aligned

    count_start = 0
    # indv_seq_align_count = 0

    ref_complement = complement(ref)

    while count_start < (len(ref) - (sRNA_length - 1)):
        query_seq_fwd = ref[count_start:(count_start + sRNA_length)]
        query_seq_rvs = ref_complement[count_start:(count_start + sRNA_length)]
        
        if query_seq_fwd in seq_dict:
            aligned_count += seq_dict[query_seq_fwd]
        if query_seq_rvs in seq_dict:
            aligned_count += seq_dict[query_seq_rvs]
        count_start += 1

    return aligned_count


def list_align_reads_to_seq_split(seq_dict, ref, sRNA_length):
    """
    Return mapped reads for a single ref_seq
    pos is 5' end of read relative to 5' end of fwd strand

    returns an integer
    alignment_list --> [sRNA,sRNA,]
    
    """
    count_start = 0
    ref_complement = complement(ref)
    alignment_list = [] #aligned sRNAs
    while count_start < (len(ref) - (sRNA_length - 1)):
        query_seq_fwd = ref[count_start:(count_start + sRNA_length)]
        query_seq_rvs = ref_complement[count_start:(count_start + sRNA_length)]
        if query_seq_fwd in seq_dict:
            alignment_list.append(query_seq_fwd)
        if query_seq_rvs in seq_dict:
            alignment_list.append(query_seq_rvs)
        count_start += 1  
    return alignment_list 


def complement(seq):
    """Provides the complement in the 5' - 3' direction

    Assumption: reference consists of A, G, C, T only

    complement(str) --> str
    """
    d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return (''.join(d[c] if c in d else c for c in reversed(seq)))



