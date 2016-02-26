'''
Created on 25 Feb 2016

@author: steve
'''
import align
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

    ref_complement = align.complement(ref)

    while count_start < (len(ref) - (sRNA_length - 1)):
        query_seq_fwd = ref[count_start:(count_start + sRNA_length)]
        query_seq_rvs = ref_complement[count_start:(count_start + sRNA_length)]
        
        if query_seq_fwd in seq_dict:
            aligned_count += seq_dict[query_seq_fwd]
        if query_seq_rvs in seq_dict:
            aligned_count += seq_dict[query_seq_rvs]
        count_start += 1

    return aligned_count


def dict_align_reads_to_seq_split(seq_dict, ref, sRNA_length):
    """
    Return mapped reads for a single ref_seq
    pos is 5' end of read relative to 5' end of fwd strand

    returns an integer
    alignment_dict --> [sRNA:align_count,sRNA:align_count]
    
    """
    count_start = 0
    ref_complement = align.complement(ref)
    split_alignment_dict = {} #aligned sRNAs
    while count_start < (len(ref) - (sRNA_length - 1)):
        query_seq_fwd = ref[count_start:(count_start + sRNA_length)]
        query_seq_rvs = ref_complement[count_start:(count_start + sRNA_length)]
        if query_seq_fwd in seq_dict and query_seq_fwd in split_alignment_dict:
            split_alignment_dict[query_seq_fwd]+=1
        elif query_seq_fwd in seq_dict and query_seq_fwd\
         not in split_alignment_dict:
            split_alignment_dict[query_seq_fwd]=1
        if query_seq_rvs in seq_dict and query_seq_rvs in split_alignment_dict:
            split_alignment_dict[query_seq_rvs]+=1
        elif query_seq_rvs in seq_dict and query_seq_rvs\
         not in split_alignment_dict:
            split_alignment_dict[query_seq_rvs]=1
        count_start += 1  
    return split_alignment_dict 

def times_read_aligns(split_alignment_dict):
    """
    Dict. of times a read aligns:
    {read: times aligned, read:times aligned} 
    """
       
    sRNA_align_counts={}
    for aligned_sRNAs in split_alignment_dict.values():
        for aligned_sRNA, count in aligned_sRNAs.iteritems():
            if aligned_sRNA in sRNA_align_counts:
                sRNA_align_counts[aligned_sRNA] += count
            else:
                sRNA_align_counts[aligned_sRNA] = count
    return sRNA_align_counts

def split_reads_for_header(split_align_dict, split_align_count_dict, seq_dict):
    """
    Dict -->Even split read aligned counts, so total reads aligned = read count
    in original seq file
    {header:split_count}
    """
    
    header_split_count = {}
    for header, sRNA_dict in split_align_dict.iteritems():
        header_split_count[header] = 0
        for sRNA in sRNA_dict:
            header_split_count[header]\
             += ((seq_dict[sRNA]/split_align_count_dict[sRNA])\
                 *split_align_dict[header][sRNA])
    return header_split_count

def header_x_y_counts(header_split_count_1, header_split_count_2, refs):
    """
    For each header in reference that has > 0 alignments in 1 file
    
    dict --> {header: (split counts 1:split counts 2)}
    """
    #construct x,y counts for each header
    counts_by_ref = {}
    for header in refs[0]:
        if header in header_split_count_1 and header in header_split_count_2:
            counts_by_ref[header] = (header_split_count_1[header], 
                                     header_split_count_2[header])
        elif header in header_split_count_1 and header \
        not in header_split_count_2:
            counts_by_ref[header] = (header_split_count_1[header], 0)
        elif header not in header_split_count_1 and header in \
        header_split_count_2:
            counts_by_ref[header] = (0, header_split_count_2[header])
    return counts_by_ref 
