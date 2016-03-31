'''
Created on 31 Mar 2016

@author: steve
'''

import ref_dict
import time
import align
import analysis_helper as ah
import numpy
import write_to_file as wtf
import post_process as pp
import plot_reads as pr


def ref_coverage(seq, seq_output, ref_file, nt, smoothWinSize, fileFig, 
                 fileName, min_read_size, max_read_size, min_read_no, 
                 onscreen, no_csv, ylim, pub):
    
    """
    Fill out
    """
    
    ref = ref_dict.load_ref_file(ref_file)
    
    if len(ref[0])>1:
        print "\nMutliple reference sequences in file.  \
        1st seq. used for alignment"

    ref_output = ah.single_file_output(ref_file)
    
    
    single_ref = ref[0][ref[1]]
    start = time.clock()
    single_alignment = align.align_reads_to_seq(seq, single_ref, nt)
    if no_csv:
        wtf.csv_output(single_alignment,
                                 nt,
                                 seq_output,
                                 ref_output)   
    if fileFig or onscreen:
        single_sorted_alignemts = align.aln_by_ref_pos(single_alignment)
        graph_processed = pp.fill_in_zeros(single_sorted_alignemts, 
            len(ref[0][ref[1]]), nt)
        x_label = ref[1][1:]
        x_ref = graph_processed[0]
        y_fwd_smoothed = pp.smooth(numpy.array(graph_processed[1]), 
            smoothWinSize, window='blackman')
        y_rvs_smoothed = pp.smooth(numpy.array(graph_processed[2]), 
            smoothWinSize, window='blackman')
        print "\n{0} nt alignment time time = {1} seconds\n"\
            .format(nt, str((time.clock() - start)))
        
        if fileName == "auto":
            fileName = ah.ref_seq_nt_output(seq_output, ref_output, nt, "pdf")
                
        pr.den_plot(x_ref, y_fwd_smoothed, y_rvs_smoothed, nt, fileFig, 
            fileName, onscreen, x_label, ylim, pub)

def coverage_21_22_24(seq, seq_output, ref_file, smoothWinSize, 
    fileFig, fileName, min_read_size, max_read_size, min_read_no,
    onscreen, no_csv,y_lim, pub):     
    """
    Fill out
    """
    ref = ref_dict.load_ref_file(ref_file)
    if len(ref[0])>1:
        print "\nMutliple reference sequences in file.  1st seq. used for \
        alignment"
 
    ref_output = ah.single_file_output(ref_file)
    
    single_ref = ref[0][ref[1]]
    
    combined_21_22_24(seq, seq_output, ref_output, single_ref, smoothWinSize, 
    fileFig, fileName, min_read_size, max_read_size, min_read_no,
    onscreen, no_csv,y_lim, pub)    
       
    
def combined_21_22_24(seq, seq_output, ref_output, single_ref, smoothWinSize, 
    fileFig, fileName, min_read_size, max_read_size, min_read_no,
    onscreen, no_csv,y_lim, pub):
    
    single_alignment_21 = align.align_reads_to_seq(seq, single_ref, 21)
    single_alignment_22 = align.align_reads_to_seq(seq, single_ref, 22)        
    single_alignment_24 = align.align_reads_to_seq(seq, single_ref, 24)

    print '\n21nt sRNAs:'
    single_sorted_alignemts_21 = align.aln_by_ref_pos(single_alignment_21)
    print '\n22nt sRNAs:'
    single_sorted_alignemts_22 = align.aln_by_ref_pos(single_alignment_22)
    print '\n24nt sRNAs:'
    single_sorted_alignemts_24 = align.aln_by_ref_pos(single_alignment_24)
    if no_csv:
        wtf.mnt_csv_output(single_alignment_21, single_alignment_22,
                                 single_alignment_24,
                                 seq_output,
                                 ref_output) 
    if fileFig or onscreen:
    
        graph_processed_21 = pp.fill_in_zeros(single_sorted_alignemts_21, 
            len(single_ref),21)
        graph_processed_22 = pp.fill_in_zeros(single_sorted_alignemts_22, 
            len(single_ref),22)
        graph_processed_24 = pp.fill_in_zeros(single_sorted_alignemts_24, 
            len(single_ref),24)
    
        x_ref = graph_processed_21[0]
        y_fwd_smoothed_21 = pp.smooth(numpy.array(graph_processed_21[1]), 
            smoothWinSize, window='blackman')
        y_rvs_smoothed_21 = pp.smooth(numpy.array(graph_processed_21[2]), 
            smoothWinSize, window='blackman')
        y_fwd_smoothed_22 = pp.smooth(numpy.array(graph_processed_22[1]), 
            smoothWinSize, window='blackman')
        y_rvs_smoothed_22 = pp.smooth(numpy.array(graph_processed_22[2]), 
            smoothWinSize, window='blackman')
        y_fwd_smoothed_24 = pp.smooth(numpy.array(graph_processed_24[1]), 
            smoothWinSize, window='blackman')
        y_rvs_smoothed_24 = pp.smooth(numpy.array(graph_processed_24[2]), 
            smoothWinSize, window='blackman')
    
        if fileName == "auto":
            fileName = ah.ref_seq_output(seq_output, ref_output, "pdf")
    
        pr.den_multi_plot_3(x_ref, y_fwd_smoothed_21, y_rvs_smoothed_21,
        y_fwd_smoothed_22, y_rvs_smoothed_22, y_fwd_smoothed_24, 
        y_rvs_smoothed_24, fileFig, fileName, onscreen, ref_output, y_lim, pub)


    