'''
Created on 22 Jan 2016

@author: steve
'''
"""
Analysis module
"""
import ref_dict
import seq_dict
import post_process
import write_to_file
import cdp
import align
import numpy
import plot_reads
import time


def single_ref_coverage(seq_file, ref_file, nt, smoothWinSize=50, 
    fileFig = False, fileName = 'plot.pdf', min_read_size = 18, 
    max_read_size = 32, min_read_no=1, onscreen = False, pub=False):
    """
    TODO:
    """

    ref = ref_dict.load_ref_file(ref_file)
    if len(ref[0])>1:
        print "\nMutliple reference sequences in file.  \
        1st seq. used for alignment"
    seq = seq_dict.load_seq_file(seq_file, 
        max_read_size, min_read_no, min_read_size)
    single_ref = ref[0][ref[1]]
    start = time.clock()
    single_alignment = align.align_reads_to_seq(seq, single_ref, nt)

    single_sorted_alignemts = align.aln_by_ref_pos(single_alignment)
    graph_processed = post_process.fill_in_zeros(single_sorted_alignemts, 
        len(ref[0][ref[1]]), nt)
    x_label = ref[1][1:]
    x_ref = graph_processed[0]
    y_fwd_smoothed = post_process.smooth(numpy.array(graph_processed[1]), 
        smoothWinSize, window='blackman')
    y_rvs_smoothed = post_process.smooth(numpy.array(graph_processed[2]), 
        smoothWinSize, window='blackman')
    print "\n{0} nt lignment time time = {0} seconds\n"\
        .format(nt, str((time.clock() - start)))
    
    write_to_file.csv_output(single_alignment,nt,seq_file,
                             ref_file.split('/')[-1])
    if fileName == "auto":
        seq_name = seq_file.split('/')[-1].split('.')[0]
        fileName = "{0}_{1}nt.pdf".format(seq_name, str(nt))

    plot_reads.den_plot(x_ref, y_fwd_smoothed, y_rvs_smoothed, nt, fileFig, 
        fileName, onscreen, x_label, pub)



def single_ref_coverage_av(seq_file_1, seq_file_2, ref_file, nt, 
    smoothWinSize=50, fileFig=False, fileName = 'plot.pdf', 
    min_read_size = 18, max_read_size = 32, min_read_no=1, 
    onscreen = False, pub=False):
    """
    TODO: fix fig output
    """

    ref = ref_dict.load_ref_file(ref_file)
    if len(ref[0])>1:
        print "\nMutliple reference sequences in file.  1st seq. used\
         for alignment"
    seq = seq_dict.load_av_seq_files(seq_file_1, seq_file_2, 
        max_read_size, min_read_no, min_read_size)
    single_ref = ref[0][ref[1]]
    start = time.clock()
    single_alignment = align.align_reads_to_seq(seq, single_ref, nt)

    single_sorted_alignemts = align.aln_by_ref_pos(single_alignment)
    graph_processed = post_process.fill_in_zeros(single_sorted_alignemts, 
        len(ref[0][ref[1]]), nt)
    x_label = ref[1][1:]
    x_ref = graph_processed[0]
    y_fwd_smoothed = post_process.smooth(numpy.array(graph_processed[1]), 
        smoothWinSize, window='blackman')
    y_rvs_smoothed = post_process.smooth(numpy.array(graph_processed[2]), 
        smoothWinSize, window='blackman')
    print "\n{0} nt alignment time = {0} seconds\n"\
        .format(nt, str((time.clock() - start)))
    plot_reads.den_plot(x_ref, y_fwd_smoothed, y_rvs_smoothed, nt, fileFig, 
        fileName, onscreen, x_label, pub)



def single_ref_coverage_21_22_24(seq_file, ref_file, smoothWinSize=50, 
    fileFig = True, fileName = 'plot.pdf', min_read_size = 18, 
    max_read_size = 32, min_read_no=1, onscreen = True, y_lim=0, pub=False):
    """
    Align reads from a single seq file to a single reference for 21,22 and 24nt
    """

    ref = ref_dict.load_ref_file(ref_file)
    if len(ref[0])>1:
        print "\nMutliple reference sequences in file.  1st seq. used for \
        alignment"
    seq = seq_dict.load_seq_file(seq_file, max_read_size, min_read_no, 
        min_read_size) 
    single_ref = ref[0][ref[1]]
    single_alignment_21 = align.align_reads_to_seq(seq, single_ref, 21)
    single_alignment_22 = align.align_reads_to_seq(seq, single_ref, 22)        
    single_alignment_24 = align.align_reads_to_seq(seq, single_ref, 24)

    print '\n21nt sRNAs:'
    single_sorted_alignemts_21 = align.aln_by_ref_pos(single_alignment_21)
    print '\n22nt sRNAs:'
    single_sorted_alignemts_22 = align.aln_by_ref_pos(single_alignment_22)
    print '\n24nt sRNAs:'
    single_sorted_alignemts_24 = align.aln_by_ref_pos(single_alignment_24)

    graph_processed_21 = post_process.fill_in_zeros(single_sorted_alignemts_21, 
        len(ref[0][ref[1]]),21)
    graph_processed_22 = post_process.fill_in_zeros(single_sorted_alignemts_22, 
        len(ref[0][ref[1]]),22)
    graph_processed_24 = post_process.fill_in_zeros(single_sorted_alignemts_24, 
        len(ref[0][ref[1]]),24)

    x_ref = graph_processed_21[0]
    x_label = ref[1][1:]
    y_fwd_smoothed_21 = post_process.smooth(numpy.array(graph_processed_21[1]), 
        smoothWinSize, window='blackman')
    y_rvs_smoothed_21 = post_process.smooth(numpy.array(graph_processed_21[2]), 
        smoothWinSize, window='blackman')
    y_fwd_smoothed_22 = post_process.smooth(numpy.array(graph_processed_22[1]), 
        smoothWinSize, window='blackman')
    y_rvs_smoothed_22 = post_process.smooth(numpy.array(graph_processed_22[2]), 
        smoothWinSize, window='blackman')
    y_fwd_smoothed_24 = post_process.smooth(numpy.array(graph_processed_24[1]), 
        smoothWinSize, window='blackman')
    y_rvs_smoothed_24 = post_process.smooth(numpy.array(graph_processed_24[2]), 
        smoothWinSize, window='blackman')


    plot_reads.den_multi_plot_3(x_ref, y_fwd_smoothed_21, y_rvs_smoothed_21,
    y_fwd_smoothed_22, y_rvs_smoothed_22, y_fwd_smoothed_24, y_rvs_smoothed_24,
    fileFig, fileName, onscreen, x_label, y_lim, pub) #fix the True

def single_ref_coverage_21_22_24_av(seq_file_1, seq_file_2, ref_file, 
    smoothWinSize=50, fileFig = True, fileName = 'plot.pdf', min_read_size = 18, 
    max_read_size = 32, min_read_no=1, onscreen = True, y_lim=0, pub=False):
    """
    Align reads from a single seq file to a single reference for 21,22 and 24nt
    """

    ref = ref_dict.load_ref_file(ref_file)
    if len(ref[0])>1:
        print "\nMutliple reference sequences in file.  1st seq. used for \
        alignment"
    seq = seq_dict.load_av_seq_files(seq_file_1, seq_file_2, 
        max_read_size, min_read_no, min_read_size) 
    single_ref = ref[0][ref[1]]
    single_alignment_21 = align.align_reads_to_seq(seq, single_ref, 21)
    single_alignment_22 = align.align_reads_to_seq(seq, single_ref, 22)        
    single_alignment_24 = align.align_reads_to_seq(seq, single_ref, 24)

    print '\n21nt sRNAs:'
    single_sorted_alignemts_21 = align.aln_by_ref_pos(single_alignment_21)
    print '\n22nt sRNAs:'
    single_sorted_alignemts_22 = align.aln_by_ref_pos(single_alignment_22)
    print '\n24nt sRNAs:'
    single_sorted_alignemts_24 = align.aln_by_ref_pos(single_alignment_24)

    graph_processed_21 = post_process.fill_in_zeros(single_sorted_alignemts_21, 
        len(ref[0][ref[1]]),21)
    graph_processed_22 = post_process.fill_in_zeros(single_sorted_alignemts_22, 
        len(ref[0][ref[1]]),22)
    graph_processed_24 = post_process.fill_in_zeros(single_sorted_alignemts_24, 
        len(ref[0][ref[1]]),24)

    x_ref = graph_processed_21[0]
    x_label = ref[1][1:]
    y_fwd_smoothed_21 = post_process.smooth(numpy.array(graph_processed_21[1]), 
        smoothWinSize, window='blackman')
    y_rvs_smoothed_21 = post_process.smooth(numpy.array(graph_processed_21[2]), 
        smoothWinSize, window='blackman')
    y_fwd_smoothed_22 = post_process.smooth(numpy.array(graph_processed_22[1]), 
        smoothWinSize, window='blackman')
    y_rvs_smoothed_22 = post_process.smooth(numpy.array(graph_processed_22[2]), 
        smoothWinSize, window='blackman')
    y_fwd_smoothed_24 = post_process.smooth(numpy.array(graph_processed_24[1]), 
        smoothWinSize, window='blackman')
    y_rvs_smoothed_24 = post_process.smooth(numpy.array(graph_processed_24[2]), 
        smoothWinSize, window='blackman')


    plot_reads.den_multi_plot_3(x_ref, y_fwd_smoothed_21, y_rvs_smoothed_21,
    y_fwd_smoothed_22, y_rvs_smoothed_22, y_fwd_smoothed_24, y_rvs_smoothed_24,
    fileFig, fileName, onscreen, x_label, y_lim, pub) #fix the True

def multi_seq_and_ref_21_22_24(seq_list, ref_file, smoothWinSize=50, 
    fileFig = True, fileName = 'plot.pdf', min_read_size = 18, 
    max_read_size = 32, min_read_no=1, onscreen = False, circos = False, 
    y_lim = 0, pub=False):
    """
    Complete for mutliple seq and ref files 
    See single_ref_coerage_21_22_24 for default values
    pairwise alignments for all seqs and refs will take place
    """
    aln_counts = [] 
    circos_results = []
    ref_seq_count = 0
    # Load all refs --> need to be in a single file I think
    refs = ref_dict.load_ref_file(ref_file)


    # Load all seqs - start with full/rel path and not just names

    seqs = seq_dict.load_seq_list(seq_list) 
    
    for single_seq in seqs:

        seq = seq_dict.load_seq_file(single_seq, max_read_size, min_read_no, 
            min_read_size) 
        ref_seq_count = 0
        for header, single_ref in refs[0].iteritems():
            ref_seq_count +=1
            single_alignment_21 = align.align_reads_to_seq(seq, single_ref, 21)
            single_alignment_22 = align.align_reads_to_seq(seq, single_ref, 22)        
            single_alignment_24 = align.align_reads_to_seq(seq, single_ref, 24)

            print '\n21nt sRNAs:'
            single_sorted_alignemts_21\
             = align.aln_by_ref_pos(single_alignment_21)
            print '\n22nt sRNAs:'
            single_sorted_alignemts_22\
             = align.aln_by_ref_pos(single_alignment_22)
            print '\n24nt sRNAs:'
            single_sorted_alignemts_24\
             = align.aln_by_ref_pos(single_alignment_24)

            #get strand counts
            FR21 = post_process.calc_alignments_by_strand\
            (single_sorted_alignemts_21)
            FR22 = post_process.calc_alignments_by_strand\
            (single_sorted_alignemts_22)
            FR24 = post_process.calc_alignments_by_strand\
            (single_sorted_alignemts_24)
            
            aln_counts.append((header, single_seq, 
                               single_sorted_alignemts_21[2],\
             FR21[0], FR21[1], single_sorted_alignemts_22[2], FR22[0], 
             FR22[1], single_sorted_alignemts_24[2], FR24[0], FR24[1]))

            graph_processed_21\
             = post_process.fill_in_zeros(single_sorted_alignemts_21, 
                 len(single_ref),21)
            graph_processed_22 \
            = post_process.fill_in_zeros(single_sorted_alignemts_22, 
                len(single_ref),22)
            graph_processed_24 \
            = post_process.fill_in_zeros(single_sorted_alignemts_24, 
                len(single_ref),24)

            x_ref = graph_processed_21[0]

            y_fwd_smoothed_21 = post_process.smooth(numpy.array\
                (graph_processed_21[1]), smoothWinSize, window='blackman')
            y_rvs_smoothed_21 = post_process.smooth(numpy.array\
                (graph_processed_21[2]), smoothWinSize, window='blackman')
            y_fwd_smoothed_22 = post_process.smooth(numpy.array\
                (graph_processed_22[1]), smoothWinSize, window='blackman')
            y_rvs_smoothed_22 = post_process.smooth(numpy.array\
                (graph_processed_22[2]), smoothWinSize, window='blackman')
            y_fwd_smoothed_24 = post_process.smooth(numpy.array\
                (graph_processed_24[1]), smoothWinSize, window='blackman')
            y_rvs_smoothed_24 = post_process.smooth(numpy.array\
                (graph_processed_24[2]), smoothWinSize, window='blackman')

            fileName = header[1:]+'_'+single_seq.split('/')[-1]+'.pdf'
            x_label = header[1:]

            if fileFig:
                plot_reads.den_multi_plot_3(x_ref, y_fwd_smoothed_21, 
                    y_rvs_smoothed_21, y_fwd_smoothed_22, y_rvs_smoothed_22, 
                    y_fwd_smoothed_24, y_rvs_smoothed_24, fileFig, fileName, 
                    onscreen, 
                    x_label, y_lim, pub) 

            if circos:
                circos_results.append((single_seq.split('/')[-1], header[1:], 
                    y_fwd_smoothed_21, y_rvs_smoothed_21, y_fwd_smoothed_22, 
                    y_rvs_smoothed_22, y_fwd_smoothed_24, y_rvs_smoothed_24))

    for i in aln_counts:
        print i[0], i[1], i[2], i[5], i[8]

    if circos:

        write_to_file.circos_output(circos_results)


def av_multi_seq_and_ref_21_22_24(seq_list, ref_file, smoothWinSize=50, 
    fileFig = True, fileName = 'plot.pdf', min_read_size = 18, 
    max_read_size = 32, min_read_no=1, onscreen = False, circos = False, 
    y_lim = 0, pub=False):
    """
    Complete for mutliple seq in replicate and ref files 
    See single_ref_coerage_21_22_24 for default values
    pairwise alignments for all seqs and refs will take place
    """
    aln_counts = [] #for getting alignment counts for each 
                    #sRNA class and each ref seq
    circos_results = []
    ref_seq_count = 0
    # Load all refs --> need to be in a single file I think
    refs = ref_dict.load_ref_file(ref_file)


    # Load all seqs - start with full/rel path and not just names

    seqs = seq_dict.load_av_seq_list(seq_list)
    
    for single_seq in seqs:

        seq = seq_dict.load_av_seq_files(single_seq[0], single_seq[1], 
            max_read_size, min_read_no, min_read_size) 
        ref_seq_count = 0
        for header, single_ref in refs[0].iteritems():
            ref_seq_count +=1
            single_alignment_21 = align.align_reads_to_seq(seq, single_ref, 21)
            single_alignment_22 = align.align_reads_to_seq(seq, single_ref, 22)        
            single_alignment_24 = align.align_reads_to_seq(seq, single_ref, 24)

            print '\n21nt sRNAs:'
            single_sorted_alignemts_21 = \
            align.aln_by_ref_pos(single_alignment_21)
            print '\n22nt sRNAs:'
            single_sorted_alignemts_22 = \
            align.aln_by_ref_pos(single_alignment_22)
            print '\n24nt sRNAs:'
            single_sorted_alignemts_24 = \
            align.aln_by_ref_pos(single_alignment_24)

            aln_counts.append((header, single_seq, 
                               single_sorted_alignemts_21[2],\
             single_sorted_alignemts_22[2], single_sorted_alignemts_24[2]))

            graph_processed_21 = post_process.fill_in_zeros\
            (single_sorted_alignemts_21, len(single_ref),21)
            graph_processed_22 = post_process.fill_in_zeros\
            (single_sorted_alignemts_22, len(single_ref),22)
            graph_processed_24 = post_process.fill_in_zeros\
            (single_sorted_alignemts_24, len(single_ref),24)

            x_ref = graph_processed_21[0]

            y_fwd_smoothed_21 = post_process.smooth(numpy.array\
                (graph_processed_21[1]), smoothWinSize, window='blackman')
            y_rvs_smoothed_21 = post_process.smooth(numpy.array\
                (graph_processed_21[2]), smoothWinSize, window='blackman')
            y_fwd_smoothed_22 = post_process.smooth(numpy.array\
                (graph_processed_22[1]), smoothWinSize, window='blackman')
            y_rvs_smoothed_22 = post_process.smooth(numpy.array\
                (graph_processed_22[2]), smoothWinSize, window='blackman')
            y_fwd_smoothed_24 = post_process.smooth(numpy.array\
                (graph_processed_24[1]), smoothWinSize, window='blackman')
            y_rvs_smoothed_24 = post_process.smooth(numpy.array\
                (graph_processed_24[2]), smoothWinSize, window='blackman')

            fileName = header[1:]+'_'+single_seq[0].split('/')[-1].\
            split('.')[0][:-2]+'.pdf'
            x_label = header[1:]

            if fileFig:
                plot_reads.den_multi_plot_3(x_ref, y_fwd_smoothed_21, 
                    y_rvs_smoothed_21, y_fwd_smoothed_22, y_rvs_smoothed_22, 
                    y_fwd_smoothed_24, y_rvs_smoothed_24, fileFig, 
                    fileName, onscreen, 
                    x_label, y_lim, pub)

            if circos:
                circos_results.append((single_seq.split('/')[-1], header[1:], 
                    y_fwd_smoothed_21, y_rvs_smoothed_21, y_fwd_smoothed_22, 
                    y_rvs_smoothed_22, y_fwd_smoothed_24, y_rvs_smoothed_24))

    if circos:
        # format --> (seq_file_name, header (chr), 21+, 21-, 22+, 22-, 24+, 24-)

        write_to_file.circos_output(circos_results)
    for i in aln_counts:
        print i[0], i[1][0], i[1][1], i[2], i[3], i[4]

def CDP(seq_file_1, seq_file_2, ref_file, nt, 
    fileFig=False, fileName = 'plot.pdf', 
    min_read_size = 18, max_read_size = 32, min_read_no=1, onscreen = False,
    pub=False):
    """
    Plots aligmnent count for each sRNA in ref file as (x,y)
    for 2 seq files.  No splitting of read count
    """  

    refs = ref_dict.load_ref_file(ref_file)
    seq_1 = seq_dict.load_seq_file(seq_file_1, 
        max_read_size, min_read_no, min_read_size)
    seq_2 = seq_dict.load_seq_file(seq_file_2, 
        max_read_size, min_read_no, min_read_size)    

    counts_by_ref = {} #header:(count1, count2)
    
    for header, single_ref in refs[0].iteritems():
        #ref_seq_count +=1
        single_alignment_1 = align.count_align_reads_to_seq(seq_1, 
                                                            single_ref, nt)
        single_alignment_2 = align.count_align_reads_to_seq(seq_2, 
                                                            single_ref, nt)        
        if single_alignment_1 != 0 or single_alignment_2 !=0:
            counts_by_ref[header] = (single_alignment_1, single_alignment_2)
    # for header, alignments in counts_by_ref.iteritems():
    #     print header, alignments
    seq1 = seq_file_1.split('/')[-1][:-5]
    seq2 = seq_file_2.split('/')[-1][:-5]
    if onscreen ==True:
        plot_reads.cdp_plot(counts_by_ref, 
                            seq1, 
                            seq2,
                            nt,
                            onscreen,
                            fileFig,
                            fileName,
                            pub)
    write_to_file.cdp_output(counts_by_ref, 
                             seq_file_1.split("/")[-1].split(".")[-2], 
                             seq_file_2.split("/")[-1].split(".")[-2],nt)


def CDP_split(seq_file_1, seq_file_2, ref_file, nt, 
    fileFig=False, fileName = 'plot.pdf', 
    min_read_size = 18, max_read_size = 32, min_read_no=1, onscreen = False,
    pub=False):

    """
    Plots aligmnent count for each sRNA in ref file as (x,y)
    for 2 seq files.  Read count split by number of times an sRNA aligns
    
    TODO: combine with avCDP and make a new function for the repeated stuff
    
    """  


    refs = ref_dict.load_ref_file(ref_file)
    seq_1 = seq_dict.load_seq_file(seq_file_1, 
        max_read_size, min_read_no, min_read_size)
    seq_2 = seq_dict.load_seq_file(seq_file_2, 
        max_read_size, min_read_no, min_read_size)    

    alignment_dict_1={} #header:aligned_sRNAs
    alignment_dict_2={}
    
    #calc aligned sRNAs for each header, duplicate if necessary
    for header, single_ref in refs[0].iteritems():

        alignment_dict_1[header] = \
        cdp.dict_align_reads_to_seq_split(seq_1, single_ref, nt)
        alignment_dict_2[header] = \
        cdp.dict_align_reads_to_seq_split(seq_2, single_ref, nt)        
    
    times_align_1=cdp.times_read_aligns(alignment_dict_1)
    times_align_2=cdp.times_read_aligns(alignment_dict_2)
    
    header_split_count_1=cdp.split_reads_for_header(alignment_dict_1, 
                                                    times_align_1, 
                                                    seq_1)
    header_split_count_2=cdp.split_reads_for_header(alignment_dict_2, 
                                                    times_align_2, 
                                                    seq_2)  
    
    
    counts_by_ref = cdp.header_x_y_counts(header_split_count_1, 
                                          header_split_count_2, 
                                          refs)  

    seq1 = seq_file_1.split('/')[-1].split('.')[0] #x_axis name
    seq2 = seq_file_2.split('/')[-1].split('.')[0] #y axis name
    if onscreen ==True:
        plot_reads.cdp_plot(counts_by_ref, 
                            seq1, 
                            seq2,
                            nt,
                            onscreen,
                            fileFig,
                            fileName,
                            pub)
    write_to_file.cdp_output(counts_by_ref, 
                             seq_file_1.split("/")[-1].split(".")[-2], 
                             seq_file_2.split("/")[-1].split(".")[-2],nt)



def avCDP(seq_file_1, seq_file_2, seq_file_3, seq_file_4, ref_file, nt, 
    fileFig=False, fileName = 'plot.pdf', 
    min_read_size = 18, max_read_size = 32, min_read_no=1, onscreen = False, 
    pub=False):
    refs = ref_dict.load_ref_file(ref_file)
    seq_1 = seq_dict.load_av_seq_files(seq_file_1, seq_file_2, 
        max_read_size, min_read_no, min_read_size)
    seq_2 = seq_dict.load_av_seq_files(seq_file_3, seq_file_4, 
        max_read_size, min_read_no, min_read_size)    
    counts_by_ref = {} #header:(count1, count2)
    
    for header, single_ref in refs[0].iteritems():
        #ref_seq_count +=1
        single_alignment_1 = align.count_align_reads_to_seq(seq_1, 
                                                            single_ref, nt)
        single_alignment_2 = align.count_align_reads_to_seq(seq_2, 
                                                            single_ref, nt)        
        if single_alignment_1 != 0 or single_alignment_2 !=0:
            counts_by_ref[header] = (single_alignment_1, single_alignment_2)
    seq1 = "{0} + {1}".format(seq_file_1.split('/')[-1].split('.')[0],
                            seq_file_2.split('/')[-1].split('.')[0])

    #y_axis name
    seq2 = "{0} + {1}".format(seq_file_3.split('/')[-1].split('.')[0],
                            seq_file_4.split('/')[-1].split('.')[0])

    if onscreen ==True:
        plot_reads.cdp_plot(counts_by_ref, 
                            seq1, 
                            seq2,
                            nt,
                            onscreen,
                            fileFig,
                            fileName,
                            pub)    
    write_to_file.cdp_output(counts_by_ref, 
                             "{0}_{1}".format(seq_file_1.split('/')[-1]\
                                              .split('.')[0],
                            seq_file_2.split('/')[-1].split('.')[0]), 
                             "{0}_{1}".format(seq_file_3.split('/')[-1]\
                                              .split('.')[0],
                            seq_file_4.split('/')[-1].split('.')[0]),
                             nt)    


def avCDP_split(seq_file_1, seq_file_2, seq_file_3, seq_file_4, ref_file, nt, 
    fileFig=False, fileName = 'plot.pdf', 
    min_read_size = 18, max_read_size = 32, min_read_no=1, onscreen = False, 
    pub=False):
 
    refs = ref_dict.load_ref_file(ref_file)
    seq_1 = seq_dict.load_av_seq_files(seq_file_1, seq_file_2, 
        max_read_size, min_read_no, min_read_size)
    seq_2 = seq_dict.load_av_seq_files(seq_file_3, seq_file_4, 
        max_read_size, min_read_no, min_read_size)    
 
    alignment_dict_1={} #header:aligned_sRNAs
    alignment_dict_2={}
    
    #calc aligned sRNAs for each header, duplicate if necessary
    for header, single_ref in refs[0].iteritems():

        alignment_dict_1[header] = \
        cdp.dict_align_reads_to_seq_split(seq_1, single_ref, nt)
        alignment_dict_2[header] = \
        cdp.dict_align_reads_to_seq_split(seq_2, single_ref, nt)        
    
    times_align_1=cdp.times_read_aligns(alignment_dict_1)
    times_align_2=cdp.times_read_aligns(alignment_dict_2)
    
    header_split_count_1=cdp.split_reads_for_header(alignment_dict_1, 
                                                    times_align_1, 
                                                    seq_1)
    header_split_count_2=cdp.split_reads_for_header(alignment_dict_2, 
                                                    times_align_2, 
                                                    seq_2)  
    
    
    counts_by_ref = cdp.header_x_y_counts(header_split_count_1, 
                                          header_split_count_2, 
                                          refs)  

    #x_axis name
    seq1 = "{0} + {1}".format(seq_file_1.split('/')[-1].split('.')[0],
                            seq_file_2.split('/')[-1].split('.')[0])

    #y_axis name
    seq2 = "{0} + {1}".format(seq_file_3.split('/')[-1].split('.')[0],
                            seq_file_4.split('/')[-1].split('.')[0])
    if onscreen ==True:
        plot_reads.cdp_plot(counts_by_ref, 
                            seq1, 
                            seq2,
                            nt,
                            onscreen,
                            fileFig,
                            fileName,
                            pub)
    write_to_file.cdp_output(counts_by_ref, 
                             "{0}_{1}".format(seq_file_1.split('/')[-1]\
                                              .split('.')[0],
                            seq_file_2.split('/')[-1].split('.')[0]), 
                             "{0}_{1}".format(seq_file_3.split('/')[-1]\
                                              .split('.')[0],
                            seq_file_4.split('/')[-1].split('.')[0]),
                             nt)  


def CDP_single_split(seq_file_1, ref_file, nt, 
    min_read_size = 18, max_read_size = 32, min_read_no=1):
    """
    Normalised split alignments for a single reference - all recorded
    including 0
    
    
    TODO: may be incorrect - check with unit tests!!!!!!
    """

    refs = ref_dict.load_ref_file(ref_file)
    seq_1 = seq_dict.load_seq_file(seq_file_1, 
        max_read_size, min_read_no, min_read_size)


    alignment_dict_1={} #header:aligned_sRNAs


    sRNA_align_count_1={} #sRNA:no_of_times_aligned


    header_split_count_1={} #header:count


    counts_by_ref = {} #header (align_count_1, align_count_2)

    #calc aligned sRNAs for each header, duplicte if necessary
    for header, single_ref in refs[0].iteritems():

        alignment_dict_1[header] = \
        align.list_align_reads_to_seq_split(seq_1, single_ref, nt)
        
    
    #calc no of times each sRNA is aligned - 1    
    for header, sRNA_list in alignment_dict_1.iteritems():
        for sRNA in sRNA_list:
            if sRNA in sRNA_align_count_1:
                sRNA_align_count_1[sRNA] +=1
            else:
                sRNA_align_count_1[sRNA] = 1


    #calc split alignment count for each header - 1
    for header, sRNA_list in alignment_dict_1.iteritems():
        header_split_count_1[header] = 0
        for sRNA in sRNA_list:
            header_split_count_1[header] +=seq_1[sRNA]/sRNA_align_count_1[sRNA]

    #counstruc x,y counts for each header
    for header in refs[0]:

        counts_by_ref[header] = header_split_count_1[header]


    write_to_file.cdp_single_output(counts_by_ref, 
                                    seq_file_1.split("/")[-1].split(".")[-2], 
                                    nt)

 