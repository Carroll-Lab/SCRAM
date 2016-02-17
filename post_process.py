'''
Created on 22 Jan 2016

@author: steve
'''
"""
Alignment post-processing module 
"""

import numpy

def fill_in_zeros(fwd_rvs_align_list, ref_len, nt):
    """
    Takes sorted alignment lists (fwd and rvs) as an imput
    Produces 3 lists:
        reference_x_axis (every nucleotid index in the ref_seq)
        fwd_alignment_y_axis (+ve y values for each ref_seq index)
        rvs_alignment_y_axis (-ve y values for each ref_seq index)

    """
    sorted_fwd_alignment=fwd_rvs_align_list[0]
    sorted_rvs_alignment=fwd_rvs_align_list[1]

    fwd_alignment_y_axis=[0]*(ref_len)
    revs_alignment_y_axis=[0]*(ref_len)
    
    reference_x_axis = range(0,ref_len)


    for i in sorted_fwd_alignment:
        fwd_alignment_y_axis[i[0]+(nt/2)]=i[1]
    for i in sorted_rvs_alignment:
        revs_alignment_y_axis[i[0]-(nt/2)]=i[1]


    return  reference_x_axis, fwd_alignment_y_axis, revs_alignment_y_axis


def fill_in_zeros_no_offset(fwd_rvs_align_list, ref_len, nt):
    """
    Takes sorted alignment lists (fwd and rvs) as an imput
    Produces 3 lists:
        reference_x_axis (every nucleotid index in the ref_seq)
        fwd_alignment_y_axis (+ve y values for each ref_seq index)
        rvs_alignment_y_axis (-ve y values for each ref_seq index)

    """
    sorted_fwd_alignment=fwd_rvs_align_list[0]
    sorted_rvs_alignment=fwd_rvs_align_list[1]

    fwd_alignment_y_axis=[0]*(ref_len)
    revs_alignment_y_axis=[0]*(ref_len)
    
    reference_x_axis = range(0,ref_len)


    for i in sorted_fwd_alignment:
        fwd_alignment_y_axis[i[0]]=i[1]
    for i in sorted_rvs_alignment:
        revs_alignment_y_axis[i[0]]=i[1]


    return  reference_x_axis, fwd_alignment_y_axis, revs_alignment_y_axis


def calc_alignments_by_strand(fwd_rvs_align_list):
    """
    FILL output
    """

    sorted_fwd_alignment=fwd_rvs_align_list[0]
    sorted_rvs_alignment=fwd_rvs_align_list[1]
    fwd_align_count = 0
    rvs_align_count = 0
    for i in sorted_fwd_alignment:
        fwd_align_count += i[1]
    for i in sorted_rvs_alignment:
        rvs_align_count -= i[1]  
    return fwd_align_count, rvs_align_count


def calc_5_prime_term_nuc_percentage(alignment_dict):
    """
    Take alignment dict for a singe reference as an input
    (alignment_dict --> read:[[pos,count]])
    
    print percentage of A,G,C,T (U) at 5' end of aligned sRNA 
    """
    a_count =0
    c_count =0
    g_count =0
    u_count =0


    for sRNA, alignment in alignment_dict.iteritems():
        if sRNA[0] == 'A': a_count+= abs(alignment[0][1])
        elif sRNA[0] == 'C': c_count+= abs(alignment[0][1])
        elif sRNA[0] == 'G': g_count+= abs(alignment[0][1])
        elif sRNA[0] == 'T': u_count+=  abs(alignment[0][1])
    total_count = a_count+c_count+g_count+u_count
    print "Percentage of 5' terminal nuc of aligned sRNAs is:"
    print "A: {0}\tC:{1}\tG:{2}\tC:{3}".format(a_count*100/total_count, 
                                               c_count*100/total_count, 
                                               g_count*100/total_count, 
                                               u_count*100/total_count)


def filter_degraded_peaks(graphed_alignments):
    """Use 25-29 sRNA peaks to filter 21,22 and 24nt sRNA peaks at the same location
    """
    sRNAs = graphed_alignments[:3]
    large_sRNAs = graphed_alignments[3:]
    offset=[2,3,4,5,6,7,8,9]
    offset_pos=0 #for calculating how many small sRNA pos to make 0
    for i in large_sRNAs:
        offset_pos_offset = 0
        max_pos_count = max(i[0])
        max_neg_count = min(i[1])

        for pos in range(len(i[0])):
            offset_pos_offset = 0
            if i[0][pos] == 0:
                pass
            elif max_pos_count/i[0][pos] < 10:
                count = 0

                while count < offset[offset_pos+offset_pos_offset]:

                    sRNAs[2][0][pos+count] = 0 #24
                    count+=1
                count=0
                offset_pos_offset+=2
                
                while count < offset[offset_pos+offset_pos_offset]:

                    sRNAs[1][0][pos+count] = 0 #22
                    count+=1

                count=0
                offset_pos_offset+=1
                
                while count < offset[offset_pos+offset_pos_offset]:

                    sRNAs[0][0][pos+count] = 0 #21
                    count+=1

        for pos in range(len(i[1])):
            try:
                offset_pos_offset = 0
                if i[1][pos] == 0:
                    pass
                elif max_neg_count/i[1][pos] < 10:
                    count = 0

                    while count < offset[offset_pos+offset_pos_offset]:

                        sRNAs[2][1][pos+count] = 0 #24
                        count+=1
                    count=0
                    offset_pos_offset+=2
                    
                    while count < offset[offset_pos+offset_pos_offset]:

                        sRNAs[1][1][pos+count] = 0 #22
                        count+=1

                    count=0
                    offset_pos_offset+=1
                    
                    while count < offset[offset_pos+offset_pos_offset]:

                        sRNAs[0][1][pos+count] = 0 #21
                        count+=1
            except:
                pass
        offset_pos+=1

    return sRNAs
    


def smooth(x,window_len,window='hamming'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2-1):-(window_len/2)]

"""End scipy cookbook code
"""

