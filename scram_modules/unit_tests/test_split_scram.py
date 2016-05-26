'''
Created on 25 Feb 2016

@author: steve
'''
import cdp

class TestAlignSplit:

    _refs = {"a":"AAAGGGG",
             "b":"AAAGGGCCC"}
    
    def test_split_align_1(self):
        """
        Same sRNA aligns twice - split reads
        """
        seq_dict_1={}
        seq_dict_1["AAA"] = 10
        seq_dict_1["GGG"] = 16
  
        seq_dict_2={}
        seq_dict_2["AAA"] = 16
        seq_dict_2["CCC"] = 4      

        alignment_dict_1={} #header:aligned_sRNAs
        alignment_dict_2={}

         
        for header, single_ref in TestAlignSplit._refs.iteritems():

            alignment_dict_1[header] = \
                cdp.dict_align_reads_to_seq_split(seq_dict_1, single_ref, 3)
            alignment_dict_2[header] = \
                cdp.dict_align_reads_to_seq_split(seq_dict_2, single_ref, 3)   
        
        assert alignment_dict_1=={"a":{"GGG":2,
                                  "AAA":1},
                                  "b":{"GGG":2,
                                  "AAA":1}}

        assert alignment_dict_2=={"a":{"CCC":2,
                                  "AAA":1},
                                  "b":{"CCC":2,
                                  "AAA":1}}
        
        times_align_1=cdp.times_read_aligns(alignment_dict_1)
        times_align_2=cdp.times_read_aligns(alignment_dict_2)
        
        assert times_align_1=={"AAA":2, "GGG":4}
        assert times_align_2=={"AAA":2, "CCC":4}
        
        header_split_count_1=cdp.split_reads_for_header(alignment_dict_1, 
                                                        times_align_1, 
                                                        seq_dict_1)
        header_split_count_2=cdp.split_reads_for_header(alignment_dict_2, 
                                                        times_align_2, 
                                                        seq_dict_2) 
        assert header_split_count_1 == {"a":13, "b":13}
        assert header_split_count_2 == {"a":10, "b":10}
        