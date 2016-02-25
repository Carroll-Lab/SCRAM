'''
Created on 25 Feb 2016

@author: steve
'''

import align

class TestAlignSimple:

    _ref = "AAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGG" 
    
    def test_align_1(self):
        """
        Simple align - sRNA present
        """
        seq_dict={}
        align_dict={}
        
        seq_dict["AAAAAAAAAAAAAAAAAAAAA"] = 10
        align_dict["AAAAAAAAAAAAAAAAAAAAA"]=[[0,10]]
        assert align.align_reads_to_seq(seq_dict, 
                                        TestAlignSimple._ref, 
                                        21) == align_dict
        
        
    def test_align_2(self):
        """
        Simple align - sRNA not present
        """
        seq_dict={}
        align_dict={}
        
        seq_dict["AAAAAAAATAAAAAAAAAAAA"] = 10
        align_dict["AAAAAAAAAAAAAAAAAAAAA"]=[[0,10]]
        assert align.align_reads_to_seq(seq_dict, 
                                        TestAlignSimple._ref, 
                                        21) != align_dict
                                        
    def test_align_3(self):
        """
        Simple align - sRNA not present - boundary
        """
        seq_dict={}
        align_dict={}
        
        seq_dict["TAAAAAAAAAAAAAAAAAAAA"] = 10
        align_dict["AAAAAAAAAAAAAAAAAAAAA"]=[[0,10]]
        assert align.align_reads_to_seq(seq_dict, 
                                        TestAlignSimple._ref, 
                                        21) != align_dict
    def test_align_4(self):
        """
        Simple align - sRNA not present - boundary
        """
        seq_dict={}
        align_dict={}
        
        seq_dict["AAAAAAAAAAAAAAAAAAAAT"] = 10
        align_dict["AAAAAAAAAAAAAAAAAAAAA"]=[[0,10]]
        assert align.align_reads_to_seq(seq_dict, 
                                        TestAlignSimple._ref, 
                                        21) != align_dict
    def test_align_5(self):
        """
        Simple align - sRNA not present - boundary
        """
        seq_dict={}
        align_dict={}
        
        seq_dict["TTTTTTTTTTTTTTTTTTTTT"] = 10
        align_dict["TTTTTTTTTTTTTTTTTTTTT"]=[[0,10]]
        assert align.align_reads_to_seq(seq_dict, 
                                        TestAlignSimple._ref, 
                                        21) != align_dict
                                        
    def test_align_double_1(self):
        """
        Same sRNA aligns twice
        """
        seq_dict={}
        align_dict={}
        seq_dict["GGGGGGGGGGGGGGGGGGGGG"] = 10
        align_dict["GGGGGGGGGGGGGGGGGGGGG"]=[[21,10],[22,10]]
        
        assert align.align_reads_to_seq(seq_dict, 
                                    TestAlignSimple._ref, 
                                    21) == align_dict   
        
        
    def test_rvs_comp_align(self):
        """
        Check position of reverse comp. align
        """
        seq_dict={}
        align_dict={}
        seq_dict["CCCCCCCCCCCCCCCCCCCCC"] = 10
        align_dict["CCCCCCCCCCCCCCCCCCCCC"]=[[21,10],[22,10]]   
        