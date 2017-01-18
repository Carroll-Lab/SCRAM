# coding=utf-8
import unittest
import scram_modules.alignedreads as ar
import scram_modules.srnaseq as srna
import scram_modules.refseq as refseq
import scram_modules.dna as dna
import os

_BASE_DIR = os.path.dirname(os.path.abspath(__file__))

class TestARMethods(unittest.TestCase):

    def test_srna_profile_1(self):
        """
        Test a single read aligning a reference once in the sense orientation

        """
        test_seq = self.load_test_read_file()

        single_ref = self.load_test_ref_file("test_ref_1.fa")

        aligned = self.align_reads(single_ref, test_seq)

        test_aligned=ar.AlignedReads()
        test_aligned[dna.DNA("ATGCGTATGGCGATGAGAGTA")]=[[0, 500000.0]]
        self.assertEqual(aligned, test_aligned)