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

    def test_srna_profile_2(self):
        """
        Test 2 reads aligning to the same reference once - sense only
        :return:
        """
        test_seq = self.load_test_read_file()

        single_ref = self.load_test_ref_file("test_ref_2.fa")

        aligned = self.align_reads(single_ref, test_seq)

        test_aligned=ar.AlignedReads()
        test_aligned[dna.DNA("GTGCGTATGGCGATGAGAGTA")]=[[27, 250000.0]]
        test_aligned[dna.DNA("ATGCGTATGGCGATGAGAGTA")]=[[0, 500000.0]]

        self.assertEqual(aligned, test_aligned)

    def test_srna_profile_3(self):
        """
        Test 2 reads aligning to the same reference once - 1 sense and 1 anti-sense
        :return:
        """
        test_seq = self.load_test_read_file()

        single_ref = self.load_test_ref_file("test_ref_3.fa")

        aligned = self.align_reads(single_ref, test_seq)

        test_aligned = ar.AlignedReads()
        test_aligned[dna.DNA("GTGCGTATGGCGATGAGAGTA")] = [[47, -250000.0]]
        test_aligned[dna.DNA("ATGCGTATGGCGATGAGAGTA")] = [[0, 500000.0]]

        self.assertEqual(aligned, test_aligned)

    def test_srna_profile_4(self):
        """
        Test a single read aligning a reference twice in the sense orientation
        No split of aligned read count

        """
        test_seq = self.load_test_read_file()

        single_ref = self.load_test_ref_file("test_ref_4.fa")

        aligned = self.align_reads(single_ref, test_seq)

        test_aligned=ar.AlignedReads()
        test_aligned[dna.DNA("ATGCGTATGGCGATGAGAGTA")]=[[0, 500000.0],[27,500000.0]]
        self.assertEqual(aligned, test_aligned)

    def test_srna_profile_5(self):
        """
        Test a single read aligning a reference twice in the sense orientation
        Split of aligned read count

        """
        test_seq = self.load_test_read_file()

        single_ref = self.load_test_ref_file("test_ref_4.fa")

        aligned = self.align_reads(single_ref, test_seq)
        aligned.split()
        test_aligned=ar.AlignedReads()
        test_aligned[dna.DNA("ATGCGTATGGCGATGAGAGTA")]=[[0, 250000.0],[27,250000.0]]
        self.assertEqual(aligned, test_aligned)

    def load_test_read_file(self):
        """
        Load test read file
        :return: SRNASeq object
        """
        seq_file = _BASE_DIR + "/test_seq.fa"
        test_seq = srna.SRNASeq()
        test_seq.load_seq_file(seq_file, 50, 1, 0)
        return test_seq


    def load_test_ref_file(self,ref_name):
        """
        Load test_ref_file
        :return: a single reference (DNA)
        """
        ref_file = "{0}/{1}".format(_BASE_DIR,ref_name)
        test_ref = refseq.RefSeq()
        test_ref.load_ref_file(ref_file)
        single_ref = ""
        for header, ref_seq in test_ref:
            single_ref = ref_seq
        return single_ref


    def align_reads(self, single_ref, test_seq):
        """
        ALign Reads
        :param single_ref:
        :param test_seq:
        :return:
        """
        aligned = ar.AlignedReads()
        aligned.align_reads_to_ref(test_seq, single_ref, 21)
        return aligned


if __name__ == '__main__':
    unittest.main()