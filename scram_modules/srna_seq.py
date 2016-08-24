'''
Created on 1 Apr 2016

@author: steve
'''

from dna import DNA
import time

class SRNA_Seq(object):
    def __init__(self):
        self._internal_dict = {}
        
    def __setitem__(self, sequence, count):
        self._internal_dict[sequence]=count
    
    def __getitem__(self, sequence):
        return self._internal_dict[sequence]
    
    def __iter__(self):
        return self._internal_dict.iteritems()
    
    def __len__(self):
        return len(self._internal_dict)

    def __contains__(self, sequence):
        return sequence in self._internal_dict

    def sRNAs(self):
        return self._internal_dict.keys()

    def counts(self):
        return self._internal_dict.values()

    def load_seq_file(self, seq_file, sRNA_max_len_cutoff, min_reads, \
        sRNA_min_len_cutoff):
        """
        load 1 seq dict in .fasta format.
    
        Calculate RPMR and apply in the function
    
        sRNA_len_cutoff --> only reads  of length <= sRNA_len_cutoff
        accepted (including for RPMR calculation)
        min_reads --> only reads >= min_reads accepted
    
        produce seq_dict --> sRNA:RPMR
        """
        start = time.clock()

        read_count_1 = 0
        loaded_seq_1 = open(seq_file, 'rU')
    
        # Complete dict for seq_1
        for line in loaded_seq_1:
            if line[0] == '>':
                count = int(line.strip().split('-')[1])
                next_line = True #ensure that the nextline has the correct seq.
            elif count >= min_reads and len(line.strip()) <= sRNA_max_len_cutoff \
            and len(line.strip()) >= sRNA_min_len_cutoff and next_line == True:
                self._internal_dict[DNA(line.strip())] = count
                read_count_1 += count
                next_line = False
            else:
                pass
    
        loaded_seq_1.close()

        # final RPMR - could simplify in future
        for sRNA, count in self._internal_dict.iteritems():
            self._internal_dict[sRNA] = count * (float(1000000) / read_count_1)
        print "\nSequence file loading time = "\
         + str((time.clock() - start)) + " seconds\n"
        print "{0} has {1} loaded reads\n".format(seq_file.split('/')[-1],
                                                  read_count_1)
        print "-"*50

    def load_seq_file_arg_list(self, seq_file_arg_list, sRNA_max_len_cutoff, min_reads, \
        sRNA_min_len_cutoff):
        """
        load seq dict from arg_list in fasta format.
    
        Calculate RPMR and apply in the function
    
        sRNA_len_cutoff --> only reads  of length <= sRNA_len_cutoff
        accepted (including for RPMR calculation)
        min_reads --> only reads >= min_reads accepted
    
        produce seq_dict --> sRNA:RPMR
        """
        start = time.clock()

        # read_count_1 = 0
        indv_seq_dict_list=[] #list of indivdual seq_dics
        indv_seq_dict_list_factor=[] #RPMR for each seq disc

        for seq_file in seq_file_arg_list:
            seq_dict={}
            read_count = 0

            loaded_seq = open(seq_file, 'rU')
        
            # Complete dict for seq_1
            for line in loaded_seq:
                if line[0] == '>':
                    count = int(line.strip().split('-')[1])
                    next_line = True #ensure that the nextline has the correct seq.
                elif count >= min_reads and len(line.strip()) <= sRNA_max_len_cutoff \
                and len(line.strip()) >= sRNA_min_len_cutoff and next_line == True:
                    seq_dict[DNA(line.strip())] = count
                    read_count += count
                    next_line = False
                else:
                    pass
        
            loaded_seq.close()

            indv_seq_dict_list.append(seq_dict)
            indv_seq_dict_list_factor.append(float(1000000)/read_count)

        for sRNA, count in indv_seq_dict_list[0].iteritems():
            if all(sRNA in d for d in indv_seq_dict_list):
                total_count = 0
                for i in range(len(indv_seq_dict_list)):
                    total_count+=(indv_seq_dict_list[i][sRNA]*indv_seq_dict_list_factor[i])


                self._internal_dict[sRNA] = total_count/len(indv_seq_dict_list)

        print "\nSequence file loading time = "\
         + str((time.clock() - start)) + " seconds\n"
        print "-"*50