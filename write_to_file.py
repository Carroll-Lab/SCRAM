'''
Created on 22 Jan 2016

@author: steve
'''
"""
Write results to file
"""

import csv



def circos_output(complete_align):
    """
    Takes dm alignemtn data and produces an output 
    file in a format that can be used for circos plotting

    complete_align format --> (seq_file_name, header 
    (chr), array[21+], array[21-], array[22+], array[22-], 
    array[24+], array[24-])
    """
    seq_file_list = []
    for alignment in complete_align:
        if alignment[0] not in seq_file_list:
            seq_file_list.append(alignment[0])

    for seq_file in seq_file_list:

        f_21=[]
        r_21=[]
        f_22=[]
        r_22=[]
        f_24=[]
        r_24=[]
        for alignment in complete_align:
            if alignment[0] == seq_file:
                chromosome = alignment[1]
                for pos in range(len(alignment[2])):
                    f_21.append((chromosome, pos+1, pos+1, alignment[2][pos]))
                    r_21.append((chromosome, pos+1, pos+1, alignment[3][pos]))
                    f_22.append((chromosome, pos+1, pos+1, alignment[4][pos]))
                    r_22.append((chromosome, pos+1, pos+1, alignment[5][pos]))
                    f_24.append((chromosome, pos+1, pos+1, alignment[6][pos]))
                    r_24.append((chromosome, pos+1, pos+1, alignment[7][pos]))
                    pos+=1
            
        with open(seq_file+'_f_21.txt', 'wb') as csvfile:
            mycsv = csv.writer(csvfile, delimiter='\t')
            for row in f_21:
                out = [row[0], row[1], row[2], round(row[3]*100, 2)]
                mycsv.writerow(out)
        csvfile.close()
        with open(seq_file+'_r_21.txt', 'wb') as csvfile:
            mycsv = csv.writer(csvfile, delimiter='\t')
            for row in r_21:
                out = [row[0], row[1], row[2], round(row[3]*100, 2)]
                mycsv.writerow(out)
        csvfile.close()
        with open(seq_file+'_f_22.txt', 'wb') as csvfile:
            mycsv = csv.writer(csvfile, delimiter='\t')
            for row in f_22:
                out = [row[0], row[1], row[2], round(row[3]*100, 2)]
                mycsv.writerow(out)
        csvfile.close()
        with open(seq_file+'_r_22.txt', 'wb') as csvfile:
            mycsv = csv.writer(csvfile, delimiter='\t')
            for row in r_22:
                out = [row[0], row[1], row[2], round(row[3]*100, 2)]
                mycsv.writerow(out)
        csvfile.close()
        with open(seq_file+'_f_24.txt', 'wb') as csvfile:
            mycsv = csv.writer(csvfile, delimiter='\t')
            for row in f_24:
                out = [row[0], row[1], row[2], round(row[3]*100, 2)]
                mycsv.writerow(out)
        csvfile.close()
        with open(seq_file+'_r_24.txt', 'wb') as csvfile:
            mycsv = csv.writer(csvfile, delimiter='\t')
            for row in r_24:
                out = [row[0], row[1], row[2], round(row[3]*100, 2)]
                mycsv.writerow(out)
        csvfile.close()


def csv_output(alignment_dict, nt, seq_file, header):
    alignment_list=[]
    seq_file_name = seq_file.split('/')[-1]
    for sRNA, alignment in alignment_dict.iteritems():
        for i in alignment:
            alignment_list.append((sRNA,i[0],i[1]))
    alignment_list.sort(key=lambda tup : tup[1])
    with open(header[1:]+'_'+seq_file_name+'_'+str(nt)+'.csv', 'wb') as csvfile:
        mycsv = csv.writer(csvfile, delimiter=',')

        for alignment in alignment_list:
            out=[alignment[0],alignment[1],alignment[2]]
            mycsv.writerow(out)

    csvfile.close()


def cdp_output(counts_by_ref, header1, header2, nt):
    out_file = "{0}_{1}_{2}nt.csv".format(header1, header2, nt)
    results_list=[]
    for header, counts in counts_by_ref.iteritems():
        results_list.append((header, counts[0], counts[1]))
    with open(out_file, 'wb') as csvfile:
        mycsv = csv.writer(csvfile, delimiter=',')
        mycsv.writerow(['',header1, header2])
        for result in results_list:
            out=[result[0],result[1], result[2]]
            mycsv.writerow(out)

    csvfile.close()

def cdp_single_output(counts_by_ref, header1, nt):
    out_file = "{0}_{1}nt.csv".format(header1, nt)
    results_list=[]
    for header, counts in counts_by_ref.iteritems():
        results_list.append((header, counts))
    with open(out_file, 'wb') as csvfile:
        mycsv = csv.writer(csvfile, delimiter=',')
        mycsv.writerow(['',header1])
        for result in results_list:
            out=[result[0],result[1]]
            mycsv.writerow(out)

    csvfile.close()

def csv_filtered_output(filtered_list, alignments, seq_file, header):
    alignment_dict_21_F={} #pos = key
    alignment_dict_21_R={}
    alignment_dict_22_F={} #pos = key
    alignment_dict_22_R={}
    alignment_dict_24_F={} #pos = key
    alignment_dict_24_R={}        
    for sRNA, alignment in alignments[0].iteritems():
        for i in alignment: #as can be more than one alignment
            if i[1]>0: #fwd_Strand_alignment
                alignment_dict_21_F[i[0]] = sRNA
            elif i[1]<0:
                alignment_dict_21_R[i[0]] = sRNA
    for sRNA, alignment in alignments[1].iteritems():
        for i in alignment: #as can be more than one alignment
            if i[1]>0: #fwd_Strand_alignment
                alignment_dict_22_F[i[0]] = sRNA
            elif i[1]<0:
                alignment_dict_22_R[i[0]] = sRNA
    for sRNA, alignment in alignments[2].iteritems():
        for i in alignment: #as can be more than one alignment
            if i[1]>0: #fwd_Strand_alignment
                alignment_dict_24_F[i[0]] = sRNA
            elif i[1]<0:
                alignment_dict_24_R[i[0]] = sRNA


    list_21_row=[]
    list_22_row=[]
    list_24_row=[]
    seq_file_name = seq_file.split('/')[-1]
    for i in range(len(filtered_list[0][0])):
        if i in alignment_dict_21_F and i not in alignment_dict_21_R:
            list_21_row.append((i,filtered_list[0][0][i], 
                                alignment_dict_21_F[i], 
                                filtered_list[0][1][i], 'NULL'))
        elif i in alignment_dict_21_R and i not in alignment_dict_21_F:
            list_21_row.append((i,filtered_list[0][0][i], 'NULL', 
                                filtered_list[0][1][i], alignment_dict_21_R[i]))
        elif i in alignment_dict_21_R and i in alignment_dict_21_F:
            list_21_row.append((i,filtered_list[0][0][i], 
                                alignment_dict_21_F[i], 
                                filtered_list[0][1][i], alignment_dict_21_R[i]))        
        if i in alignment_dict_22_F and i not in alignment_dict_22_R:
            list_22_row.append((i,filtered_list[1][0][i], 
                                alignment_dict_22_F[i], 
                                filtered_list[1][1][i], 'NULL'))
        elif i in alignment_dict_22_R and i not in alignment_dict_22_F:
            list_22_row.append((i,filtered_list[1][0][i], 
                                'NULL', filtered_list[1][1][i], 
                                alignment_dict_22_R[i]))
        elif i in alignment_dict_22_R and i in alignment_dict_22_F:
            list_22_row.append((i,filtered_list[1][0][i], 
                                alignment_dict_22_F[i], filtered_list[1][1][i], 
                                alignment_dict_22_R[i]))    
        if i in alignment_dict_24_F and i not in alignment_dict_24_R:
            list_24_row.append((i,filtered_list[2][0][i], 
                                alignment_dict_24_F[i], 
                                filtered_list[2][1][i], 'NULL'))
        elif i in alignment_dict_24_R and i not in alignment_dict_24_F:
            list_24_row.append((i,filtered_list[2][0][i], 
                                'NULL', filtered_list[2][1][i], 
                                alignment_dict_24_R[i]))
        elif i in alignment_dict_24_R and i in alignment_dict_24_F:
            list_24_row.append((i,filtered_list[2][0][i], 
                                alignment_dict_24_F[i], 
                                filtered_list[2][1][i], alignment_dict_24_R[i]))    
    with open(header[1:]+'_'+seq_file_name+'_21.csv', 'wb') as csvfile:
        mycsv = csv.writer(csvfile, delimiter=',')

        for i in list_21_row:
            out=[i[0],i[1],i[2],i[3],i[4]]
            mycsv.writerow(out)
    csvfile.close()

    with open(header[1:]+'_'+seq_file_name+'_22.csv', 'wb') as csvfile:
        mycsv = csv.writer(csvfile, delimiter=',')

        for i in list_22_row:
            out=[i[0],i[1],i[2],i[3],i[4]]
            mycsv.writerow(out)
    csvfile.close()

    with open(header[1:]+'_'+seq_file_name+'_24.csv', 'wb') as csvfile:
        mycsv = csv.writer(csvfile, delimiter=',')

        for i in list_24_row:
            out=[i[0],i[1],i[2],i[3],i[4]]
            mycsv.writerow(out)
    csvfile.close()

def tospo_circos_output(complete_align, nt, out_file):
    """
    Takes dm alignemtn data and produces an output file 
    in a format that can be used for circos plotting

    complete_align input format --> (([x_ax],[y+],[y-])) L, M, S order

    Output format --> RNA_L    1    1    10.85
    """
    chr_list = ['RNA_L','RNA_M','RNA_S']

    with open(out_file+'_f_'+str(nt)+'.txt', 'wb') as csvfile:
        mycsv = csv.writer(csvfile, delimiter='\t')
        count = 0
        for row in complete_align[0][0]:
            out = [chr_list[0], count, count, round(row*100, 2)]
            mycsv.writerow(out)
            count +=1
        count = 0
        for row in complete_align[1][0]:
            out = [chr_list[1], count, count, round(row*100, 2)]
            mycsv.writerow(out)
            count +=1            
        count = 0
        for row in complete_align[2][0]:
            out = [chr_list[2], count, count, round(row*100, 2)]
            mycsv.writerow(out)
            count +=1

    with open(out_file+'_r_'+str(nt)+'.txt', 'wb') as csvfile:
        mycsv = csv.writer(csvfile, delimiter='\t')
        count = 0
        for row in complete_align[0][1]:
            out = [chr_list[0], count, count, round(row*100, 2)]
            mycsv.writerow(out)
            count +=1
        count = 0
        for row in complete_align[1][1]:
            out = [chr_list[1], count, count, round(row*100, 2)]
            mycsv.writerow(out)
            count +=1            
        count = 0
        for row in complete_align[2][1]:
            out = [chr_list[2], count, count, round(row*100, 2)]
            mycsv.writerow(out)
            count +=1



