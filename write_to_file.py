'''
Created on 22 Jan 2016

@author: steve
'''
"""
Write results to file
"""

import csv

def csv_output(alignment_dict, nt, seq_file, header):
    """
    Write to file --> csv
    sRNA,pos,count
    """
    alignment_list=[]
    seq_file_name = seq_file.split('/')[-1].split('.')[0]
    for sRNA, alignment in alignment_dict.iteritems():
        for i in alignment:
            alignment_list.append((sRNA,i[0],i[1]))
    alignment_list.sort(key=lambda tup : tup[1])
    with open(header.split('.')[0]+'_'+seq_file_name+'_'+str(nt)\
              +'.csv', 'wb') as csvfile:
        mycsv = csv.writer(csvfile, delimiter=',')
        mycsv.writerow(['sRNA','Position', 'Count'])
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


