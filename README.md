SCRAM
---
___

SCRAM is lightweight Python package for aligning small RNA reads to one or more reference sequences and producing publication-quality images.

Developed by Stephen Fletcher @ the laboratory of Prof. Bernie Carroll, University of Queensland

___
*Installation:*

Scram is written in Python 3.5.  Install scram and its dependencies with `pip`:

`pip install scram`

Or download and extract the tarball, then run:

`python setup.py install` 

___
*Input File Format:*

Reference File : DNA nucleotides only (AGCT) - FASTA format

Sequence File : Collapsed reads - DNA nucleotides only (AGCT) - FASTA format

Post-processing of FASTQ reads to collapsed FASTA format can be carried out using the [FASTX-Toolkit from the Hannon Lab](http://hannonlab.cshl.edu/fastx_toolkit/). Collapsed reads are unique, and contain the read count in the header.

An example of the required read file format:
    
`head seq1.fa`
    
    >1-607041
    TCGGACCAGGCATCATTCCCC
    >2-202886
    TCGGACCAGGCTTCATACCCC
    >3-71446
    TCCCAAATATAGACAAAGCA
___

*Usage:*

`scram analysis_type reference_file [-h] [-s1 SEQ_FILE_1] [-s2 SEQ_FILE_2] [-s3 SEQ_FILE_3] [-s4 SEQ_FILE_4] [-nt SRNA_LEN] [-f FILE_NAME] [-seq_list SEQ_LIST] [-min_read MIN_READ_SIZE] [-max_read MAX_READ_SIZE] [-min_count MIN_READ_COUNT] [-win SMOOTH_WIN_SIZE] [-ylim YLIM] [-no_display] [-split] [-pub] [-V]`

Analysis types

* **den** : align reads of a single sRNA class (eg. 21 nt) from a single sequence file to a single reference sequence (-s1 and -nt required)
* **mnt3dm** : align 21, 22 and 24 nt reads from a single sequence file to a single reference sequence (-s1 required)
* **CDP** : count aligned reads of a single sRNA class (eg. 21 nt) to multiple reference sequences. Counts for two sequence files are plotted as (x,y) coordinates for each reference (-s1, -s2 and -nt required)

Flags

* **-h** : Help message
* **-s1** : Sequence file/s 1 (if more than 1 file, read count is averaged for each read (if present in all files)
* **-s2** : Sequence file/s 2 (if more than 1 file, read count is averaged for each read (if present in all files)
* **-nt** : sRNA length to analyse
* **-f** : Figure output file name (if not auto-generated)
* **-p** : No of cores (processors) to use (default=4)
* **-min_read** : Minimum length of sRNA reads used for normalisation (default=18)
* **-max_read** : Maximum length of sRNA reads used for normalisation (default=32)
* **-min_count** : Minimum read count for an sRNA to be aligned and used for normalisation (default=1)
* **-win** : Window size for smoothing of den plots (default=50)
* **-ylim** : +/- y-axis limit on plots
* **-no_display** : Do not display plot on screen
* **-no_csv** : Do not generate .csv alignment file
* **-split** : Split CDP read alignment counts based on no. of alignments
* **-pub** : Remove all labels from density maps for publication
* **-bokeh** : For Jupyter notebook inline plotting when scram started using magic run. No figure output to file.
* **-V** : Show program's version number and exit

___

den example:

`scram den ./ref.fa -s1 seq1.fa -nt 24 -win 30 -f fig1.pdf`

mnt3dm Example:

`scram mnt3dm ./ref.fa -s1 seq1.fa -win 20 -ylim 110 -f fig2.pdf`

CDP Example:

`scram CDP ./cDNAs.fa -s1 seq1.fa -s2 seq2.fa -nt 21 -f fig3.pdf -split`

___

###### (c) 2016 - Stephen Fletcher. MIT License