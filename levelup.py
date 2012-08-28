#!/usr/bin/python

""" Usage:
LevelUp.py  (--marker-file=<marker-file>) ...
            (--probe-file=<probe-file>) ...
            [--sample-name=<sample-name>] ...
            [--cbs-output-file=<cbs-output-file>]
            [ (--gistic-exec=<path> --gistic-options=<options>) ]
            [--help | -h]

Pipeline for getting from MIT Broad's Level_2 data to Level_4 data.  Level_2
data gets converted to Level_3 data by the Circular Binary Segmentation (CBS)
algorithm.  Level_3 data gets converted to Level_4 data by the Genomic
Identification of Significant Targets in Cancer (GISTIC).

Arguments: CBS
  probe-file            The copy number data that you are going to segment, the
                        columns are something like [probe_name signal_strength]
  marker-file           The marker files to be used as reference the columns are
                        something like [probe_name chr position]
Arguments: GISTIC
  gistic-exec           The gistic executable
  gistic-options        Options to get passed directly to the gistic program.
                        Enclose in quotes ' .  For details see the gistic docs

Options:
  -h --help                     Print this message
  --sample-name=NAME            Uses NAME instead of the sample's filename
                                for the CBS algorithm
                                NB: order matters.  The first listed name gets
                                matched to the first probe file, etc.  Probably
                                best to add them together in pairs, i.e.
                                "--probe-file SNP-probe-file --sample-name=SNP"
                                probe-file and its name form a pair
  --cbs-output-file=OUTPUT      Appends the output of CBS to a file named
                                OUTPUT.  [default: cbs.out]
  --gistic-exec=<path>          Path to executable GISTIC program
  --gistic-options=<options>    These commandline options are passed directly
                                to the gistic program.  If omitted, gistic will
                                not be run.
"""

import sys, os
import rpy2.robjects as robjects
import rpy2.robjects.lib.ggplot2 as ggplot2
from rpy2.robjects.packages import importr
from docopt import docopt

# useful links and dependencies:

    # quick and dirty:  import os   os.system(' ... ')
    # http://stackoverflow.com/questions/450285/executing-command-line-programs-from-within-python

    # rpy :
    # http://pypi.python.org/pypi/rpy2/2.2.6

    # CNV file
    # http://www.broadinstitute.org/mpr/publications/projects/CNS/100K_CNVs_080423.txt
    # todo : what is the meaning of this file?  Isn't it redundant to the probesFile?
    # Boris knows

    # provide the following link to updated versions of GISTIC:
    # http://www.broadinstitute.org/cgi-bin/cancer/publications/pub_paper.cgi?mode=view&paper_id=162

class levelup:

    def marker_position_hash(self, markerPos_files):
    # takes a list of marker files (mappings of probe -> chr locus)
    # and turns it into a hash table in memory { mark : [chr, locus] }
        print "...loading marker files..."

        hash = {}
        for file in markerPos_files:

            # open the file and read in lines from it
            file = open(file)
            f = file.readlines()

            # parse out chr, pos, and make a hash to them
            for line in f:
                line = line.split("\t")
                mark = str(line[0]).strip()
                chr = line[1]
                pos = line[2].replace('\n','')
                hash[mark] = [chr, pos]

        file.close()

        print "done!"
        return hash

    def probes_to_chrLocus(self, probesFile_name, hash):
    # go through the signal data,
    # match markers to chr loci,
    # match with the corresponding signal level,
    # and return a list of [signal, chr#, position]

        probesFile = open(probesFile_name)

        f = probesFile.readlines()

        list = []

        if (f[1].find('Signal') == -1):
            print f[0]
            print f[1]
            print f[2]
            print "It appears that your data file does not have a Signal Column"

        # *** ignore the first two lines of the probe file ***
        # they are usually column names and whatnot
        for line in f[2:]:
            line = line.split('\t')

            mark = line[0]
            signal = line[1].replace('\n', '')

            try:
                map = hash[mark]
            except KeyError:
                print "The following probe appears to be unmapped in the marker files: <" + mark + ">"

            map.insert(0, signal)
            list.append(map)

        probesFile.close()

        return list


    def CBS(self, list, name_of_sample):
    # run the CBS algorithm on a *python list* of [signal, chr#, position]
    # returns an R object of segmented data
    # name_of_sample eg. secondary_GBM_6, primary_GBM_30, etc.

    # todo : make the output silent.
    # todo: make data_type an option for the user?
    # Question : does CBS always output logratio?
        data_type = "logratio"          
        # convert data_type to R
        data_type = robjects.r('data.type="' + data_type + '"')

        id = name_of_sample
        # convert id to R
        id = robjects.r('sampleid="' + id + '"')

        # create a CNA object (signal, chr, pos, data.type, sampleid)
        signal = []
        chr = []
        pos = []
        for locus in list:
            signal.append(locus[0])
            chr.append(locus[1])
            pos.append(locus[2])

        if (len(signal) != len(chr) or len(chr) != len(pos) or len(signal) != len(pos)):
            print '''there has to be the same number of signals, chromosomes, and
            positions'''
            sys.exit(-1)

        # convert to R
        signal = robjects.FloatVector(signal)
        chr = robjects.StrVector(chr)       # StrVector because of 'X' and 'Y' chromosomes
        pos = robjects.IntVector(pos)

        # import DNAcopy
        DNAcopy = importr("DNAcopy")

        # create CNA object
        CNA = DNAcopy.CNA(signal, chr, pos, data_type, id)

        # smooth
        smooth = robjects.r['smooth.CNA']
        smoothed = smooth(CNA)

        # segment
        segment = robjects.r['segment']
        verbose = robjects.r('verbose=1')

        segmented = segment(smoothed)
        #segmented = segment(smoothed, verbose)     # todo: why doesn't this work?

        return segmented

    def write_to_file(self, R_obj, filename):
    # takes an R object and writes it to a file
    # e.g. takes the output of CBS and writes it to a file for GISTIC to run

       R_obj = robjects.r('print')(R_obj)
       R_obj = robjects.r('write.table')(R_obj, filename, sep='\t', append='TRUE')
       # tab delimited

    def clean_up_cbs_output(self, filename):
        tmp_file = filename + ".tmp"

        cmd = "sed \'s/\"//g\' " + filename             # replace all " with ''
        cmd += "| cut -f 2- "                           # delete first column of row #s
    #    cmd += "| sed \'s/ /     /g\' | sed \'1d\'"    # make tab-delimited
        cmd += ">>" + tmp_file

        os.system(cmd)
        os.system('mv ' + tmp_file + ' ' + filename)


    if __name__ == "__main__":

        # use docopt to parse our args
        args = docopt(__doc__)
        print args
        sys.exit(0)

        # get out the args and deal with them!
        cbs_output_filename = args['--cbs-output-file']
        markerPos_files = args['--marker-file']
        probe_files = args['--probe-file']
        sample_name = args['--sample-name']

        # set by docopt
        #if not cbs_output_filename:
        #    cbs_output_filename = 'cbs.out'

        if sample_name != None:
            probes_and_names = zip(probe_files, sample_name)
            probe_files = probes_and_names + probe_files[len(probes_and_names):]

        print markerPos_files
        hash = marker_position_hash(markerPos_files)

        for probe_file in probe_files:
        # this is a concatenation of proble_file_name their corresponding
        # sample_name pairs with a list of probe_file_names.
        # [(probe_file_name, sample_name), probe_file_names]
            if len(probe_file) == 2:
                filename = probe_file[0]
                sample_name = probe_file[1]
                print probe_file
            elif type(probe_file) == 'str' and len(probe_file) > 2:
                raise IndexError('''more than two things got zipped together with a
                probe file and its sample name''')
            else:
                filename = probe_file
                sample_name = probe_file
                # todo : some regexs or module to get a proper name

            cbs_input = probes_to_chrLocus(filename, hash)
            cbs = CBS(cbs_input, sample_name)
            write_to_file(cbs, cbs_output_filename)

        clean_up_cbs_output(cbs_output_filename)

        gistic = args['--gistic-exec']
        gistic_args = args['--gistic-options']
        gistic_cmd = ''
        if gistic:
            print '--- running gistic ---'

            if gistic_args == True:
            # gistic_args is the empty string
                os.system(gistic)
                sys.exit(0)
            else:
                gistic_cmd = gistic + ' ' + gistic_args
                print gistic_cmd

        os.system(gistic_cmd)
