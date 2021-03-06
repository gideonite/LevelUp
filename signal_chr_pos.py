#!/usr/bin/python

""" Usage:
signal_chr_pos_toCBS.py   --m=marker_file ...
                    --p=probe_file ...
                    [--O=output_file ]
                    [ (--gistic=gistic-options) ]
                    [--help | -h]

Options:
  --O=output_file                   [default: cbs.out]
  --gistic=gistic-options           make sure to include single quotes around your gistic
                                    options.  The output of cbs is automatically added to these
                                    options.

"""

# todo: should I give users the choice to name the sample?
# [--sample_name=<sample_name>]

# todo: improve the name of this output file

# for probe mapping files, go here:
# http://www.broadinstitute.org/igv/book/export/html/36
# agilentCgh1x1m.txt

from docopt import docopt
import os
import sys
import csv

def marker_position_hash(markerPos_files):
    # turns a marker file [probe    chromosome   position]
    # into a hash table in memory { probe : [chr#, locus] }

    hash = {}

    files = [open(f) for f in markerPos_files]

    for file in files:
        # open the file and read in lines from it
        f = file.readlines()

        # parse out chr, pos, and make a hash to them
        for line in f:
            line = line.split("\t")
            mark = str(line[0]).strip()
            chr = line[1]
            pos = line[2].replace('\n','')

            hash[mark] = [chr, pos]

    file.close()
    return hash


def signal_chr_pos(probes_f, hash):
# go through the signal data,
# match markers to chr loci,
# match with the corresponding signal level,
# and return the name of the sample
# and a list of [signal, chr#, locus]

    list, name, unmapped = [], '', []
    probes_f_open = open(probes_f)
    probe_csv = csv.reader(probes_f_open, delimiter='\t')
    for line in probe_csv:
        if (line[0] == 'Hybridization REF'):
            # todo : return this name somehow
            name = line[1]
            print "<" + name + ">"
        elif (line[0] == 'CompositeElement REF' or line[0] == 'Composite Element REF'):
            assert(line[1] == "normalizedLog2Ratio" or line[1] == "Signal")
            # debug:
            # print f[0], f[1], f[2]
        else:
            mark, signal = line[0], line[1]
            try:
                chr_loc = hash[mark]
                list.append([signal] + chr_loc)
            except KeyError:
                unmapped.append("<" + mark + ">")

    probes_f_open.close()

    if (len(unmapped) != 0):
        if len(unmapped) < 10:
            print "The following probes appear to be unmapped in the marker files: " + " ".join(unmapped)
        print "Total unmapped probes: ", len(unmapped)

        while True:
            yn = raw_input("Would you like to continue anyway?(y/n) ")
            if (yn == 'n'):
                sys.exit(1)
            elif (yn != 'y'):
                yn = raw_input("Would you like to continue anyway?(y/n) ")
            else:
                break

    return (name, list)

if __name__ == '__main__':

    args = docopt(__doc__)

    #print args
    #boom

    #marker_files = ['test_data/marker.na30.lst', 'test_data/marker.na31.lst']
    marker_files = args['--m']
    #probe_f = 'test_data/SN6.CN.test.data.txt'
    probe_files = args['--p']

    SCP_OUT = 'signal_chr_pos.out.tmp'
    ## todo : this is bad.  factor it out

    print "...loading marker files..."
    hash = marker_position_hash(marker_files)
    print "done!"

    # -- CBS -- #
    for probe_f in probe_files:

        sys.stdout.write("...mapping probe signals to chr positions for: ")     # hack, no newline character at end

        s_c_p = signal_chr_pos(probe_f, hash)
        sample_name, scp_list = s_c_p[0], s_c_p[1]
        print "done!"

        # write this to a file to pass to R
        s_c_p_out = open(SCP_OUT, 'w')
        s_c_p_out.write('signal\tchr\tpos\n');    # write the column names (a.k.a header)
        for row in scp_list:
            s_c_p_out.write("%s\n" % '\t'.join(row))
            # write the row to the file
            # tab-deliminited
            # rows separated by new line

        s_c_p_out.close()

        # run cbs
        cbs_cmd = [SCP_OUT, args['--O'], sample_name]      # args
        cbs_cmd = 'Rscript cbs.r ' + ' '.join(cbs_cmd)
        # debug : don't run cbs
        print cbs_cmd
        os.system(cbs_cmd)

        # debug : don't remove this file
        os.system('rm ' + SCP_OUT)

    # -- GISTIC -- #

    if args['--gistic'] != None:
        cmd = './lib/gp_gistic2_from_seg' \
                + ' ' + args['--gistic'] \
                + ' ' + '-seg ' + args['--O']
        os.system(cmd)
