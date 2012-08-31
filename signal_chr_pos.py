#!/usr/bin/python

""" Usage:
signal_chr_pos_toCBS.py   --m=marker_file ...
                    --p=probe_file
                    [--O=output_file ]
                    [--sample_name=<sample_name>]

Options:
  --O=output_file   [default: cbs.out]

"""

from docopt import docopt
import os

def marker_position_hash(markerPos_files):
    # turns a marker file [probe    chromosome   position]
    # into a hash table in memory { probe : [chr, locus] }

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


def signal_chr_pos(probesFile_name, hash):
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

if __name__ == '__main__':

    args = docopt(__doc__)

    #marker_files = ['test_data/marker.na30.lst', 'test_data/marker.na31.lst']
    marker_files = args['--m']
    #probe_f = 'test_data/SN6.CN.test.data.txt'
    probe_f = args['--p']

    SCP_OUT = 'signal_chr_pos.out.tmp'
    ## todo : this is bad.  factor it out

    print "...loading marker files..."
    hash = marker_position_hash(marker_files)
    print "done!"

    print "...mapping probes to chromosome positions..."
    s_c_p = signal_chr_pos(probe_f, hash)
    print "done!"

    # write this to a file to pass to R
    s_c_p_out = open(SCP_OUT, 'w')
    s_c_p_out.write('signal\tchr\tpos\n');    # write the column names (a.k.a header)
    for row in s_c_p:
        s_c_p_out.write("%s\n" % '\t'.join(row))
        # write the row to the file
        # tab-deliminited
        # rows separated by new line

    # todo : remove this file
    s_c_p_out.close()

    # run cbs
    cbs_cmd = [SCP_OUT, args['--O'], args['--sample_name']]      # args
    cbs_cmd = 'Rscript cbs.r ' + ' '.join(cbs_cmd)
    os.system(cbs_cmd)
