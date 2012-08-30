#todo :
args <- commandArgs(TRUE);      # grab trailing args only

# import data, headers [signal chr pos]
input <- read.table("signal_chr_pos.out.tmp", header=TRUE);

#todo: change this
sample_name <- 'test';
output_f <- 'cbs.out';

library(DNAcopy);
## DNAcopy library available at
## http://www.bioconductor.org/packages/2.10/bioc/html/DNAcopy.html#

# run CBS algorithm
CNA.object <- CNA(cbind(input$signal),
    input$chr,input$pos, data.type="logratio",sampleid=sample_name);

smoothed.CNA.object <- smooth.CNA(CNA.object);

segmented <- segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1);

# get the output out of the segmented CNA object
# and write to file
write.csv(segmented$output, output_f);
