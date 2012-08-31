#todo :
args <- commandArgs(TRUE);      # grab trailing args only

if (length(args) != 3) {
    print("Usage: cbs.r  <inputfile [signal chr pos]> <output_file> <sample_name>");
    q();        # exit (quit)
}

# import data, headers [signal chr pos]
cat("...CBS is importing data...\n");
input <- read.table(args[1], header=TRUE);
output_f <- args[2];
sample_name <- args[3];
cat("done!\n");

library(DNAcopy);
## DNAcopy library available at
## http://www.bioconductor.org/packages/2.10/bioc/html/DNAcopy.html#

# run CBS algorithm
CNA.object <- CNA(as.numeric(input$signal),
        as.numeric(input$chr),
        as.numeric(input$pos),
        data.type="logratio",
        sampleid=sample_name);

cat("...CBS is smoothing the probe level data...\n");
smoothed.CNA.object <- smooth.CNA(CNA.object);
print("done!");

cat("..CBS is segmenting...\n");
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1);
cat("done!\n");

# get the output out of the segmented CNA object
# and write to file
# don't print the row names
# append to end of file
cat("CBS is writing results to: ", output_f, '\n');
write.table(segment.smoothed.CNA.object$output, output_f,
append=TRUE, row.names=FALSE, eol='\n', sep='\t');
