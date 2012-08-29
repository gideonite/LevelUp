library(DNAcopy);      # DNAcopy library available at
                       # http://www.bioconductor.org/packages/2.10/bioc/html/DNAcopy.html
data(coriell);

CNA.object <- CNA(cbind(coriell$Coriell.05296),
    coriell$Chromosome,coriell$Position, data.type="logratio",sampleid="c05296");

CNA.object <- CNA(cbind(coriell$Coriell.05296),
    coriell$Chromosome,coriell$Position, data.type="logratio",sampleid="c05296");

smoothed.CNA.object <- smooth.CNA(CNA.object);

segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1);
