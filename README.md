# CNV-detection-using-read-depths
The script CNV.py detects Copy Number Variations in the genome by generating read depths across the genome using a gaussian distribution. Genomic regions with very high read depths are filtered out to correct for biases resulting from high GC content or repeatitive regions. CNVs like deletions and duplications are then identified in regions having higher/lower than average read depths.

The script Plotting_read_depths.py is used for plotting read depths across specific genomic intervals for a given scaffold in a genome.
