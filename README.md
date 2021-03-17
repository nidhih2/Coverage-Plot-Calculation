# Coverage-Plot-Calculation

Calculates kmer-coverage plot and returns the number of kmers after removal of infrequent kmers using a calculated threshold value. 

## Syntax
python3 coverage_calc.py -reads <read file> -kmer_len <insert desired length of kmer> -plot_name <kmer-coverage plot name>
  
## Example usage
Input:
ATGCTAG
TGCTAGT
GCTAGTC

Output:
Based on the chosen kmer_len,
1) Total count of kmers 
2) kmer-coverage plot 
3) Threshold value for removal of sequencing errors
4) Total count of kmers after applying threshold value
