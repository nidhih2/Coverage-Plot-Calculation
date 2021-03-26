# Coverage-Plot-Calculation

Calculates kmer-coverage plot and returns the number of kmers after removal of infrequent kmers using a calculated threshold value (local minimum). 

# Syntax
python3 coverage_calc.py -reads <read file> -kmer_len <insert desired length of kmer> -plot_name <kmer-coverage plot name>
  
# Arguments
```-reads``` : Provide a read file\
```-kmer_len``` : Provide an integer value\
```-plot_name``` : Provide a name for the kmer-coverage plot
  
# Assumption
Takes the minimum of calculated local minimum and plots on the sns.distplot graph
  
# Example usage
## Input:
ATGCTAG\
TGCTAGT\
GCTAGTC

## Output:
Based on the chosen kmer_len,
1) Total count of kmers 
2) kmer-coverage plot with threshold value (indicated by a vertical line)
3) Threshold value for removal of sequencing errors
4) Total count of kmers after applying threshold value

## Run time: 
Takes about 118 to 130s for a desired 31 to 101 kmer length.
Used time module to estimate run time 

