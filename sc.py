'''
Author : Jay Sorathiya
Running this file finds the position of each modified amino acid residue in a
input.tsv from a fasta file of Protein Sequesnces 
'''

#### Import libraries and dependencies #########
import sys
import pandas as pd
from itertools import groupby

#### FUNCTIONS #########
def fasta_iter(fasta_name):
    """
    Parses a Fasta File, given a fasta file. 
    Yield tuples of header, sequence

    Parameters:
    -----------
    fasta_name : a file with .fasta extension
    Note: The Fasta file should be in the same directory
    """

    # Open the file outside
    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)
        
#### Function to clean data
def data_clean(inputf):
    """
    Function to clean input Data
    Takes a dataframe as input

    Function to perform the search operation

    Parameters:
    -----------
    inputlist : Python list of input sequences
    sequences : Python list of Protein sequences
    arg_3 : Name of the output file

    """
    inputlist = []
    exact_input = []

    for i in inputf.FullPeptideName:
        inputlist.append(i)
        exact_input.append(i)

    #Removing UniMod
    inputlist = [x.replace('UniMod','') for x in inputlist]

    #Removing '('
    inputlist = [x.replace('(','') for x in inputlist]

    #Removing ')'
    inputlist = [x.replace(')','') for x in inputlist]

    #Removing ':'
    inputlist = [x.replace(':','') for x in inputlist]

    #Removing Digits
    for i in range(len(inputlist)):
        inputlist[i] = ''.join(i for i in inputlist[i] if not i.isdigit())

    return inputlist, exact_input

def join_similars(out):
    """
    Function to clean output by joining similars
    Takes a dataframe with three cloumns 'ProteinName','FullPeptideName',
    'Modified_amino_acid_Positions' as input

    Parameters:
    -----------
    out : Pandas dataframe of output        
    """
    count = 0
    for i in range(len(out.FullPeptideName)-1):
        if out.FullPeptideName[i] == out.FullPeptideName[i+1]:
            out.Modified_amino_acid_Positions[i+1] = ';'.join(map(str,out.Modified_amino_acid_Positions[i:i+2]))
            out = out.drop(labels = int(i), axis=0)   ### Drops the duplicate entries
            count+=1
        
    print('{} Modified amino acid positions were similar, hence adjusted'.format(count))
    return out

#### Main function to search indexes
def Index_Search(inputlist, heads, exact_input, sequences, arg_3):
    """
    Function to perform the search operation

    Parameters:
    -----------
    inputlist : Python list of input sequences
    sequences : Python list of Protein sequences
    arg_3 : Name of the output file

    """
    Prot_name = []
    Fpept_name = []
    Position = []

    for i in range(len(inputlist)):
    
        for j in range(len(sequences)):
            
            if inputlist[i] in sequences[j]:
                Prot_name.append(heads[j])
                Fpept_name.append(exact_input[i])
                Position.append(sequences[j].index(inputlist[i]))  ## Finding the Index of Full peptide in Protien Sequence
    
    #Saving the results to a dataFrame    
    data = {'ProteinName': Prot_name ,'FullPeptideName': Fpept_name ,
           'Modified_amino_acid_Positions': Position}
    outputfile = pd.DataFrame(data)
    
    #### Removing the duplicates and joining them together , semi-colon separated
    outputfile = join_similars(outputfile)  ### Pre-defined function join_similars
    
    #### Reset Index after removal and saving the output to a csv file
    outputfile = outputfile.reset_index(drop=True)
    outputfile.to_csv(arg_3, sep='\t')
        
    return outputfile


#### MAIN #########
#### Save all the arguments to a variable 
arg_1 = sys.argv[1]
arg_2 = sys.argv[2]
arg_3 = str(sys.argv[3])


#### Read input data and Fasta file for Protein sequences

####################### Argument 1
### Input file containing FullPeptideName
std_input = pd.read_csv(arg_1,sep='\t')

# Filtering the Modified Protein sequences
inputf = std_input[std_input.FullPeptideName.str.contains('UniMod', regex= True, na=False)]

inputlist, exact_input = data_clean(inputf) ## Function data_clean to clean the data


####################### Argument 2
### Fasta file of protein sequences
fiter = fasta_iter(arg_2)  ## Function fasta_iter 
sequences = []
heads = []
for ff in fiter:
    headerStr, seq = ff
    sequences.append(seq)
    heads.append(headerStr)

    
#### Function to run program
outputfile = Index_Search(inputlist, heads, exact_input, sequences, arg_3)


#### Program complete statement
print('Required output file {} has been generated in the current directory'.format(arg_3))