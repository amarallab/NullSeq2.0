
import collections
from Bio.Seq import Seq


def get_ATCG_Count(seq):
    '''Determines the number of each nucleotide in the sequence

    Parameters:
    -----------
    seq : str
        nucleotide sequence to be evaluated

    Returns:
    --------
    ATCGcount : dictionary
        Keys : str
            A, T, C, G
        Values : float
            number of corresponding nucleotide in the sequence

    Note:
    -----
    input sequence should include start and stop Codons
    calculated values do not take into accout the start and stop codons
    '''

    ATCGcount = collections.Counter(list(seq)[3:-3])

    return ATCGcount


def get_ATCG_Freq(seq):
    '''Determines the frequence each nucleotide is used in the SequenceGenerato

    Parameters:
    -----------
    seq : str
        nucleotide sequence to be evaluated

    Returns:
    --------
    ATCGfeq : dictionary
        Keys : str
            A, T, C, G
        Values : float
            frequency of the corresponding nucleotide

    Note:
    -----
    input sequence should include start and stop Codons
    calculated values do not take into accout the start and stop codons
    '''

    ATCGfreq = get_ATCG_Count(seq)
    for nuc in ATCGfreq:
        ATCGfreq[nuc] = ATCGfreq[nuc]/(len(seq)-6)

    return ATCGfreq

def get_GC_Count(seq):
    '''Determines the number GC nucleotides in the sequence

    Parameters:
    -----------
    seq : str
        nucleotide sequence to be evaluated

    Returns:
    --------
    GCcount : float
        the number of GC nucleotides in the sequence

    Note:
    -----
    input sequence should include start and stop Codons
    calculated values do not take into accout the start and stop codons
    '''

    ATCGcount = get_ATCG_Count(seq)
    GCcount = ATCGcount['G'] + ATCGcount['C']

    return GCcount

def get_GC_Freq(seq):
    '''Determines the frequency of GC nucleotide in the sequence

    Parameters:
    -----------
    seq : str
        nucleotide sequence to be evaluated

    Returns:
    --------
    GCfreq : float
        the frequency of GC nucleotides in the sequence

    Note:
    -----
    input sequence should include start and stop Codons
    calculated values do not take into accout the start and stop codons
    '''

    GCfreq = get_GC_Count(seq)/(len(seq)-6)

    return GCfreq

def get_AA_Count(seq, n, Nuc=True):
    '''Determines the number of each Amino acids used in the sequence

    Parameters:
    -----------
    seq : str
        The sequence to be evaluated
        if Nuc is True:
            seq is the nucleotide sequence
        if Nuc is False:
            seq is the amino acid sequence

    n : int
        Codon table to be used

    Nuc : bool (optional)
        True:
            evaulates seq as a nucleotide sequence
        False:
            evaulates seq as a Amino acid sequence

    Returns:
    --------
    AAUsageDict : dictionary
        Keys : str
            Amino acid symbol
        Values : float
            the number of the corresponding AA in the sequence
    '''

    if Nuc:
        AA = str(Seq(seq).translate(table=n)[1:-1])
    else:
        AA = seq[1:]

    AAUsageDict = collections.Counter(AA)

    return AAUsageDict

def get_AA_Freq(seq, n, nucleotide=True):
    '''Determines the frequency of each Amino acids used in the sequence

    Parameters:
    -----------
    seq : str
        The sequence to be evaluated
        if nucleotide is True:
            seq is the nucleotide sequence
        if nucleotide is False:
            seq is the amino acid sequence

    n : int
        Codon table to be used

    nucleotide : bool (optional)
        True:
            evaulates seq as a nucleotide sequence
        False:
            evaulates seq as a Amino acid sequence

    Returns:
    --------
    AAUsageDict : dictionary
        Keys : str
            Amino acid symbol
        Values : float
            the frequency of the corresponding amino acid
    '''

    AAUsageDict = get_AA_Count(seq, n, Nuc=nucleotide)

    for AA in AAUsageDict.keys():
        if nucleotide:
            AAUsageDict[AA] = AAUsageDict[AA]/((len(seq)-6)/3)
        else:
            AAUsageDict[AA] = AAUsageDict[AA]/(len(seq)-1)

    return AAUsageDict
