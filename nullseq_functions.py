'''
File: NullSeq_Functions.py
Author: Sophia Liu
Date: May 16, 2016
Description: This script will generate random nucleotide sequences with a
             given GC content either based on a nucleotide or amino acid
             sequence or amino acid usage probabitliy
'''

import collections
from Bio.Seq import Seq
from Bio.Data import CodonTable
import random
import math
import bisect
import numpy as np
import pandas as pd
import os.path

def Codons_for_AA(n):
    ''' Determines the codons used to code an amino acid given a codon table

    Parameters
    ----------
    n : int
        the codon table index according to NCBI

    Returns
    -------
    CodonforAA : dictionary
        Key : str
            Amino acids in Single letter notation (caps)
        Values : str
            Three letter codon (in caps)
    '''


    GivenCodonTable = CodonTable.unambiguous_dna_by_id[n]
    nucleotides = ['A', 'T', 'C', 'G']
    CodonforAA = {}
    for first in nucleotides:
        for second in nucleotides:
            for third in nucleotides:
                Codon = first + second + third
                if Codon not in CodonTable.unambiguous_dna_by_id[n].stop_codons:
                    if GivenCodonTable.forward_table[Codon] in CodonforAA.keys():
                        CodonforAA[GivenCodonTable.forward_table[Codon]].append(Codon)
                    else:
                        CodonforAA[GivenCodonTable.forward_table[Codon]] = [Codon]
                else:
                    pass
    return CodonforAA

def AA_for_Codon(n):
    ''' Determines the AA for a codon given a codon table

    Parameters
    ----------
    n : int
        the codon table index according to NCBI

    Returns
    -------
    dictionary
        Key : str
            Three letter codon (in caps)
        Values : str
            Amino acids in Single letter notation (caps)
    '''
    CodonforAA = Codons_for_AA(n)
    codonlist = []
    AAlist = []
    for AA, clist in CodonforAA.items():
        codonlist += clist
        AAlist += [AA]*len(clist)

    return dict(zip(codonlist, AAlist))

def create_codon_table(ncodons, n):
    ''' Randomly creates abridged codon table of
        ncodon number of codons

    Parameters
    ----------
    ncodons : int, [20, 61]
        the number of codons in the new codon table

    n : int
        the codon table index according to NCBI

    Returns
    -------
    CofAA : dict
        new codon table dictinoary with ncodon number of
        codons
        Key : str
            Amino acids in Single letter notation (caps)
        Values : str
            Three letter codon (in caps)
    '''

    CodonforAA = Codons_for_AA(n)
    ListofCodons = []
    ListofAA = []
    for AA in CodonforAA.keys():
        ListofAA.append(AA)
        ListofCodons.append(CodonforAA[AA])
    number_ListofCodons = len([item for sublist in ListofCodons for item in sublist])
    for n in range(number_ListofCodons-ncodons):
        possibleindexlist = []
        for k in range(len(ListofCodons)):
            if len(ListofCodons[k]) > 1:
                possibleindexlist.append(k)
            else:
                pass
        removeAAindex = random.choice(possibleindexlist)
        ListofCodons[removeAAindex].pop(random.randint(0, len(ListofCodons[removeAAindex])-1))

    CofAA = {}
    for l in range(len(ListofAA)):
        CofAA[ListofAA[l]] = ListofCodons[l]

    return CofAA

def check_correct_pAA(pAA, n):
    ''' Corrects imput p(AA) to ensure
        sum(p(AA)) = 1

    Parameters
    ----------
    AA : dict
        AA probability to be used in each sequence
        Keys : str
            AA
        Values : float
            p(AA), [0,1]

    n : int
        The codon table index according to NCBI

    Returns
    -------
    dictionary
        Key : str
            Amino acids in Single letter notation (caps)
        Values : float
            p(AA), [0,1]
    '''

    pAAlist = []
    nAAlist = []
    AAlist = list(Codons_for_AA(n).keys())
    oAAlist = list(pAA.keys())
    for AA in AAlist:
        if AA in oAAlist:
            nAAlist.append(AA)
            pAAlist.append(pAA[AA])
        else:
            nAAlist.append(AA)
            pAAlist.append(0)
    pAAlist = np.array(pAAlist)/sum(pAAlist)
    return dict(zip(nAAlist, pAAlist))

def minmax_GC_for_AA(codontable):
    ''' Determines the maximum or minimum number
    of G/C nucleotide for a given AA

    Parameters:
    -----------
    n : int
        NCBI translation table

    Returns:
    --------
    minGCDict : dictionary
        Keys : str
            AA
        Values : int
            smallest number of G/C nucleotides

    maxGCDict : dictionary
        Keys : str
            AA
        Values : int
            largest number of G/C nucleotides
    '''

    CodonGCContent = GC_Content_for_Codons(codontable)
    maxGCDict = {}
    minGCDict = {}
    for AA in codontable.keys():
        templist = [CodonGCContent[codon] for codon in codontable[AA]]
        maxGCDict[AA] = max(templist)
        minGCDict[AA] = min(templist)

    return (minGCDict, maxGCDict)

def GC_Content_for_Codons(codontable):
    ''' Determines the number of G/C nucletides for all codons in the
    codon table

    Parameters
    ----------
    codontable : dict
        AA to codon table
        Key : str
            Amino acids in Single letter notation (caps)
        Values : str
            Three letter codon (in caps)

    Returns
    -------
    CodonGCContent : dictionary
        Key : str
            Three letter codon (in caps)
        Values : int
            Number of G/C nucleotides in the codon
    '''

    CodonGCContent = collections.defaultdict(int)
    for AA in codontable:
        for codon in codontable[AA]:
            codonusage = collections.Counter(codon)
            CodonGCContent[codon] = codonusage['G'] + codonusage['C']

    return CodonGCContent

def get_maxmin_GC_count(AAfreq, codontable):
    '''Determines the maximum and minimum GC ratio for
    specified AA compostion

    Parameters:
    -----------
    AAfreq : dictionary
        Keys : str
            AA
        Values : float
            frequency of AA in sequence

    codontable : dict
        AA to codon table
        Key : str
            Amino acids in Single letter notation (caps)
        Values : str
            Three letter codon (in caps)

    Returns:
    --------
    low : float
        lowest possible GC ratio

    high : float
        highest possible GC ratio
    '''

    (minGCDict, maxGCDict) = minmax_GC_for_AA(codontable)
    high = 0
    low = 0
    for AA in AAfreq:
        high += AAfreq[AA]*maxGCDict[AA]/3
        low += AAfreq[AA]*minGCDict[AA]/3

    return (low, high)

def evaluate_possibility(AAfreq, codontable, GC):
    '''Determines whether is it possible to obtain the GC content
    for the given AA compostion

    Parameters:
    -----------
    AAfreq : dictionary
        Keys : str
            AA
        Values : float
            frequency of AA in sequence

    GC : float
        desired GC content [0, 1]

    codontable : dictionary
        Keys : str
            AA
        Values : list
            list of synonymous codons

    Returns:
    --------
    bool
        True : if the GC content is obtainable
        False: if the GC content is not obtainable
    '''

    (L, H) = get_maxmin_GC_count(AAfreq, codontable)
    if GC > L and GC < H:
        return True
    else:
        return False

def pcodon_wGC(codontable, GC):
    '''creates nucleotide sequences with desired Gc content

    Parameters:
    -----------
    codontable : dictionary
        Keys : str
            AA
        Values : list
            list of synonymous codons

    GC : float
        Desired GC content [0,1]

    Return:
    -------
    dict
        Key : str
            three letter codon
        Values : float
            p(codon)
    '''
    flatcodonlist = [item for key, value in codontable.items() for item in value]
    beta = get_beta(GC, codontable, True)
    p = Probability_Given_Beta(beta, codontable, True)
    plist = [p[codon] for codon in flatcodonlist]
    return dict(zip(flatcodonlist, plist))

def pcodon_wGCAA(codontable, GC, AA):
    '''creates nucleotide sequences with desired Gc content

    Parameters:
    -----------
    codontable : dictionary
        Keys : str
            AA
        Values : list
            list of synonymous codons

    GC : float
        Desired GC content [0,1]

    AAprob : dictionary, optional
        Keys : str
            AA
        Values : float
            p(AA)

    Return:
    -------
    dict
        Key : str
            three letter codon
        Values : float
            p(codon)
    '''
    beta = get_beta(GC, codontable, False, AAprob=AA)
    p = Probability_Given_Beta(beta, codontable, False)
    codonlist = []
    plist = []
    for A, clist in codontable.items():
        codonlist += clist
        plist += list(np.array([p[codon] for codon in clist])*AA[A])
    return dict(zip(codonlist, plist))

def get_beta(given, codontable, noAA, AAprob=None):
    '''Determines the value of beta given the GC content of the sequence

    Parameters:
    -----------
    given : float
        GC content can be speficied [0,1]

    codontable : dictionary
        Keys : str
            AA
        Values : list
            list of synonymous codons
    noAA : bool
        True if AA constraints are considered
        False if no AA constraints are considered

    AAprob : dictionary, optional
        Keys : str
            AA
        Values : float
            frequency of AA in sequence

    Returns:
    -------
    beta : float
        the constant in the equation P = exp(-beta*N)/Z
    '''

    m = given

    rl = -20
    rr = 20

    avl = compute_average_GC(rl, codontable, noAA, AAprob=AAprob) - m
    avr = compute_average_GC(rr, codontable, noAA, AAprob=AAprob) - m

    if avl*avr > 0:
        if avl < 0:
            return rl
        else:
            return rr
    else:

        rm = (rl+rr)/2
        avv = compute_average_GC(rm, codontable, noAA, AAprob=AAprob) - m

        while math.fabs(avv) > 1e-7:
            if avv*avl > 0:
                rl = rm
            else:
                rr = rm
            rm =  0.5*(rl+rr)
            avv = compute_average_GC(rm, codontable, noAA, AAprob=AAprob) - m

    return rm

def compute_average_GC(b, codontable, noAA, AAprob=None):
    ''' Computes the expected value of the GC content of the sequences
        generated with a given beta

    Parameters:
    -----------
    b : float
        the constant in the equation P = exp(-beta*N)/Z

    codontable : dictionary
        Keys : str
            AA
        Values : list
            list of synonymous codons
    noAA : bool
        True if AA constraints are considered
        False if no AA constraints are considered

    AAprob : dictionary, optional
        Keys : str
            AA
        Values : float
            frequency of AA in sequence

    Returns:
    -------
    float
    expected value of the GC content of the random sequences
    '''

    CodonGC = GC_Content_for_Codons(codontable)
    GCnumber = []
    Probability = Probability_Given_Beta(b, codontable, noAA)

    if noAA:
        flatcodonlist = [item for key, value in codontable.items() for item in value]
        for c in flatcodonlist:
            GCnumber.append(Probability[c]*CodonGC[c]/3)
        return sum(GCnumber)

    else:
        for AA, codons in codontable.items():
            for codon in codons:
                GCnumber.append(Probability[codon]*CodonGC[codon]/3*AAprob[AA])
        return sum(GCnumber)

def Probability_Given_Beta(beta, codontable, noAA):
    ''' Computes p(codon) for a given beta and codon table

    Parameters:
    -----------
    beta : float
        the constant in the equation P = exp(-beta*N)/Z

    codontable : dictionary
        Keys : str
            AA
        Values : list
            list of synonymous codons
    noAA : bool
        True if AA constraints are considered
        False if no AA constraints are considered

    Returns:
    -------
    probs : dict
        Key : str
            theree letter codon
        Values : float
            p(codon)
    '''

    CodonGC = GC_Content_for_Codons(codontable)
    probs = {}

    if noAA:
        flatcodonlist = [item for key, value in codontable.items() for item in value]
        probs = np.array([math.exp(-beta*CodonGC[codon]) for codon in flatcodonlist])
        probs = probs/sum(probs)
        probs = dict(zip(flatcodonlist, probs))
        return probs

    else:
        for AA, codons in codontable.items():
            Z = sum([math.exp(-beta*CodonGC[codon]) for codon in codons])
            for codon in codons:
                probs[codon] = math.exp(-beta*CodonGC[codon])/Z
        return probs
