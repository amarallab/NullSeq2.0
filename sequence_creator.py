'''
File: sequence_creator.py
Author: Sophia Liu
Description: The contains a sequence creator object. Creates sequences of known CUB
given AA, GC, and length contraints.
'''

import numpy as np
import nullseq_functions as nsf

class InputError(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return (repr(self.msg))

class Sequence_Creator(object):

    def __init__(self, codonnumber, length, size, AA=None, GC=None, n=11):
        ''' Initialization of the object

        Parameters
        ----------
        codonnumber : int, [20, 61]
            the number of codons to be used in each sequence

        length : int
            length of sequences
            elements : int

        size : int
            number of random sequences to create

        AA : dict, optional
            AA probability to be used in each sequence
            Keys : str
                AA
            Values : float
                p(AA), [0,1]

        GC : float, optional, [0,1]
            GC content to be used in each sequence



        n : int
            The codon table index according to NCBI
        '''

        self.n = n
        self.pcodonlist = None
        self.codontablelist = []
        self.codonlistlist = []

        if codonnumber <= 61 and codonnumber >= 20:
            self.codonnumber = codonnumber
        else:
            raise InputError('codonnumber must be in the range [20, 61]')

        if not isinstance(length, int):
            raise InputError ("sequence length must be interger value")
        elif length < 1:
            raise InputError ("sequence length must be greater than 0")
        else:
            self.length = length

        if not isinstance(size, int):
            raise InputError ("number of sequences must be interger value")
        elif size < 1:
            raise InputError ("number of sequences must be greater than 0")
        else:
            self.size = size

        if AA is None:
            self.AA = AA
        else:
            self.AA = nsf.check_correct_pAA(AA, n)

        if GC is None:
            self.GC = GC

        else:
            if GC >=0 and GC <= 1:
                self.GC = GC
            else:
                raise InputError('GC content must be in the range [0,1]')

        if self.AA is not None and self.GC is not None:
            iteration = 0
            while len(self.codontablelist) < self.size:
                if iteration > self.size * 1000 * (62-codonnumber):
                    raise InputError('GC content and AA content not compatible')
                else:
                    codontable = nsf.create_codon_table(self.codonnumber, n=self.n)
                    if nsf.evaluate_possibility(self.AA, codontable, self.GC):
                        self.codontablelist.append(codontable)
                        self.codonlistlist.append([])
                        for AA in self.codontablelist[-1].keys():
                            self.codonlistlist[-1] += self.codontablelist[-1][AA]
                    else:
                        pass
                iteration += 1
        else:
            while len(self.codontablelist) < self.size:
                self.codontablelist.append(nsf.create_codon_table(self.codonnumber, n=self.n))
                self.codonlistlist.append([])
                for AA in self.codontablelist[-1].keys():
                    self.codonlistlist[-1] += self.codontablelist[-1][AA]

        self.calculate_pcodon()

    def add_pAA(self, AA):
        ''' add p(AA) dictionary to object

        Parameters
        ----------

        AA : dict
            AA probability to be used in each sequence
            Keys : str
                AA
            Values : float
                p(AA), [0,1]
        '''

        self.AA = nsf.check_correct_pAA(AA, self.n)
        self.calculate_pcodon()

    def add_GC(self, GC):
        ''' add GC content parameter to object

        Parameters
        ----------
        GC : float, [0,1]
            GC content to be used in each sequence
        '''

        if GC >=0 and GC <= 1:
            self.GC = GC
            self.calculate_pcodon()
        else:
            raise InputError('GC content must be in the range [0,1]')

    def calculate_pcodon(self):
        ''' creates list of p(codon) used to create random sequences
        '''
        AAforCodon = nsf.AA_for_Codon(self.n)
        self.pcodonlist = []
        if self.AA is None:
            if self.GC is None:
                # no AA or GC constraints
                self.pcodonlist = [[1/self.codonnumber]*self.codonnumber]*self.size
            else:
                #no AA constaints, only GC
                for i in range(self.size):
                    pcodon = nsf.pcodon_wGC(self.codontablelist[i], self.GC)
                    self.pcodonlist.append([pcodon[codon] for codon in self.codonlistlist[i]])
        else:
            if self.GC is None:
                ## AA constaints no GC constraints
                for  i in range(self.size):
                    self.pcodonlist.append([])
                    for codon in self.codonlistlist[i]:
                        self.pcodonlist[-1].append(self.AA[AAforCodon[codon]]/\
                            len(self.codontablelist[i][AAforCodon[codon]]))
            else:
                # both GC and AA constriants
                for i in range(self.size):
                    pcodon = nsf.pcodon_wGCAA(self.codontablelist[i], self.GC, self.AA)
                    self.pcodonlist.append([pcodon[codon] for codon in self.codonlistlist[i]])

    def make_random_sequences(self):
        ''' Creates list of random sequences according to criteria

        Returns:
        --------
        seqlist : list
            string, random sequence
        '''
        seqlist = []

        for i in range(self.size):
            seqlist.append('ATG' + ''.join(np.random.choice(self.codonlistlist[i],
                                                            replace=True,
                                                            p=self.pcodonlist[i],
                                                            size=self.length)) + 'TGG')
        return seqlist





