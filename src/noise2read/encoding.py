# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-02-16 11:01:06
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-05-14 11:00:24
# The codes in this section for features construction are modified from https://github.com/Bonidia/MathFeature, if you re-use these codes or methods, please also cite the original publication of Robson P Bonidia, Douglas S Domingues, Danilo S Sanches, André C P L F de Carvalho, MathFeature: feature extraction package for DNA, RNA and protein sequences based on mathematical descriptors, Briefings in Bioinformatics, 2021; bbab434, https://doi.org/10.1093/bib/bbab434.


import numpy as np
import statistics
from scipy.fftpack import fft
import math

class FourierTransform():
    def __init__(self, read):
        self.read = read.upper()

    def feature_extraction(self, spectrum, spectrumTwo):
        features = []
        average = sum(spectrum)/len(spectrum)
        features.append(average)
        ###################################
        median = np.median(spectrum)
        features.append(median)
        ###################################
        maximum = np.max(spectrum)
        features.append(maximum)
        ###################################
        minimum = np.min(spectrum)
        features.append(minimum)
        ###################################
        if average != 0:
            peak = (len(spectrum)/3)/(average)
        else:
            peak = 0
        features.append(peak)
        ###################################
        if np.mean(spectrumTwo) != 0:
            peak_two = (len(spectrumTwo)/3)/(np.mean(spectrumTwo))
        else:
            peak_two = 0
        features.append(peak_two)
        ###################################
        standard_deviation = np.std(spectrum) # standard deviation
        features.append(standard_deviation)
        ###################################
        if len(spectrum) >= 2:
            standard_deviation_pop = statistics.stdev(spectrum) # population sample standard deviation 
        else:
            standard_deviation_pop = 0
        features.append(standard_deviation_pop)
        ###################################
        percentile15 = np.percentile(spectrum, 15)
        features.append(percentile15)
        ###################################
        percentile25 = np.percentile(spectrum, 25)
        features.append(percentile25)
        ###################################
        percentile50 = np.percentile(spectrum, 50)
        features.append(percentile50)
        ###################################
        percentile75 = np.percentile(spectrum, 75)
        features.append(percentile75)
        ###################################
        amplitude = maximum - minimum
        features.append(amplitude)
        ###################################
        # mode = statistics.mode(spectrum)
        ###################################
        if len(spectrum) >= 2:
            variance = statistics.variance(spectrum)
        else:
            variance = 0
        features.append(variance)
        ###################################
        interquartile_range = np.percentile(spectrum, 75) - np.percentile(spectrum, 25)
        features.append(interquartile_range)
        ###################################
        semi_interquartile_range = (np.percentile(spectrum, 75) - np.percentile(spectrum, 25))/2 
        features.append(semi_interquartile_range)
        ###################################
        if average != 0:
            coefficient_of_variation = standard_deviation/average
        else:
            coefficient_of_variation = 0
        features.append(coefficient_of_variation)
        ###################################
        if standard_deviation != 0:
            skewness = (3 * (average - median))/standard_deviation
        else:
            skewness = 0
        features.append(skewness)   
        ###################################
        if (2 * (np.percentile(spectrum, 90) - np.percentile(spectrum, 10))) == 0:
            kurtosis = 0
        else:
            kurtosis = (np.percentile(spectrum, 75) - np.percentile(spectrum, 25)) / (2 * (np.percentile(spectrum, 90) - np.percentile(spectrum, 10))) 
        # print(kurtosis)
        features.append(kurtosis)
        ###################################
        return features

    def binary_fourier(self):
        # if len(self.read) > 0:
        spectrum = []
        spectrumTwo = []
        A = []
        C = []
        T = []
        G = []
        for nucle in self.read:
            if nucle == "A":
                A.append(1)
            else:
                A.append(0)
            if nucle == "C":
                C.append(1)
            else:
                C.append(0)
            if nucle == "T" or nucle =="U":
                T.append(1)
            else:
                T.append(0)
            if nucle == "G":
                G.append(1)
            else:
                G.append(0)
        FA = fft(A)
        FC = fft(C)
        FT = fft(T)
        FG = fft(G)
        for i in range(len(self.read)):
            specTotal = (abs(FA[i])**2) + (abs(FC[i])**2) + (abs(FT[i])**2) + (abs(FG[i])**2)
            specTwo = (abs(FA[i])) + (abs(FC[i])) + (abs(FT[i])) + (abs(FG[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        features = self.feature_extraction(spectrum, spectrumTwo)
        return features

    def zcurve_fourier(self):
        # if len(self.read) > 0:
        spectrum = []
        spectrumTwo = []
        ###################################
        ###################################
        R = 0 # x[n] = (An + Gn) − (Cn + Tn) ≡ Rn − Yn
        Y = 0
        M = 0 # y[n] = (An + Cn) − (Gn + Tn) ≡ Mn − Kn
        K = 0
        W = 0 # z[n] = (An + Tn) − (Cn + Gn) ≡ Wn − Sn
        S = 0
        ###################################
        ###################################
        x = []
        y = []
        z = []
        for nucle in self.read:
            if nucle == "A" or nucle == "G":
                R += 1
                x.append((R)-(Y))
            else:
                Y += 1
                x.append((R)-(Y))
            if nucle == "A" or nucle == "C":
                M += 1
                y.append((M)-(K))
            else:
                K += 1
                y.append((M)-(K))
            if nucle == "A" or nucle == "T" or nucle == "U":
                W += 1
                z.append((W)-(S))
            else:
                S += 1
                z.append((W)-(S))
        # if len(x) == 0:
        #     x = [0]
        # if len(y) == 0:
        #     y = [0]
        # if len(z) == 0:
        #     z = [0]

        FX = fft(x)
        FY = fft(y)
        FZ = fft(z)
        for i in range(len(self.read)):
            specTotal = (abs(FX[i])**2) + (abs(FY[i])**2) + (abs(FZ[i])**2)
            specTwo = (abs(FX[i])) + (abs(FY[i])) + (abs(FZ[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        features = self.feature_extraction(spectrum, spectrumTwo)
        return features

    def integer_fourier(self):
        # if len(self.read) > 0:
        spectrum = []
        spectrumTwo = []
        integer = []
        for nucle in self.read:
            if nucle == "T" or nucle == "U":
                integer.append(0)
            elif nucle == "C":
                integer.append(1)
            elif nucle == "A":
                integer.append(2)
            else:
                integer.append(3)
        # if len(integer) == 0:
        #     integer = [0]
        FI = fft(integer)
        for i in range(len(self.read)):
            specTotal = (abs(FI[i])**2)
            specTwo = (abs(FI[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        features = self.feature_extraction(spectrum, spectrumTwo)
        return features

    def real_fourier(self):
        # if len(self.read) > 0:
        spectrum = []
        spectrumTwo = []
        real = []
        for nucle in self.read:
            if nucle == "T" or nucle == "U":
                real.append(1.5)
            elif nucle == "C":
                real.append(0.5)
            elif nucle == "A":
                real.append(-1.5)
            else:
                real.append(-0.5)
        # if len(real) == 0:
        #     real = [0]
        FR = fft(real)
        for i in range(len(self.read)):
            specTotal = (abs(FR[i])**2)
            specTwo = (abs(FR[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        features = self.feature_extraction(spectrum, spectrumTwo)
        return features

    def eiip_fourier(self):
        # if len(self.read) > 0:
        spectrum = []
        spectrumTwo = []
        eiip = []
        for nucle in self.read:
            if nucle == "T" or nucle == "U":
                eiip.append(0.1335)
            elif nucle == "C":
                eiip.append(0.1340)
            elif nucle == "A":
                eiip.append(0.1260)
            else:
                eiip.append(0.0806)
        # if len(eiip) == 0:
        #     eiip = [0]
        Feiip = fft(eiip)
        for i in range(len(self.read)):
            specTotal = (abs(Feiip[i])**2)
            specTwo = (abs(Feiip[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        features = self.feature_extraction(spectrum, spectrumTwo)
        return features

    def complex_number_fourier(self):
        # if len(self.read) > 0:
        spectrum = []
        spectrumTwo = []
        complexNumber = []
        for nucle in self.read:
            if nucle == "T" or nucle == "U":
                complexNumber.append(1-1j)
            elif nucle == "C":
                complexNumber.append(-1+1j)
            elif nucle == "A":
                complexNumber.append(1+1j)
            else:
                complexNumber.append(-1-1j)
        # if len(complexNumber) == 0:
        #     complexNumber = [0]
        FcomplexNumber = fft(complexNumber)
        for i in range(len(self.read)):
            specTotal = (abs(FcomplexNumber[i])**2)
            specTwo = (abs(FcomplexNumber[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        features = self.feature_extraction(spectrum, spectrumTwo)
        return features

    def atomic_number_fourier(self):
        # if len(self.read) > 0:
        spectrum = []
        spectrumTwo = []
        atomicNumber = []
        for nucle in self.read:
            if nucle == "T" or nucle == "U":
                atomicNumber.append(66)
            elif nucle == "C":
                atomicNumber.append(58)
            elif nucle == "A":
                atomicNumber.append(70)
            else:
                atomicNumber.append(78)
        # if len(FatomicNumber) == 0:
        #     FatomicNumber = [0]
        FatomicNumber = fft(atomicNumber)
        for i in range(len(self.read)):
            specTotal = (abs(FatomicNumber[i])**2)
            specTwo = (abs(FatomicNumber[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        features = self.feature_extraction(spectrum, spectrumTwo)
        return features

    def Encoding(self):
        # ed1 = self.binary_fourier()
        ed2 = self.atomic_number_fourier()
        # ed3 = self.complex_number_fourier()
        # ed4 = self.eiip_fourier()
        # ed5 = self.integer_fourier()
        # ed6 = self.zcurve_fourier()
        # ed7 = self.real_fourier()
        # encoding_vector = ed1 + ed2 + ed3 + ed4 + ed5 + ed6 + ed7
        # return encoding_vector
        return ed2


#######################################################################################################
class ChaosGame():
    def __init__(self, read, max_len, kmer_freq):
        '''
        kmer_freq: Frequency of k-mer (e.g., 3, 4)
        '''
        self.read = read.upper()
        self.max_len = max_len
        # self.padding_mode = padding_mode
        self.kmer_freq = kmer_freq
        
    def chunksTwo(self, seq, win):
        seqlen = len(seq)
        for i in range(seqlen):
            j = seqlen if i+win>seqlen else i+win
            yield seq[i:j]
            if j==seqlen: break
        return

    def frequency_chaos(self, k):
        # if len(self.read) >= k:
        mapping = []
        kmer = {}
        totalWindows = (len(self.read) - k) + 1 # (L - k + 1)
        for subseq in self.chunksTwo(self.read, k):
            # print(subseq)
            if subseq in kmer:
                kmer[subseq] = kmer[subseq] + 1
            else:
                kmer[subseq] = 1
        for subseq in self.chunksTwo(self.read, k):

            mapping.append(kmer[subseq]/totalWindows)
 
        padding = ((self.max_len - k + 1) - len(mapping))

        # mapping = np.pad(mapping, (0, padding), self.padding_mode)
        return mapping + [0]*padding

    def classifical_chaos(self):
        # if len(self.read) > 0:
        Sx = []
        Sy = []
        for nucle in self.read:
            if nucle == "A":
                Sx.append(1)
                Sy.append(1)
            elif nucle == "C":
                Sx.append(-1)
                Sy.append(-1)
            elif nucle == "T" or nucle =="U":
                Sx.append(-1)
                Sy.append(1)
            else:
                Sx.append(1)
                Sy.append(-1)
        CGR_x = [] 
        CGR_y = []
        for i in range(0, len(Sx)):
            if i == 0:
                CGR_x.append(0.5 * Sx[i])
                CGR_y.append(0.5 * Sy[i])
            else:
                CGR_x.append(0.5 * Sx[i] + 0.5 * CGR_x[i - 1])
                CGR_y.append(0.5 * Sy[i] + 0.5 * CGR_y[i - 1])           
        concat = CGR_x + CGR_y
        # if len(concat) == 0:
        #     concat = [0]
        padding = (self.max_len - len(Sx)) * 2
        # mapping = np.pad(concat, (0, padding), self.padding_mode)
        return concat + [0]*padding

    def classifical_chaos_fourier(self):
        # if len(self.read) > 0:
        spectrum = []
        spectrumTwo = []
        Sx = []
        Sy = []
        for nucle in self.read:
            if nucle == "A":
                Sx.append(1)
                Sy.append(1)
            elif nucle == "C":
                Sx.append(-1)
                Sy.append(-1)
            elif nucle == "T" or nucle =="U":
                Sx.append(-1)
                Sy.append(1)
            else:
                Sx.append(1)
                Sy.append(-1)
        CGR_x = [] 
        CGR_y = []
        for i in range(0, len(Sx)):
            if i == 0:
                CGR_x.append(0.5 * Sx[i])
                CGR_y.append(0.5 * Sy[i])
            else:
                CGR_x.append(0.5 * Sx[i] + 0.5 * CGR_x[i - 1])
                CGR_y.append(0.5 * Sy[i] + 0.5 * CGR_y[i - 1])      
        Fx = fft(CGR_x)
        Fy = fft(CGR_y)
        for i in range(len(Fx)):
            specTotal = (abs(Fx[i])**2) + (abs(Fy[i])**2)
            specTwo = (abs(Fx[i])) + (abs(Fy[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        features = FourierTransform(self.read).feature_extraction(spectrum, spectrumTwo)
        # else:
        #     features = np.zeros(19)
        return features

    def frequency_chaos_fourier(self, k):
        # if len(self.read) > 0:
        spectrum = []
        spectrumTwo = []
        mapping = []
        kmer = {}
        totalWindows = (len(self.read) - k) + 1 # (L - k + 1)
        for subseq in self.chunksTwo(self.read, k):
            # print(subseq)
            if subseq in kmer:
                kmer[subseq] = kmer[subseq] + 1
            else:
                kmer[subseq] = 1
        for subseq in self.chunksTwo(self.read, k):
            # print(kmer[subseq])
            # print(kmer[subseq]/totalWindows)
            mapping.append(kmer[subseq]/totalWindows)
        # if len(mapping) == 0:
        #     mapping = [0]
        Fmap = fft(mapping)
        for i in range(len(mapping)):
            specTotal = (abs(Fmap[i])**2)
            specTwo = (abs(Fmap[i]))
            spectrum.append(specTotal)
            spectrumTwo.append(specTwo)
        features = FourierTransform(self.read).feature_extraction(spectrum, spectrumTwo)
        # else:
        #     features = np.zeros(19)
        return features

    def Encoding(self):
        # ed1 = self.classifical_chaos()
        ed2 = self.classifical_chaos_fourier()
        # ed3 = self.frequency_chaos(self.kmer_freq)
        ed4 = self.frequency_chaos_fourier(self.kmer_freq)
        # encoding_vector = np.concatenate((ed1, ed2, ed3, ed4), axis=None)
        # return encoding_vector
        # return ed1 + ed2 + ed3 + ed4
        return ed2 + ed4

#######################################################################################################
class Entropy():
    def __init__(self, read, kmer, q=2):
        '''
        kmer: 'Range of k-mer, E.g., 1-mer (1) or 2-mer (1, 2) ...'
        q: 'Tsallis - entropic parameter q'
        '''
        self.read = read.upper()
        self.kmer = kmer
        self.q = q

    def chunks(seq, win, step):
        seqlen = len(seq)
        for i in range(0,seqlen,step):
            j = seqlen if i+win>seqlen else i+win
            yield seq[i:j]
            if j==seqlen: break
        return         

    def chunks_two(self, seq, win):
        seqlen = len(seq)
        for i in range(seqlen):
            j = seqlen if i+win>seqlen else i+win
            yield seq[i:j]
            if j==seqlen: break
        return

    def entropy_equation(self, entropy_type):
        information_entropy = []
        for k in range(1, self.kmer+1):
            probabilities = []
            kmer = {}
            total_windows = (len(self.read) - k) + 1 # (L - k + 1)
            for subseq in self.chunks_two(self.read, k):
                if subseq in kmer:
                    # print(subseq)
                    kmer[subseq] = kmer[subseq] + 1
                else:
                    kmer[subseq] = 1
            for key, value in kmer.items():
                # print(key)
                # print(value)
                probabilities.append(value/total_windows)
            if entropy_type == "Shannon" or entropy_type == "shannon":
                entropy_equation = [(p * math.log(p, 2)) for p in probabilities]
                entropy = -(sum(entropy_equation))
                information_entropy.append(entropy)
            elif entropy_type == "Tsallis" or entropy_type == "tsallis":
                # q = 2
                entropy_equation = [(p ** self.q) for p in probabilities]
                entropy =  (1/(self.q - 1)) * (1 - sum(entropy_equation))
                information_entropy.append(entropy)
        return information_entropy  

    def Encoding(self):
        ed1 = self.entropy_equation(entropy_type = "tsallis")
        ed2 = self.entropy_equation(entropy_type = "shannon")
        return ed1 + ed2

###############################################################################
class FickettScore():
    def __init__(self, read, read_type = 'DNA'):
        '''
        kmer_freq: Frequency of k-mer (e.g., 3, 4)
        '''
        self.read = read.upper()
        self.read_type = read_type

    def look_up_position_prob(self, value, base, position_para, position_prob, position_weight):
        """look up positional probability by base and value"""

        if float(value) < 0:
            return None
        for idx, val in enumerate(position_para):
            if float(value) >= val:
                return float(position_prob[base][idx]) * float(position_weight[base])

    def look_up_content_prob(self, value, base, content_para, content_prob, content_weight):

        """look up content probability by base and value"""

        if float(value) < 0:
            return None
        for idx, val in enumerate(content_para):
            if float(value) >= val:
                return float(content_prob[base][idx]) * float(content_weight[base])

    def fickett_value_orf(self, seq):

        """calculate Fickett value. Input is DNA sequence"""

        position_prob = {
            'A': [0.94, 0.68, 0.84, 0.93, 0.58, 0.68, 0.45, 0.34, 0.20, 0.22],
            'C': [0.80, 0.70, 0.70, 0.81, 0.66, 0.48, 0.51, 0.33, 0.30, 0.23],
            'G': [0.90, 0.88, 0.74, 0.64, 0.53, 0.48, 0.27, 0.16, 0.08, 0.08],
            'T': [0.97, 0.97, 0.91, 0.68, 0.69, 0.44, 0.54, 0.20, 0.09, 0.09]}

        position_weight = {'A': 0.26, 'C': 0.18, 'G': 0.31, 'T': 0.33}
        position_para = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0]

        content_prob = {
            'A': [0.28, 0.49, 0.44, 0.55, 0.62, 0.49, 0.67, 0.65, 0.81, 0.21],
            'C': [0.82, 0.64, 0.51, 0.64, 0.59, 0.59, 0.43, 0.44, 0.39, 0.31],
            'G': [0.40, 0.54, 0.47, 0.64, 0.64, 0.73, 0.41, 0.41, 0.33, 0.29],
            'T': [0.28, 0.24, 0.39, 0.40, 0.55, 0.75, 0.56, 0.69, 0.51, 0.58]}

        content_weight = {'A': 0.11, 'C': 0.12, 'G': 0.15, 'T': 0.14}
        content_para = [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.17, 0]

        if len(seq) < 2:
            return 0
        fickett_score = 0
        seq = seq.upper()
        total_base = len(seq)
        A_content = float(seq.count('A')) / total_base
        C_content = float(seq.count('C')) / total_base
        G_content = float(seq.count('G')) / total_base
        if self.read_type == 'DNA':
            T_content = float(seq.count('T')) / total_base
        elif self.read_type == 'RNA':
            T_content = float(seq.count('U')) / total_base

        phase_0 = [seq[i] for i in range(0, len(seq)) if i % 3 == 0]
        phase_1 = [seq[i] for i in range(0, len(seq)) if i % 3 == 1]
        phase_2 = [seq[i] for i in range(0, len(seq)) if i % 3 == 2]
        
        A_position = max(phase_0.count('A'), phase_1.count('A'), phase_2.count('A')) / (min(phase_0.count('A'), phase_1.count('A'), phase_2.count('A')) + 1.0)

        C_position = max(phase_0.count('C'), phase_1.count('C'), phase_2.count('C')) / (min(phase_0.count('C'), phase_1.count('C'), phase_2.count('C')) + 1.0)

        G_position = max(phase_0.count('G'), phase_1.count('G'), phase_2.count('G')) / (min(phase_0.count('G'), phase_1.count('G'), phase_2.count('G')) + 1.0)

        if self.read_type == 'DNA':
            T_position = max(phase_0.count('T'), phase_1.count('T'), phase_2.count('T')) / (min(phase_0.count('T'), phase_1.count('T'), phase_2.count('T')) + 1.0)
        elif self.read_type == 'RNA':
            T_position = max(phase_0.count('U'), phase_1.count('U'), phase_2.count('U')) / (min(phase_0.count('U'), phase_1.count('U'), phase_2.count('U')) + 1.0)

        fickett_score += self.look_up_content_prob(A_content, 'A', content_para, content_prob, content_weight)
        fickett_score += self.look_up_content_prob(C_content, 'C', content_para, content_prob, content_weight)
        fickett_score += self.look_up_content_prob(G_content, 'G', content_para, content_prob, content_weight)
        fickett_score += self.look_up_content_prob(T_content, 'T', content_para, content_prob, content_weight)
        
        fickett_score += self.look_up_position_prob(A_position, 'A', position_para, position_prob, position_weight)
        fickett_score += self.look_up_position_prob(C_position, 'C', position_para, position_prob, position_weight)
        fickett_score += self.look_up_position_prob(G_position, 'G', position_para, position_prob, position_weight)
        fickett_score += self.look_up_position_prob(T_position, 'T', position_para, position_prob, position_weight)  
        return fickett_score

    def fickett_value_full_sequence(self, seq):

        """calculate Fickett from full sequence - CPC2"""

        position_para = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0]
        content_para = [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19, 0.17, 0]

        position_prob = {
            'A': [0.51, 0.55, 0.57, 0.52, 0.48, 0.58, 0.57, 0.54, 0.50, 0.36],
            'C': [0.29, 0.44, 0.55, 0.49, 0.52, 0.60, 0.60, 0.56, 0.51, 0.38],
            'G': [0.62, 0.67, 0.74, 0.65, 0.61, 0.62, 0.52, 0.41, 0.31, 0.17],
            'T': [0.51, 0.60, 0.69, 0.64, 0.62, 0.67, 0.58, 0.48, 0.39, 0.24]}
        
        position_weight = {'A': 0.062, 'C': 0.093, 'G': 0.205, 'T': 0.154}
        content_weight = {'A': 0.084, 'C': 0.076, 'G': 0.081, 'T': 0.055}

        content_prob = {
            'A': [0.40, 0.55, 0.58, 0.58, 0.52, 0.48, 0.45, 0.45, 0.38, 0.19],
            'C': [0.50, 0.63, 0.59, 0.50, 0.46, 0.45, 0.47, 0.56, 0.59, 0.33],
            'G': [0.21, 0.40, 0.47, 0.50, 0.52, 0.56, 0.57, 0.52, 0.44, 0.23],
            'T': [0.30, 0.49, 0.56, 0.53, 0.48, 0.48, 0.52, 0.57, 0.60, 0.51]}

        if len(seq) < 2:
            return 0

        fickett_score = 0
        seq = seq.upper()
        total_base = len(seq)

        phase_0 = seq[::3]
        phase_1 = seq[1::3]
        phase_2 = seq[2::3]

        phase_0_A = phase_0.count('A')
        phase_1_A = phase_1.count('A')
        phase_2_A = phase_2.count('A')
        phase_0_C = phase_0.count('C')
        phase_1_C = phase_1.count('C')
        phase_2_C = phase_2.count('C')
        phase_0_G = phase_0.count('G')
        phase_1_G = phase_1.count('G')
        phase_2_G = phase_2.count('G')
        if self.read_type == 'DNA':
            phase_0_T = phase_0.count('T')
            phase_1_T = phase_1.count('T')
            phase_2_T = phase_2.count('T')
        elif self.read_type == 'RNA':
            phase_0_T = phase_0.count('U')
            phase_1_T = phase_1.count('U')
            phase_2_T = phase_2.count('U')

        A_content = float(phase_0_A + phase_1_A + phase_2_A) / total_base
        C_content = float(phase_0_C + phase_1_C + phase_2_C) / total_base
        G_content = float(phase_0_G + phase_1_G + phase_2_G) / total_base
        T_content = float(phase_0_T + phase_1_T + phase_2_T) / total_base
        A_position = max([phase_0_A, phase_1_A, phase_2_A]) / (min([phase_0_A, phase_1_A, phase_2_A]) + 1.0)
        C_position = max([phase_0_C, phase_1_C, phase_2_C]) / (min([phase_0_C, phase_1_C, phase_2_C]) + 1.0)
        G_position = max([phase_0_G, phase_1_G, phase_2_G]) / (min([phase_0_G, phase_1_G, phase_2_G]) + 1.0)
        T_position = max([phase_0_T, phase_1_T, phase_2_T]) / (min([phase_0_T, phase_1_T, phase_2_T]) + 1.0)

        fickett_score += self.look_up_content_prob(A_content, 'A', content_para, content_prob, content_weight)
        fickett_score += self.look_up_content_prob(C_content, 'C', content_para, content_prob, content_weight)
        fickett_score += self.look_up_content_prob(G_content, 'G', content_para, content_prob, content_weight)
        fickett_score += self.look_up_content_prob(T_content, 'T', content_para, content_prob, content_weight)

        fickett_score += self.look_up_position_prob(A_position, 'A', position_para, position_prob, position_weight)
        fickett_score += self.look_up_position_prob(C_position, 'C', position_para, position_prob, position_weight)
        fickett_score += self.look_up_position_prob(G_position, 'G', position_para, position_prob, position_weight)
        fickett_score += self.look_up_position_prob(T_position, 'T', position_para, position_prob, position_weight)

        return fickett_score

    def calculate_sequences(self):
        measure_orf = self.fickett_value_orf(self.read)
        measure_full = self.fickett_value_full_sequence(self.read)
        return [measure_orf, measure_full]

class EncodeScheme():
    def __init__(self, max_len, entropy_kmer = 3, entropy_q=2, kmer_freq = 3, read_type='DNA'):
        super(EncodeScheme, self).__init__()
        self.max_len = max_len
        self.entropy_kmer = entropy_kmer
        self.entropy_q = entropy_q
        self.kmer_freq = kmer_freq        
        self.read_type = read_type

    def atomic_number(self, read):
        mapping = []
        for nucle in read:
            if nucle == "T" or nucle == "U":
                mapping.append(66)
                # mapping.append(10)
            elif nucle == "C":
                mapping.append(58)
                # mapping.append(20)
            elif nucle == "A":
                mapping.append(70)
                # mapping.append(30)
            elif nucle == "G":
                mapping.append(78)
                # mapping.append(40)
            else:
                mapping.append(-1)
        padding = self.max_len + 1 - len(read)
        mapping = mapping + [0]*padding
        return mapping

    def binary(self, read):
        mapping = []
        for nucle in read:
            if nucle == "T" or nucle == "U":
                mapping.extend([1, 0, 0, 0])
                # mapping.append(10)
            elif nucle == "C":
                mapping.extend([0, 1, 0, 0])
                # mapping.append(20)
            elif nucle == "A":
                mapping.extend([0, 0, 1, 0])
                # mapping.append(30)
            elif nucle == "G":
                mapping.extend([0, 0, 0, 1])
                # mapping.append(40)
            else:
                mapping.extend([-1, -1, -1, -1])
        padding = self.max_len + 1 - len(read)
        mapping = mapping + [0, 0, 0, 0]*padding
        return mapping

    def descriptors(self, method, read):
        if method == "FourierTransform":
            features = FourierTransform(read).Encoding()
        elif method == "ChaosGame":
            features = ChaosGame(read, self.max_len, self.kmer_freq).Encoding()
        elif method == "Entropy":
            features = Entropy(read, self.entropy_kmer, self.entropy_q).Encoding()
        elif method == "FickettScore": # dim = 2
            features = FickettScore(read, self.read_type).calculate_sequences()
        elif method == "binary":
            features = self.binary(read)
        elif method == "atomic_number":
            features = self.atomic_number(read)
        return features