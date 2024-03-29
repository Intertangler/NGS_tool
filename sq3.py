__author__ = 'ian'
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline

import datetime
import os as os

import os, io, random
import string
import numpy as np
from tqdm import tqdm

from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import panel as pn
import panel.widgets as pnw


from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
pn.extension()

import Levenshtein
# import psutil
# import random
# import ast
# from Bio.Blast import NCBIWWW





class Search_Protocol:
    #attributes of the pattern are its name and coordinates
    def __init__(self,name):
        # ipdb.set_trace()
        self.name = name
        self.super_directory = "output"#+ str(start_t)#str(today) + str(self.name) + str(start_t)
        if not os.path.exists(self.super_directory):
            os.makedirs(self.super_directory)
        print( 'analysis results will be saved to: ', self.super_directory )

    def configure_score_matrix(self):
        self.score_matrix = {
        ('A','A'):4,('A','C'):-4,('A','G'):-4,('A','T'):-4,('A','R'):3,('A','Y'):-4,('A','S'):-4,('A','W'):3,('A','K'):-4,('A','M'):3,('A','B'):-4,('A','D'):2,('A','H'):2,('A','V'):2,('A','N'):0,('A','-'):0,
        ('C','A'):-4,('C','C'):4,('C','G'):-4,('C','T'):-4,('C','R'):-4,('C','Y'):3,('C','S'):3,('C','W'):-4,('C','K'):-4,('C','M'):3,('C','B'):2,('C','D'):-4,('C','H'):2,('C','V'):2,('C','N'):0,('C','-'):0,
        ('G','A'):-4,('G','C'):-4,('G','G'):4,('G','T'):-4,('G','R'):3,('G','Y'):-4,('G','S'):3,('G','W'):-4,('G','K'):3,('G','M'):-4,('G','B'):2,('G','D'):2,('G','H'):-4,('G','V'):2,('G','N'):0,('G','-'):0,
        ('T','A'):-4,('T','C'):-4,('T','G'):-4,('T','T'):4,('T','R'):-4,('T','Y'):3,('T','S'):-4,('T','W'):3,('T','K'):3,('T','M'):-4,('T','B'):2,('T','D'):2,('T','H'):2,('T','V'):-4,('T','N'):0,('T','-'):0,
        ('R','A'):3,('R','C'):-4,('R','G'):3,('R','T'):-4,('R','R'):3,('R','Y'):-4,('R','S'):0,('R','W'):0,('R','K'):0,('R','M'):0,('R','B'):-1,('R','D'):1,('R','H'):-1,('R','V'):1,('R','N'):0,('R','-'):0,
        ('Y','A'):-4,('Y','C'):3,('Y','G'):-4,('Y','T'):3,('Y','R'):-4,('Y','Y'):3,('Y','S'):0,('Y','W'):0,('Y','K'):0,('Y','M'):0,('Y','B'):1,('Y','D'):-1,('Y','H'):1,('Y','V'):-1,('Y','N'):0,('Y','-'):0,
        ('S','A'):-4,('S','C'):3,('S','G'):3,('S','T'):-4,('S','R'):0,('S','Y'):0,('S','S'):3,('S','W'):-4,('S','K'):0,('S','M'):0,('S','B'):1,('S','D'):-1,('S','H'):-1,('S','V'):1,('S','N'):0,('S','-'):0,
        ('W','A'):3,('W','C'):-4,('W','G'):-4,('W','T'):3,('W','R'):0,('W','Y'):0,('W','S'):-4,('W','W'):3,('W','K'):0,('W','M'):0,('W','B'):-1,('W','D'):1,('W','H'):1,('W','V'):-1,('W','N'):0,('W','-'):0,
        ('K','A'):-4,('K','C'):-4,('K','G'):3,('K','T'):3,('K','R'):0,('K','Y'):0,('K','S'):0,('K','W'):0,('K','K'):3,('K','M'):-4,('K','B'):1,('K','D'):1,('K','H'):-1,('K','V'):-1,('K','N'):0,('K','-'):0,
        ('M','A'):3,('M','C'):3,('M','G'):-4,('M','T'):-4,('M','R'):0,('M','Y'):0,('M','S'):0,('M','W'):0,('M','K'):-4,('M','M'):3,('M','B'):-1,('M','D'):-1,('M','H'):1,('M','V'):1,('M','N'):0,('M','-'):0,
        ('B','A'):-4,('B','C'):2,('B','G'):2,('B','T'):2,('B','R'):-1,('B','Y'):1,('B','S'):1,('B','W'):-1,('B','K'):1,('B','M'):-1,('B','B'):1,('B','D'):.5,('B','H'):.5,('B','V'):.5,('B','N'):0,('B','-'):0,
        ('D','A'):2,('D','C'):-4,('D','G'):2,('D','T'):2,('D','R'):1,('D','Y'):-1,('D','S'):-1,('D','W'):1,('D','K'):1,('D','M'):-1,('D','B'):.5,('D','D'):1,('D','H'):.5,('D','V'):.5,('D','N'):0,('D','-'):0,
        ('H','A'):2,('H','C'):2,('H','G'):-4,('H','T'):2,('H','R'):-1,('H','Y'):1,('H','S'):-1,('H','W'):1,('H','K'):-1,('H','M'):1,('H','B'):.5,('H','D'):.5,('H','H'):1,('H','V'):.5,('H','N'):0,('H','-'):0,
        ('V','A'):2,('V','C'):2,('V','G'):2,('V','T'):-4,('V','R'):1,('V','Y'):-1,('V','S'):1,('V','W'):-1,('V','K'):-1,('V','M'):1,('V','B'):.5,('V','D'):.5,('V','H'):.5,('V','V'):1,('V','N'):0,('V','-'):0,
        ('N','A'):0,('N','C'):0,('N','G'):0,('N','T'):0,('N','R'):0,('N','Y'):0,('N','S'):0,('N','W'):0,('N','K'):0,('N','M'):0,('N','B'):0,('N','D'):0,('N','H'):0,('N','V'):0,('N','N'):0,('N','-'):0,
        ('-','A'):0,('-','C'):0,('-','G'):0,('-','T'):0,('-','R'):0,('-','Y'):0,('-','S'):0,('-','W'):0,('-','K'):0,('-','M'):0,('-','B'):0,('-','D'):0,('-','H'):0,('-','V'):0,('-','N'):0,('-','-'):0
        }
        # print( 'Score Matrix: ', self.score_matrix)

    def get_fastq_files(self):
        self.file_list = []
        count_files = 0
        self.file_counter = 0
        for file in os.listdir(os.getcwd()):
            if file[-6:] == '.fastq':
                self.file_list.append(os.getcwd() + '/' + file)
                count_files += 1
        print('number of files to process: ', count_files)

    def get_single_reference_sequence(self):
        with open('reference_sequence.txt', 'r') as myfile:
            self.reference=myfile.read().replace('\n', '')

    def get_multiple_reference_sequences(self):
        self.reference_sequences = []
        with open("reference_sequences.fsa", "r") as file:
            for seq_record in SeqIO.parse(file, "fasta"):
                self.reference_sequences.append(seq_record)

    def initialize_user_parameters(self):
        default_threshold = 0.7
        self.cs_thresh = input("Please enter reference sequence mismatch tolerance (btw 0 and 1) or leave blank (default = " +str(default_threshold)+"): ") or default_threshold
    def initialize_fixed_parameters(self, cs_thresh= 0.75):
        self.cs_thresh = cs_thresh

    def initialize_file_data(self):
        self.reference_positions = []
        self.discarded_cs_position = []
        self.empty_tally = 0
        good_tally = 0
        self.discard_tally = 0
        self.full_tally = 0
        discard_by_barcode_problem_tally = 0
        self.reads_processed_counter = 0
        self.file_counter += 1
        self.directory_name = self.super_directory + '/' + str(self.file_counter)
        # os.makedirs(self.directory_name)

    def structure_data(self, sample_size, dimensions, landmarks):
        apriori_dimensions = len(self.reference_sequences)
        f_i = 0
        self.primary_sample = [SeqRecord(seq=Seq('AAAAAAAAAAAAAAAAAAAA'),
                                        id='artificial ctrl sequence - not data',
                                        name='artificial ctrl sequence - not data',
                                        description="included as a control for the analysis")]
        while len(self.primary_sample) < sample_size:
            file = self.file_list[f_i]
            self.primary_sample = self.primary_sample + [i for i in SeqIO.parse(file, "fastq")][:sample_size-len(self.primary_sample)]
            if len(self.primary_sample) >= sample_size:
                break
            f_i += 1
        self.distance_matrix = np.zeros((dimensions + apriori_dimensions, sample_size))
        self.initialize_file_data()
        # random_landmarks = np.random.randint(0, len(self.primary_sample), size=dimensions)
        for j in tqdm(range(0, sample_size)):
            # if j % 100 == 0: print(f'reads processed: {j}')
            for i in range(0, dimensions):
                self.distance_matrix[i, j] = Levenshtein.distance(str(self.primary_sample[landmarks[i]].seq),
                                                                  str(self.primary_sample[j].seq))
            for i in range(0, apriori_dimensions):
                self.distance_matrix[dimensions + i, j] = Levenshtein.distance(str(self.reference_sequences[i].seq),
                                                                          str(self.primary_sample[j].seq))
    def structure_data_single(self,sequence, dimensions,landmark_sequences):
        apriori_dimensions = len(self.reference_sequences)
        distance_array = np.zeros((dimensions + apriori_dimensions))
        for i in range(0, dimensions):
            distance_array[i] = Levenshtein.distance(str(landmark_sequences[i].seq),
                                                              str(sequence.seq))
        for i in range(0, apriori_dimensions):
            distance_array[dimensions + i] = Levenshtein.distance(str(self.reference_sequences[i].seq),
                                                                      str(sequence.seq))
        return distance_array


    def get_principal_components(self,ncomp,dim=4):
        from sklearn.decomposition import PCA
        from scipy.stats import gaussian_kde

        self.PCs = PCA(n_components=ncomp)
        self.PCs.fit(self.distance_matrix)
        print(f'Shape of the reduced data {self.PCs.components_.shape}')
        fig, ax = plt.subplots(nrows=int(ncomp / dim / 2), ncols=dim, figsize=(24, 12))
        y = -1
        x = 0
        for k in range(int(ncomp / 2)):
            x = k % dim
            if x == 0:  ## True in first round
                y += 1
            xy = np.vstack([self.PCs.components_[(2 * k), :], self.PCs.components_[(2 * k + 1), :]])
            z = gaussian_kde(xy)(xy)
            ax[y, x].scatter(self.PCs.components_[(2 * k), :], self.PCs.components_[(2 * k + 1), :], s=5, c=z, edgecolor='', cmap="cividis")
            ax[y, x].set_title(f'Component {2 * k} vs Component {2 * k + 1}')
        plt.savefig(self.name + "_principal_components.svg")
        plt.show()
        plt.close()

    def generate_3d_TSNE(self):
        from sklearn.manifold import TSNE
        from mpl_toolkits.mplot3d import Axes3D
        self.TSNE3 = TSNE(n_components=3, init='pca')
        self.TSNE3.fit(self.PCs.components_.T)
        print(f'TSNE embedding has shape: {self.TSNE3.embedding_.shape}')
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.TSNE3.embedding_[:, 0], self.TSNE3.embedding_[:, 1], self.TSNE3.embedding_[:, 2], alpha=.5, c="k", s=1)
        ax.set_title('tSNE 3d')
        plt.savefig(self.name+"_tsne3.svg")
        plt.close()

    def generate_2d_TSNE(self):
        from sklearn.manifold import TSNE
        self.TSNE2 = TSNE(n_components=2,
                     init='pca')  ## TSNE uses a random seed to initiate, meaning that the results don't always look the same!
        self.TSNE2.fit(self.PCs.components_.T)
        print(f'TSNE embedding has shape: {self.TSNE2.embedding_.shape}')
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111)
        ax.scatter(self.TSNE2.embedding_[:, 0], self.TSNE2.embedding_[:, 1], alpha=.5, c="k", s=2)
        ax.set_title('tSNE 2d')
        plt.savefig(self.name + "_tsne2.svg")
        plt.close()

    def cluster_agglomerative(self,distance_thresh = 0.5):
        from sklearn.cluster import AgglomerativeClustering
        self.AC = AgglomerativeClustering(distance_threshold=distance_thresh, n_clusters=None)
        self.AC = self.AC.fit(self.PCs.components_.T)
        self.nclusts = len(np.unique(self.AC.labels_))
        print(f'Identified {self.nclusts} clusters')
        from mpl_toolkits.mplot3d import Axes3D
        try:
            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(self.TSNE3.embedding_[:, 0], self.TSNE3.embedding_[:, 1], self.TSNE3.embedding_[:, 2], c=self.AC.labels_, s=3,
                       cmap='Spectral', alpha=.9)
            plt.savefig(self.name + "agglomerative_clust3.svg")
            plt.show()
            plt.close()
        except:pass
        try:
            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(111)
            ax.scatter(self.TSNE2.embedding_[:, 0], self.TSNE2.embedding_[:, 1], c=self.AC.labels_, s=10, cmap='Spectral', alpha=.9)
            ax.set_title('Agglomerative clustering')
            plt.savefig(self.name + "agglomerative_clust2.svg")
            plt.show()
            plt.close()
        except:pass

    def cluster_spectral(self):
        from sklearn.cluster import SpectralClustering

        self.spectral = SpectralClustering(n_clusters=self.nclusts)
        self.spectral.fit(self.PCs.components_.T)

        from mpl_toolkits.mplot3d import Axes3D
        try:
            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(self.TSNE3.embedding_[:, 0], self.TSNE3.embedding_[:, 1], self.TSNE3.embedding_[:, 2], c=self.spectral.labels_, s=3,
                       cmap='Spectral', alpha=.9)
            ax.set_title('Spectral clustering')
            plt.savefig(self.name + "spectral_clust3.svg")
            plt.show()
            plt.close()
        except:pass
        try:
            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(111)
            ax.scatter(self.TSNE2.embedding_[:, 0], self.TSNE2.embedding_[:, 1], c=self.spectral.labels_, s=10, cmap='Spectral', alpha=.9)
            ax.set_title('Spectral clustering')
            plt.savefig(self.name + "spectral_clust2.svg")
            plt.show()
            plt.close()
        except:pass

    def initialize_sorting_directories(self, nclusts):
        for clust in range(0, nclusts):
            if not os.path.exists(self.super_directory + '/' + str(clust)):
                os.makedirs(self.super_directory + '/' + str(clust))


class Alignment():
    def __init__(self, parent_search_protocol):
        self.psp = parent_search_protocol
        # invoking the __init__ of the parent class
        # Search_Protocol.__init__(self,name)
    def initialize_reference_threshold(self, reference):
        self.max_score = pairwise2.align.localds(reference.seq, reference.seq,
                                                                           self.psp.score_matrix, -15, -1,
                                                                           penalize_end_gaps=False,
                                                                           one_alignment_only=True)[0][2]
        self.reference_score_threshold = int(float(self.psp.cs_thresh) * self.max_score)
        self.reference_id = reference.id

    def align_read_to_ref(self, seq_record, reference):
        # ipdb.set_trace()
        self.psp.full_tally += 1
        #                 iteration_start_t = time.time()
        self.psp.reads_processed_counter += 1
        # print 'sequence ' + seq_record.id + ' length: ' + str(len(seq_record.seq))
        #                 #locate hybridization region by pw alignment with reference sequence
        plus_strand_alignment = pairwise2.align.localds(reference, seq_record.seq[:], self.psp.score_matrix, -15, -1,
                                                        penalize_end_gaps=False, one_alignment_only=True)
        minus_strand_alignment = pairwise2.align.localds(reference, seq_record.seq.reverse_complement()[:],
                                                         self.psp.score_matrix, -15, -1, penalize_end_gaps=False,
                                                         one_alignment_only=True)
        if max([i[2] for i in plus_strand_alignment]) > max([i[2] for i in minus_strand_alignment]):
            #                     #then we have a plus strand
            self.main_sequence = seq_record.seq
            main_alignments = plus_strand_alignment
        elif max([i[2] for i in plus_strand_alignment]) < max([i[2] for i in minus_strand_alignment]):
            #                     # then we have a minus strand
            self.main_sequence = seq_record.seq.reverse_complement()
            main_alignments = minus_strand_alignment
        else:
            #                     #then there is a tie, and this is strange
            if max([i[2] for i in plus_strand_alignment]) > self.reference_score_threshold:
                pass
                # print( 'error - plus and minus strands have equal hybridization region scores of ' + str(
                #     max([i[2] for i in plus_strand_alignment])))
            self.main_sequence = seq_record.seq
            main_alignments = plus_strand_alignment
        #                 #now we have a winning read that should be oriented correctly, if its scores isnt good though, we discard
        self.oriented_reference_alignment = main_alignments[0]

    def align_read_to_ref_unisense(self, seq_record, reference):
        self.psp.full_tally += 1
        self.psp.reads_processed_counter += 1
        plus_strand_alignment = pairwise2.align.localds(reference, seq_record.seq[:], self.psp.score_matrix, -15, -1,
                                                        penalize_end_gaps=False, one_alignment_only=True)
        self.oriented_reference_alignment = plus_strand_alignment[0]

    def discard_by_reference_fail(self):
        flipped_record_full = SeqRecord(self.main_sequence,
                                        id='flipped' + seq_record.id,
                                        name='flipped' + seq_record.id,
                                        description="flipped raw reads - both bad and good scoring")

        # SeqIO.write(flipped_record_full, e, "fasta")
        self.psp.discarded_cs_position += [self.oriented_reference_alignment[3] + 11]
        self.psp.reference_positions += [self.oriented_reference_alignment[3] + 11]
        discard_record_full = SeqRecord(self.main_sequence,
                                        id='_discard_bycommonsequence' + seq_record.id,
                                        name='_discard_bycommonsequence' + seq_record.id,
                                        description="reference alignment below threshold")
        self.psp.discard_tally += 1
        if not discard_record_full.seq:
            self.psp.empty_tally += 1
        else:
            # SeqIO.write(discard_record, b, "fasta")
            # SeqIO.write(discard_record, bb, "fasta")
            # SeqIO.write(discard_record_full, c, "fasta")
            # SeqIO.write(discard_record_full, cc, "fasta")
            # discard_records += [discard_record]
            if self.psp.reads_processed_counter % 1000 == 0:
                # iteration_time = iteration_end_t - iteration_start_t
                # average_iteration_time = (average_iteration_time * counter - 1 + iteration_time) / counter
                print( 'reads processed so far: ' + str(self.psp.reads_processed_counter))
                # print average_iteration_time
                # process = psutil.Process(os.getpid())
                # print(process.memory_info().rss)
            # ipdb.set_trace()
        self.discarded_record_entry = '-,-,discarded because reference alignment is below threshold'

    #
    #                 #the winning orientation has a high enough common sequence alignment score to permit barcode check

    def check_reference_alignment(self):
        if self.oriented_reference_alignment[2] < self.reference_score_threshold:
            # self.discard_by_reference_fail()
            proceed = False
        elif self.oriented_reference_alignment[2] >= self.reference_score_threshold:
            proceed = True
        else:
            proceed = False
        return proceed

    def extract_5p_upstream(self):
        # ipdb.set_trace()
        print('hit! score:' + str(self.oriented_reference_alignment[2]) + " versus threshold: " + str(self.reference_score_threshold) )
        self.five_prime_extract = SeqRecord(self.main_sequence[0:self.oriented_reference_alignment[3]+10],
                                        id='flipped' + seq_record.id,
                                        name='flipped' + seq_record.id,
                                        description="Strand-oriented read. Harvest because alignment score of " + str(self.oriented_reference_alignment[2]) + " exceeded the threshold " + str(self.reference_score_threshold) + ". Reference "+str(self.reference_id)+" aligns with read positions " + str(self.oriented_reference_alignment[3]) + " - " +str(self.oriented_reference_alignment[4]))

        # ipdb.set_trace()
        # harvested_alignments.write(str(flipped_record_full.id) + str(flipped_record_full.description) + "\n%s\n" % str(format_alignment(*self.oriented_reference_alignment)))
        # ipdb.set_trace()
        # umi_candidate = ''

    def check_5p_flanking_alignment(self,seq_record, reference):
        if self.oriented_reference_alignment[2] < self.reference_score_threshold:
            self.discard_by_reference_fail()
            proceed = False
        elif self.oriented_reference_alignment[2] >= self.reference_score_threshold:
            proceed = True
        else:
            proceed = False
        return proceed

    def harvest_hit_basic(self, harvested_reads, harvested_alignments):
        # ipdb.set_trace()
        print('hit! score:' + str(self.oriented_reference_alignment[2]) + " versus threshold: " + str(self.reference_score_threshold) )
        flipped_record_full = SeqRecord(self.main_sequence,
                                        id='flipped' + seq_record.id,
                                        name='flipped' + seq_record.id,
                                        description="Strand-oriented read. Harvest because alignment score of " + str(self.oriented_reference_alignment[2]) + " exceeded the threshold " + str(self.reference_score_threshold) + ". Reference "+str(self.reference_id)+" aligns with read positions " + str(self.oriented_reference_alignment[3]) + " - " +str(self.oriented_reference_alignment[4]))
        # ipdb.set_trace()
        SeqIO.write(flipped_record_full, harvested_reads, "fasta")
        # ipdb.set_trace()
        harvested_alignments.write(str(flipped_record_full.id) + str(flipped_record_full.description) + "\n%s\n" % str(format_alignment(*self.oriented_reference_alignment)))
        # ipdb.set_trace()

    def harvest_hit_parallel(self, parallel_harvest_record, parallel_record):
        # ipdb.set_trace()
        # print('hit! score:' + str(self.oriented_reference_alignment[2]) + " versus threshold: " + str(self.reference_score_threshold) )
        parallel_entry = SeqRecord(parallel_record.seq,
                                        id='flipped' + parallel_record.id,
                                        name='flipped' + parallel_record.id,
                                        description="paired end harvest")
        # ipdb.set_trace()
        SeqIO.write(parallel_entry, parallel_harvest_record, "fasta")

        # harvested_alignments.write(str(flipped_record_full.id) + str(flipped_record_full.description) + "\n%s\n" % str(format_alignment(*self.oriented_reference_alignment)))
        # ipdb.set_trace()

    def harvest_hit_insert(self, harvested_reads, harvested_alignments):
        # ipdb.set_trace()
        print('hit! score:' + str(self.oriented_reference_alignment[2]) + " versus threshold: " + str(self.reference_score_threshold) )
        flipped_record_full = SeqRecord(self.main_sequence,
                                        id='flipped' + seq_record.id,
                                        name='flipped' + seq_record.id,
                                        description="Strand-oriented read. Harvest because alignment score of " + str(self.oriented_reference_alignment[2]) + " exceeded the threshold " + str(self.reference_score_threshold) + ". Reference "+str(self.reference_id)+" aligns with read positions " + str(self.oriented_reference_alignment[3]) + " - " +str(self.oriented_reference_alignment[4]))
        SeqIO.write(flipped_record_full, harvested_reads, "fasta")
        # ipdb.set_trace()
        harvested_alignments.write(str(flipped_record_full.id) + str(flipped_record_full.description) + "\n%s\n" % str(format_alignment(*self.oriented_reference_alignment)))
        # ipdb.set_trace()

    def harvest_read(self, harvested_reads):
        # ipdb.set_trace()
        # print('hit! score:' + str(self.oriented_reference_alignment[2]) + " versus threshold: " + str(self.reference_score_threshold) )
        flipped_record_full = SeqRecord(self.main_sequence,
                                        id='flipped' + seq_record.id,
                                        name='flipped' + seq_record.id,
                                        description="Strand-oriented read. Harvest because alignment score of " + str(self.oriented_reference_alignment[2]) + " exceeded the threshold " + str(self.reference_score_threshold) + ". Reference "+str(self.reference_id)+" aligns with read positions " + str(self.oriented_reference_alignment[3]) + " - " +str(self.oriented_reference_alignment[4]))
        SeqIO.write(flipped_record_full, harvested_reads, "fasta")
        # ipdb.set_trace()
        # harvested_alignments.write(str(flipped_record_full.id) + str(flipped_record_full.description) + "\n%s\n" % str(format_alignment(*self.oriented_reference_alignment)))
        # ipdb.set_trace()
    def harvest_alignment(self, harvested_alignments):
        # ipdb.set_trace()
        print('hit! score:' + str(self.oriented_reference_alignment[2]) + " versus threshold: " + str(self.reference_score_threshold) )
        flipped_record_full = SeqRecord(self.main_sequence,
                                        id='flipped' + seq_record.id,
                                        name='flipped' + seq_record.id,
                                        description="Strand-oriented read. Harvest because alignment score of " + str(self.oriented_reference_alignment[2]) + " exceeded the threshold " + str(self.reference_score_threshold) + ". Reference "+str(self.reference_id)+" aligns with read positions " + str(self.oriented_reference_alignment[3]) + " - " +str(self.oriented_reference_alignment[4]))
        # SeqIO.write(flipped_record_full, harvested_reads, "fasta")
        # ipdb.set_trace()
        harvested_alignments.write(str(flipped_record_full.id) + str(flipped_record_full.description) + "\n%s\n" % str(format_alignment(*self.oriented_reference_alignment)))
        # ipdb.set_trace()

    def harvest_insert(self, harvested_inserts, insert_sequence):
        # ipdb.set_trace()
        print('INSERT HARVESTED!')
        # insert_record = SeqRecord(insert_sequence,
        #                                 id='flipped' + seq_record.id,
        #                                 name='flipped' + seq_record.id,
        #                                 description="")
        # ipdb.set_trace()
        SeqIO.write(insert_sequence, harvested_inserts, "fasta")
        # ipdb.set_trace()
        # harvested_alignments.write(str(insert_record.id) + str(insert_record.description) + "\n%s\n" % str(format_alignment(*self.oriented_reference_alignment)))
        # ipdb.set_trace()

def make_seq(length=40):
    return ''.join([random.choice(['A','C','T','G']) for i in range(length)])

def get_colors(seqs):
    """make colors for bases in sequence"""
    text = [i for s in list(seqs) for i in s]
    clrs =  {'A':'red','T':'green','G':'orange','C':'blue','-':'white'}
    colors = [clrs[i] for i in text]
    return colors

def view_alignment(aln, fontsize="9pt", plot_width=800):
    """Bokeh sequence alignment view"""

    #make sequence and id lists from the aln object
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]
    text = [i for s in list(seqs) for i in s]
    colors = get_colors(seqs)
    N = len(seqs[0])
    S = len(seqs)
    width = .4

    x = np.arange(1,N+1)
    y = np.arange(0,S,1)
    #creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    #flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    #use recty for rect coords with an offset
    recty = gy+.5
    h= 1/S
    #now we can create the ColumnDataSource with all the arrays
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs)*15+50
    x_range = Range1d(0,N+1, bounds='auto')
    if N>100:
        viewlen=100
    else:
        viewlen=N
    #view_range is for the close up view
    view_range = (0,viewlen)
    tools="xpan, xwheel_zoom, reset, save"

    #entire sequence view (no text, with zoom)
    p = figure(title=None, plot_width= plot_width, plot_height=50,
               x_range=x_range, y_range=(0,S), tools=tools,
               min_border=0, toolbar_location='below')
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                 line_color=None, fill_alpha=0.6)
    p.add_glyph(source, rects)
    p.yaxis.visible = False
    p.grid.visible = False

    #sequence text view with ability to scroll along x axis
    p1 = figure(title=None, plot_width=plot_width, plot_height=plot_height,
                x_range=view_range, y_range=ids, tools="xpan,reset",
                min_border=0, toolbar_location='below')#, lod_factor=1)
    glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                text_font="monospace",text_font_size=fontsize)
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                line_color=None, fill_alpha=0.4)
    p1.add_glyph(source, glyph)
    p1.add_glyph(source, rects)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0

    p = gridplot([[p],[p1]], toolbar_location='below')
    return p






################ BLAST SIMPLE ##########################
#Code is borrowed based on that of Rocha and Ferreira 2018

#first step is to create a database
#the database is a simple list of sequences (strings)
#we define a function that loads the list from a text file
def read_database(filename):
    f = open(filename)
    db = []
    for line in f:
        db.append(line.rstrip())
    f.close()
    return db

#next function creates a dictionary - a map that will identify all words of size w
#we keep track of the positions where it occurs
#this technique is known as hashing
def build_map(query, w, subseqprint = False):
    res = {}
    for i in range(len(query) - w + 1):
        subseq = query[i:i+w]
        if subseq in res:
            res[subseq].append(i)
        else:
            res[subseq] = [i]
        if subseqprint == True: print(subseq)
    return res

#in our simplified version of blast, we will not use a substitution matrix
#we will not consider gaps
#only scores of 1 or 0, match of mismatch will be considered
def get_hits(seq,m,w):
    res = []
    for i in range(len(seq)-w+1):
        subseq = seq[i:i+w]
        if subseq in m:
            l = m[subseq]
            for ind in l:
                res.append((ind,i))
    return res


def extends_hit(seq, hit, query, w):
    stq, sts = hit[0], hit[1]
    ## move forward
    matfw = 0
    k = 0
    bestk = 0
    while 2 * matfw >= k and stq + w + k < len(query) and sts + w + k < len(seq):
        if query[stq + w + k] == seq[sts + w + k] or query[stq + w + k] == "N"   :
            matfw += 1
            bestk = k + 1
        k += 1
    size = w + bestk
    ## move backwards
    k = 0
    matbw = 0
    bestk = 0
    while 2 * matbw >= k and stq > k and sts > k:
        if query[stq - k - 1] == seq[sts - k - 1] or query[stq - k - 1] == "N"  :
            matbw += 1
            bestk = k + 1
        k += 1
    size += bestk
    # if sts - bestk == 2946: import pdb; pdb.set_trace()
    return (stq - bestk, sts - bestk, size, w + matfw + matbw)

def extends_hitbackup(seq, hit, query, w):
    stq, sts = hit[0], hit[1]
    ## move forward
    matfw = 0
    k = 0
    bestk = 0
    while 2 * matfw >= k and stq + w + k < len(query) and sts + w + k < len(seq):
        if query[stq + w + k] == seq[sts + w + k]:
            matfw += 1
            bestk = k + 1
        k += 1
    size = w + bestk
    ## move backwards
    k = 0
    matbw = 0
    bestk = 0
    while 2 * matbw >= k and stq > k and sts > k:
        if query[stq - k - 1] == seq[sts - k - 1]:
            matbw += 1
            bestk = k + 1
        k += 1
    size += bestk
    # if sts - bestk == 2946: import pdb; pdb.set_trace()
    return (stq - bestk, sts - bestk, size, w + matfw + matbw)

def hit_best_score(seq,query,m,w):
    hits = get_hits(seq,m,w)
    bestScore = -1.0
    best = ()
    for h in hits:
        ext = extends_hit(seq,h,query,w)
        score = ext[3]
        # import pdb; pdb.set_trace()
        if score > bestScore or (score == bestScore and ext[2] < best[2]):
            bestScore = score
            best = ext


    return best
def best_alignment(db,query,w,subseqprint = True):
    m = build_map(query,w,subseqprint)
    bestScore = -1.0
    res = (0,0,0,0,0)
    for k in range(0,len(db)):
        bestSeq = hit_best_score(db[k],query,m,w)
        if bestSeq != ():
            score = bestSeq[3]
            if score > bestScore or (score == bestScore and bestSeq[2] <res[2]):

                bestScore = score
                res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3], k
    if bestScore < 0: return ()
    else: return res

def threshold_alignments(db,db_ids, query,w,threshold = 0.5, subseqprint = False):
    threshold_score = len(query)*threshold
    m = build_map(query,w,subseqprint)
    res = []
    for k in range(0,len(db)):
        bestSeq = hit_best_score(db[k],query,m,w)
        if bestSeq != ():
            score = bestSeq[3]
            if score >= threshold_score:
                rec = SeqRecord(Seq(db[k]),
                                                           id=db_ids[k],
                                                           name=db_ids[k],
                                                           description=f"blast hit with score: {score/len(query)}")
                res.append(rec)
                # res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3], k
    # if bestScore < 0: return ()
    if res:
        return res

def read_database_fasta(filename):
    from Bio import SeqIO
    db = []
    with open(filename, "r") as file:
        for seq_record in SeqIO.parse(file, "fasta"):
            db.append(seq_record)
    return db

def best_alignment_fasta(db,query,w, subseqprint = False):
    from Bio import SeqIO
    m = build_map(query,w,subseqprint)
    bestScore = -1.0
    res = []#(0,0,0,0,0)
    for k in range(0,len(db)):
        bestSeq = hit_best_score(str(db[k].seq),query,m,w)
        if bestSeq != ():
            score = bestSeq[3]
            if score > bestScore or (score == bestScore and bestSeq[2] <res[2]):
                bestScore = score
                res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3], k


    if bestScore < 0: return ()
    else: return res

def batch_iterator(iterator, batch_size):
    entry = True
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

def split_files(reads_per_file=5000, directory_name="ngs_data", source_directory=None):
    file_list = []
    if not source_directory:
        source_directory = os.getcwd()
    for file in os.listdir():
        if file[-6:] == '.fastq':
            file_list.append(os.getcwd() + '/' + file)
    if not os.path.exists(directory_name):
        os.mkdir(directory_name)
        print("Directory ", directory_name, " Created ")
    else:
        print("Directory ", directory_name, " already exists")

    # os.makedirs(directory_name)
    directory_path = str(os.getcwd()) + '/' + str(directory_name) + '/'
    # A = True

    filename = file_list[0]
    master_count = 0
    record_iter = SeqIO.parse(filename, "fastq")
    for i, batch in enumerate(batch_iterator(record_iter, reads_per_file)):
        filename = directory_path + "group_%i.fastq" % (i + 1)
        with open(filename, "w") as handle:
            count = SeqIO.write(batch, handle, "fastq")
        master_count += reads_per_file
        if master_count % 10000 == 0 or master_count % (reads_per_file * 10) == 0: print(
            "%i reads processed" % master_count)