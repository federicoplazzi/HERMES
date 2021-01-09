# coding=utf-8
#######################################################################################################################
#                                                                                                                     #
# HERMES is a straightforward index which tries to summarize the mitochondrial evolution pace using a single number.  #
# Several mitogenomic features are evaluated in a factor analysis framework; namely, in the current version, they are #
# the %URs, the Amount of Mitochondrial Identical Gene Arrangements (AMIGA) index, the absolute value of the SU skew, #
# the root-to-tip distance, the ML pairwise distance from a given outgroup, the %AT, the AT skew, the GC skew, the    #
# number of (annotated) genes, the length of the mtDNA, the CAI, the topology, the UR AT content, and the AT median   #
# length.                                                                                                             #
#                                                                                                                     #
# Copyright (C) 2021 Guglielmo Puccio, Federico Plazzi                                                                #
#                                                                                                                     #
# This program is free software: you can redistribute it and/or modify                                                #
# it under the terms of the GNU General Public License as published by                                                #
# the Free Software Foundation, either version 3 of the License, or	                                              #
# (at your option) any later version.                                                                                 #
#                                                                                                                     #
# This program is distributed in the hope that it will be useful,                                                     #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                                      #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                       #
# GNU General Public License for more details.                                                                        #
#                                                                                                                     #
# You should have received a copy of the GNU General Public License                                                   #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.                                               #
#                                                                                                                     #
#######################################################################################################################
import warnings
from statistics import median
from math import log, exp
from Bio import AlignIO
from Bio import Entrez
from Bio import BiopythonWarning
from sys import exit
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils import GC
from ete3 import Tree
import os, argparse, subprocess
import Bio.Data.CodonTable

Entrez.email = 'mozoolab@gmail.com'

warnings.simplefilter('ignore', BiopythonWarning)  # we silence biopython warnings

# arguments parsing
parser = argparse.ArgumentParser()
# change the optional title
parser._optionals.title = "Arguments"
# input file
parser.add_argument('-I', dest='gb_file', required=True,
                    help='Input File containing the records in Genebank format (.gb extension) or a list of Genbank accessions (.list extension)')
# list of tree names and NCBI ids
parser.add_argument('-D', dest='entry_names', required=True,
                    help='File listing all entry names with the corresponding NCBI ids')
# tree in Newick format
parser.add_argument('-L', dest='tree', required=True, help='Tree file in Newick format')
# outgroup in the tree
parser.add_argument('-O', dest='outgroup', required=False, help='Name of the outgroup')
# genetic code
parser.add_argument('-G', dest='code', required=False,
                    help='NCBI Genetic Code ID. Default is 5 (The Invertebrate Mitochondrial Code)', default="5")
# arguments for RAxML
parser.add_argument('-s', dest='alignment', required=True, help='Alignment file')
parser.add_argument('-m', dest='model', required=True, help='Molecular evolution model (for RAxML)')
parser.add_argument('-q', dest='partition', help='Partition file (for RAxML)')
parser.add_argument('-t', dest='threads', help='Number of CPU cores to be used (for RAxML multithreading)', default="1")
# alpha for RMSEA
parser.add_argument('-a', dest='alpha', type=float, help='Confidence level for RMSEA CI (for factor analysis)',
                    default="0.05")
# custom annotation comma separated
parser.add_argument('-A', dest='annotation', type=str,
                    help='Custom Annotation for all the genomes not available in genbank (comma separated)')
parser.add_argument('-F', dest='fasta', type=str,
                    help='Fasta file with the custom annotated sequences (only needed if a custom annotation is given)')
args = parser.parse_args()

output = "HERMES"
version = "3.0"
try:
    os.mkdir("./Results")
except:
    pass

file_for_R = open("./Results/HERMES_variables.txt", "w")

try:
    os.remove("RAxML_distances." + output)
except:
    pass
try:
    os.remove("RAxML_info." + output)
except:
    pass
try:
    os.remove("RAxML_parsimonyTree." + output)
except:
    pass


def append_gene(name, feature):  # it takes the NCBI name and the feature from genbank
    Q = feature.qualifiers
    global f, D
    if 'gene' in list(Q.keys()):
        if Q['gene'][0].lower() in list(D.keys()):  # here we check if the dictionary contain the gene
            if feature.strand == 1:
                f[name].append(["+", D[Q['gene'][0].lower()]])
            else:
                f[name].append(["-", D[Q['gene'][0].lower()]])
        else:  # if it isn't in our dictionary we quit and print the name and the instructions for the user
            exit(
                "The name of a gene in your annotation or genbank file is not clear, please change the name according to HERMES manual, or update the Genes.dict file in the HERMES folder (both steps covered inside the HERMES manual)\n The gene is the " + str(
                    Q['gene'][0].lower()))
    elif 'product' in list(Q.keys()):
        if Q['product'][0].lower() in list(D.keys()):  # here we check if the dictionary contain the product
            if feature.strand == 1:
                f[name].append(["+", D[Q['product'][0].lower()]])
            else:
                f[name].append(["-", D[Q['product'][0].lower()]])
        else:
            exit(
                "The name of a gene in your annotation or genbank file is not present in our dictionary, please change the name according to HERMES manual, or update the Genes.dict file in the HERMES folder (both steps covered inside the HERMES manual)\n The gene is the " + str(
                    Q['product'][0].lower()))
    else:  # if the feature isn't annotated as gene or product
        if feature.strand == 1:
            f[name].append(["+", "orfan"])
        else:
            f[name].append(["-", "orfan"])
    return


def change_sign(gene):
    if gene[0] == '-':
        gene[0] = '+'
    else:
        gene[0] = '-'
    return


# if a custom annotation is given, we check for the fasta (mandatory)
if args.annotation:
    if args.fasta:
        pass
    else:
        exit(
            "ERROR: If you provide a custom annotation file, you have to provide a FASTA file containing the annotated sequences with the -F option (check the manual for more informations)")

# GB entries download
# if the user gives a gb accessions list we directly download them
if args.gb_file.endswith(".gb"):  # if the user gives the genbank file we skip this part
    args.gb_file = open(args.gb_file, "r")  # CONTROLLA: SE FUNZIONA COSÌ ANCHE DANDO I GB COME FILES
    pass
elif args.gb_file.endswith(".list"):  # if the user provide a list with Genbank accessions
    with open(args.gb_file) as f:
        gb_list = ",".join(f.read().splitlines())
    args.gb_file = Entrez.efetch(db="nuccore", id=gb_list, rettype="gb",
                                 retmode="text")  # we directly download them and replace the args.gb_file with the genbank records

else:
    exit("Input file must have a .gb or .list extension. Check out the manual or the help for detailed info")

# CUSTOM_ANNOTATION_PARSING
# first we check if user provided a custom annotation and split the file

ANN = {}
if args.annotation:
    a = 0
    ann = open(args.annotation).readlines()
    for i in range(len(ann)):
        ann[i] = ann[i].strip().split(",")
        # then we extract the annotation for each genome that later will be added to the annotations extracted from genebank
        if len(ann[i]) == 2:
            ANN[ann[i][0]] = [ann[i][1], []]  # list of lists (2) that are genome_length and genes_annotation
            a = i
        else:
            ANN[ann[a][0]][1].append(ann[i])

# Here we create the general mitochondrial Dictionary from the file
DD = open("./Genes.dict").readlines()
D = {}
for i in range(len(DD)):
    DD[i] = DD[i].strip().split(",")
    D[DD[i][0]] = DD[i][1].upper()

# Here we create the dictionary to translate colors into R color index (colors table in HERMES folder)
col = {}  # dictionary with names for keys and numbers for values
colors = open("./colors").readlines()
for i in range(len(colors)):
    colors[i] = colors[i].split()
    col[colors[i][1]] = int(colors[i][0])

# Here we create the dictionary to translate the genbank ids with the names used in the tree (usually species)
NN = open(args.entry_names).readlines()
N = {}  # Dictionary with species as keys and NCBI IDs as values
M = {}  # Dictionary with NCBI IDs as keys and species as values (it is used in the end to pass the species name to the variables set)
O = {}  #dictionary for the color coded version of the script
P = {}  #dictionary with NCBI IDs as keys and nucleotide range as values

#these files are needed for the multiple chromosomes check
gb_file = SeqIO.to_dict(SeqIO.parse(args.gb_file, "genbank"))   #genbank file is parsed as a dictionary
sequence_string = ""
sequence_object = Seq(sequence_string)
length=0
length1=0
counter=0
delete_list = []    #list with ncbi id's to be deleted (multiple chromosomes except the first)
delete_spec = []    #the same as before but with species ids

#universal split for the NN file
for i in range(len(NN)):
    NN[i] = NN[i].strip().split(",")

#This chunk of code is for the multi-chromosome parsing. It checks if the user provided multiple chromosome (gb entry) using ":" in the names file and merge them together
for i in range(len(NN)):
    # this part of code is used to check if any genome has multiple chromosomes. If so the genbanks are merged together updating the annotation and the single chromosomes gb are deleted
    if ":" in NN[i][1]:  # given the multiple ncbi id's we rename (later) this entry by the first ncbi id creating a new gb record
        for record in list(gb_file.keys()):
            if gb_file[record].id.split(".")[0] in NN[i][1].split(":"):
                if counter == 0:  # if it's the first match we create a copy of the genbank entry
                    first_record = record  # we save the key of the first record to use it as the main to be updated. the other chromosomes will be deleted once concatenated
                for feat in gb_file[record].features:
                    if feat.type == "source":  # we store the length of the genome because is the first feature but we use it after the for
                        length1 = feat.location.end
                        if counter == 0:
                            counter = 1
                            break  # like this we don't add new features to the first gbrecord (it would create a double copy)
                    newstart = feat.location.start + length
                    newend = feat.location.end + length
                    new_feature_location = FeatureLocation(newstart, newend, strand=feat.location.strand)
                    feat.location = new_feature_location
                    gb_file[first_record].features.append(
                        feat)  # here we add the updated feature to the first chromosome
                length = length + length1  # here we sum the length of the previous genomes with ne one we just analysed
                sequence_string = sequence_string + str(gb_file[record].seq)
                sequence_object = Seq(sequence_string)
                delete_list.append(record)  # we save the key of the chromosomes to be deleted

                # print(len(new_record.seq),len(record.seq))
        gb_file[first_record].seq = sequence_object
        #NN[i][1] = NN[i][1].split(":")[0]  # we replace the names with ":" with a single entry (the first) OLD it created errors because the order in the dictionary (gb file) may
        #differ from that in the name file (in the : column).
        NN[i][1] = first_record.split(".")[0]   #here we replace the ncbi id's (multiple) with the first record
    sequence_string = ""
    length = 0

if counter == 1:    #we use counter to check if there was at least one multi-chromosome entry
    delete_list.remove(first_record)    #then we remove the first entry from the list
    for el in delete_list:  #and we delete all the other records (chromosomes) from the gb_file
        del(gb_file[el])



if len(NN[0]) > 2:
    for i in range(len(NN)):
        if NN[i][1] == "-":  # if the taxa has a custom annotation just use the name
            N[NN[i][0]] = NN[i][0]
            M[NN[i][0]] = NN[i][0]
            # here we check if the user provided a number or the name for a color
            if NN[i][2].isdigit():
                O[NN[i][0]] = int(NN[i][2])
            else:
                O[NN[i][0]] = col[NN[i][2]]
        else:
            N[NN[i][0]] = NN[i][1]  # Dictionary with species as keys and NCBI IDs as values
            M[NN[i][1]] = NN[i][0]  # Dictionary with NCBI IDs as keys and species as values
            if NN[i][2].isdigit():
                O[NN[i][1]] = int(NN[i][2])  # Dictionary with NCBI IDs as keys and color number as value
            else:
                O[NN[i][1]] = col[NN[i][2]]
else:
    for i in range(len(NN)):
        NN[i] = NN[i].strip().split(",")
        if NN[i][1] == "-":
            N[NN[i][0]] = NN[i][0]
            M[NN[i][0]] = NN[i][0]
        else:
            N[NN[i][0]] = NN[i][1]  # Dictionary with species as keys and NCBI IDs as values
            M[NN[i][1]] = NN[i][0]  # Dictionary with NCBI IDs as keys and species as values



# Here we check that entries in the alignment, entrieslist and in the tree are exactly the same
# First we extract the names for the three files
S_names = list(N.keys())  # species in the entries list

T_leaves = t = Tree(args.tree, format=1).get_leaves()  # species in the tree file
T_names = []  # species in the tree
for leaf in T_leaves:
    T_names.append(leaf.name)


# alignment = AlignIO.read(args.alignment, "phylip") #AlignIO strangely can't open the alignment, we have to obtain the id by splitting with empty spaces
alignment = open(args.alignment).readlines()

A_names = []  # species in the alignment file

for line in alignment[1:]:
    line = line.split(" ")[0]
    if line != "" and line != "\n":
        A_names.append(line)
    else:
        continue
# here we create sets to be able to compare them

A_names = set(A_names)
T_names = set(T_names)
S_names = set(S_names)

if A_names == T_names == S_names:
    pass
else:
    print(A_names, T_names, S_names)
    diff = (A_names | T_names | S_names) - (A_names & T_names & S_names)
    print(
        "FORMAT ERROR\nThe following species names differ in at least one of the files: Tree, Entry_names, Alignment\n")
    print("\n".join(diff))
    quit()  # i've removed the quit that could block the script for some mistakes that may not be critical

# -----------------------------------------------MAIN--------------------------------------------------#

f = {}  # it is used by the function append gene
l = {}  # it is used later for the AT% in URs
q = {}  # it is used later for the URs median length
g = {}  # it is used later for the URs percentage
e = {}  # it is used later for the Circularity
h = {}  # it is used later for the SU skew
s = {}  # here we store all the remaining variables (AT_content,AT_skew,GC_skew,Genes,Length,Cai_index)

# codons
CodonsDict = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
              'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
              'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
              'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
              'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
              'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
              'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
              'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
              'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
              'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
              'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
              'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
              'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

# Here we generate the dictionary using CodonTable
t = Bio.Data.CodonTable.unambiguous_dna_by_id[int(args.code)]
bt = dict()
for a1 in "ATCG":  # Here we generate all the possible codons and append them to the corresponding amminoacid entry in the dictionary
    for a2 in "ATCG":
        for a3 in "ATCG":
            codon = a1 + a2 + a3
            try:
                amino = t.forward_table[codon]
            except KeyError:
                assert codon in t.stop_codons
                continue
            try:
                bt[amino].append(codon)
            except KeyError:
                bt[amino] = [codon]
bt2 = dict(bt)
for amino in bt:  # here we split the groups of synonymous codons that has more than one family
    aminolist = []
    for i in bt[amino]:
        aminolist.append(i[:2])  # we use a list with only the two initial bases of each codon
    if len(set(aminolist)) == 1:  # if the set (uniq) is longer than 1 means that more than one family is present
        continue
    else:
        for i in range(0, len(set(aminolist))):
            bt2[amino + "_" + str(i + 1)] = []
            for j in bt[amino]:
                if j[:2] == list(set(aminolist))[i]:  # this works for any given number of (sub)families
                    bt2[amino + "_" + str(i + 1)].append(j)
        bt2.pop(amino, None)  # after splitting the aminoacid we delete the old key
# here we delete from the dict aminoacids with only one synonymous codon (they would have score of 1 but they would count in the denominator)
for i in bt:
    if len(bt[i]) == 1:
        print(i, bt[i])
        try:
            bt2.pop(i, None)
        except:
            continue
SynonymousCodons = dict(bt2)  # here i use the same old variable

# here is the old dictionry that we used. Now we automatically generate a dictionary based on the genetic code selected with the -G argument.

# SynonymousCodons = {
# 	'CYS': ['TGT', 'TGC'],
# 	'ASP': ['GAT', 'GAC'],
# 	'SER_1': ['TCT', 'TCG', 'TCA', 'TCC'],
# 	'SER_2': ['AGC', 'AGT','AGA','AGG'],
# 	'GLN': ['CAA', 'CAG'],
# 	'MET': ['ATG','ATA'],
# 	'ASN': ['AAC', 'AAT'],
# 	'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
# 	'LYS': ['AAG', 'AAA'],
# 	'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
# 	'PHE': ['TTT', 'TTC'],
# 	'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
# 	'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
# 	'ILE': ['ATC', 'ATT'],
# 	'LEU_1': ['TTA', 'TTG'],
# 	'LEU_2': ['CTC', 'CTT', 'CTG', 'CTA'],
# 	'HIS': ['CAT', 'CAC'],
# 	'ARG': ['CGA', 'CGC', 'CGG', 'CGT'],
# 	'TRP': ['TGG','TGA'],
# 	'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
# 	'GLU': ['GAG', 'GAA'],
# 	'TYR': ['TAT', 'TAC']
# }

# 'STOP': ['TAG', 'TAA'],

CodonsFreq = {}

# here we parse and evaluate everything for the genbank datasets

# all results will be stored in various dcitionaries composed by:
# key = genome (from genbank or from custom annotation)
# value = variable (one of the variables evaluated like CAI or AtSKEW)
# so all dictionaries shares the same keys but values differs so that in the end we can use the same key to obtain
# different variables from each dictionary

# for gb_record in SeqIO.parse(open(args.gb_file, "r"), "genbank"):
H=0
for gb_record in gb_file:
    H=H+1
    #for circularity
    if gb_file[gb_record].annotations["topology"] == "linear":   #linear chromosomes will be 1 and circular 0
        e[(gb_file[gb_record].id).split(".")[0]] = 1
    else:
        e[(gb_file[gb_record].id).split(".")[0]] = 0
    Cai_list = []
    CodonsDict_2 = dict(CodonsDict)  # we reset the codon count for each gb record
    # for AT content
    AT = 100 - GC(gb_file[gb_record].seq)
    # for AT SKEW
    As = int(gb_file[gb_record].seq.count("A"))
    Ts = int(gb_file[gb_record].seq.count("T"))
    ATskew = float(As - Ts) / (As + Ts)
    # for GC
    Gs = int(gb_file[gb_record].seq.count("G"))
    Cs = int(gb_file[gb_record].seq.count("C"))
    # for GC skew
    GCskew = float(Gs - Cs) / (Gs + Cs)
    # for genome length
    Glength = len(gb_file[gb_record].seq)
    F = gb_file[gb_record].features  # F now are the genbank features
    LL = len(gb_file[gb_record].seq)  # equal to Glength, but for logical purpose i will keep two different variables
    L = [1 for i in range(len(gb_file[gb_record].seq))]
    f[(gb_file[gb_record].id).split(".")[
        0]] = []  # here we initialize the dictionaries with the gb_file[gb_record].id without the version
    h[(gb_file[gb_record].id).split(".")[0]] = [[], []]  # h will have [ [abs_su_skew (evaluated later)] , [strand]]
    NG = 0
    for i in range(1, len(F)):  # for each feature in the gb_file[gb_record] starting from 1 (0 is the total genome)
        if F[i].type in ["CDS", "rRNA", "tRNA"]:
            # for number of genes
            NG = NG + 1
            # for SU skew
            h[(gb_file[gb_record].id).split(".")[0]][1].append(F[i].strand)
            # for URs
            #if gb_file[gb_record].id == "NC_018373.1":
                #print(F[i].location.parts)
            for x in F[i].location.parts:
                for y in x:
                    if L[y] == 1:
                        L[y] = 0  # all 0s will represent coding regions, instead 1 will be non coding

        # for CAI
        if F[i].type == "CDS":
            append_gene((gb_file[gb_record].id).split(".")[0], F[i])  # here we update the f list with strand and the gene name (translated with our dictionary)
            Seq = F[i].extract(gb_file[gb_record].seq)  # here we extract each feature's sequence
            for j in range(3, len(Seq), 3):  # we start from 3 (the fourth nucleotide) to skip the start codon
                # these exceptions WERE used to detect the truncated stop codons but now we actually use them to skip them and also the full stop codons TAG and TAA (old version of cai is commented)
                if Seq[j:j + 3] == "T":  # these exceptions are for the truncated mitochondrial stop codons
                    pass
                # CodonsDict_2["TAA"] = CodonsDict_2["TAA"] +1
                elif Seq[j:j + 3] == "TA":
                    pass
                # CodonsDict_2["TAA"] = CodonsDict_2["TAA"] +1
                # elif Seq[j:j+3] == "TAA" or Seq[j:j+3] == "TAG":	#we skip also the full stop codons
                #	pass
                # ACTUALLY i've commented this part because i removed the STOP codons from the table, so they should be included in the following elif statement
                elif Seq[j:j + 3] not in CodonsDict_2:  # if a codon is not in our table we skip it
                    continue
                else:
                    CodonsDict_2[Seq[j:j + 3]] = CodonsDict_2[Seq[j:j + 3]] + 1  # and we count the occurrence of each codon
    #for % AT in URs
    aa = 0  # number of As in URs
    tt = 0  # number of Ts in URs
    for i in range(len(L)):
        if L[i] == 1:
            if gb_file[gb_record].seq[i] == "A":
                aa = aa + 1
            elif gb_file[gb_record].seq[i] == "T":
                tt = tt + 1
    #and Median length of URs
    list_string = map(str, L)   #first we remap the integers to strings
    list_string_2 = "".join(list_string).split("0") #then i split the joined string to obtain empty strings instead of 0s
    list_string = list(filter(None, list_string_2)) #and finally i filter out the empty strings
    urs_len = []
    for i in list_string:  #to evaluate median length of URs i use split to remove the 0s (coding portions)
        urs_len.append(len(i))
    #print(gb_file[gb_record].id)
    #print(L)
    #print(urs_len)
    #print(median(urs_len))
    #print(sum(urs_len) / len(urs_len))
    #print(list_string_2)  #questo comando è bellissimo. Fa vedere proprio le regioni non codificanti come 1 e le codificanti come stringhe vuote!
    l[(gb_file[gb_record].id).split(".")[0]] = ((aa + tt) / float(sum(L))) * 100
    q[(gb_file[gb_record].id).split(".")[0]] = median(urs_len)


    # now for CAI we find the frequencies of each codon for a residue
    CAI_dict = {}
    for i in list(
            SynonymousCodons.keys()):  # HERMES1 version (commented) had three steps using gruppo but now we just divide the codon count by the max count in the synonymous group
        # gruppo=0
        # for j in SynonymousCodons[i]:
        #	gruppo = gruppo+CodonsDict_2[j]
        Lista_gruppo = []
        for j in SynonymousCodons[i]:
            Lista_gruppo.append(CodonsDict_2[j])  # /float(gruppo))
        # print Lista_gruppo
        max_lista = max(Lista_gruppo)  # this is the highest count in the codon count table
        for j in SynonymousCodons[i]:  # here j is the codon and i the synonymous group
            if max_lista == 0:
                CAI_dict[j] = [CodonsDict_2[j], 0]	#the two elements are the codon count and the count divided by the max
            else:
                CAI_dict[j] = [CodonsDict_2[j], CodonsDict_2[j]/float(max_lista)]	#the two elements are the codon count and the count divided by the max
                #Cai_list.append(Lista_gruppo[j])

    somma_num = 0
    somma_den = 0
    for i in list(CAI_dict.keys()):
        if CAI_dict[i][0] == 0:  # skip codons with 0 frequencies
            continue
        else:
            somma_num = somma_num + (CAI_dict[i][0] * log(CAI_dict[i][1]))  # here we sum all te (Fij * log(Wij))
            somma_den = somma_den + CAI_dict[i][0]  # here all the Fij are summed up
    Cai_index = exp(somma_num / somma_den)  # the new CAI index is now the division
    # prodotto**(float(1)/len(Cai_list))
    # devo salvare il cai, vedi se si può inserire in i insieme alle altre variabili
    # print prodotto**(float(1)/len(Cai_list)),prodotto
    # print Cai_list,len(Cai_list)
    # print CAI_dict,"\n\n\n\n\n", Cai_index

    L = sum(L)  # the sum of all nonconding will give the number of bases of noncoding regions
    g[(gb_file[gb_record].id).split(".")[0]] = (float(L) / LL) * 100  # and then dividing by the length we obtain the % of URs
    s[(gb_file[gb_record].id).split(".")[0]] = [str(AT), str(ATskew), str(GCskew), str(NG), str(Glength), str(Cai_index)]
# then we do the same for the custom annotation datasets. Results from this block will be added to the same dictionary
# used for the genbank file so that at the end we have several dictionaries (one for each variable) with
# key = genome (from genbank or from custom annotation)
# value = variable (one of the variables evaluated like CAI or AtSKEW)

for key in ANN:  # ANN is a dictionary with genomes as key and [genome_length , [genes_annotations]] as value
    #DEVO FINIRE LA PARTE DELLA CUSTOM ANNOT PERCHe' DOVREMMO FARE INDICARE IL TIPO DI CROMOSOMA NEL CUSTOM FILE E MODIFICARE TUTTO QUI:
    #PER ORA LA PRIMA RIGA E' GENOME E GENOMELIENGTH; DOVREBBE DIVENTARE A 3 CON ANCHE IL CIRCULARITY
    #if gb_file[gb_record].annotations["topology"] == "linear":   #linear chromosomes will be 1 and circular 0
    #    e[key] = 1
    #else:
    #    e[key] = 0
    e[key] = 0   #for now we assume topology is circular for the custom annotated genomes
    CDS = 0  # this is used to count how mani CDS were found in the custom annotation
    CodonsDict_2 = dict(CodonsDict)  # we reset the codon count for each custom annotated genome
    fasta = SeqIO.to_dict(
        SeqIO.parse(open(args.fasta), "fasta"))  # here we open the fasta file to obtain the genomes sequences
    F = ANN[key][1]  # F now is the number of annotated genes (instead of the features in genbank)
    LL = int(ANN[key][0])  # genome's length
    L = [1 for i in range(LL)]
    f[key] = []  # here we update f and h dictionaries with the values from the custom annotated genomes
    h[key] = [[], []]
    NG = 0
    # for AT content
    AT = 100 - GC(fasta[key].seq)
    # for AT skew
    As = int(fasta[key].seq.upper().count("A"))  # we count the As and Ts from the alignment file
    Ts = int(fasta[key].seq.upper().count("T"))  # we use upper because some fasta can be lower case
    ATskew = float(As - Ts) / (As + Ts)
    # for GC
    Gs = int(fasta[key].seq.upper().count("G"))  # we count the Gs and Cs from the alignment file
    Cs = int(fasta[key].seq.upper().count("C"))
    # for GC skew
    GCskew = float(Gs - Cs) / (Gs + Cs)
    # for genome length
    Glength = int(ANN[key][0])  # equal to LL, as before
    for i in range(0, len(F)):  # here is from 0 and not from 1 like in the genbank part
        # here we use all the annotated features (instead of selecting CDS,rRNA and tRNA like for genbanks)
        # for number of genes
        NG = NG + 1
        # for SU skew
        h[key][1].append(F[i][1])
        # for URs
        for y in range(int(F[i][2]) - 1, int(F[i][3])):  # we count 0 the regions with coding genes
            if L[y] == 1:
                L[y] = 0
        # for CAI
        if F[i][0].lower() in list(D.keys()):  # first we check if the name is in the dictionary (CDS and RRNA)
            if D[F[i][0].lower()] != "16S" and D[
                F[i][0].lower()] != "12S":  # we remove also the RRNA from the CDS count
                CDS = CDS + 1
            # can't use the append_gene function so we check for the strand
            if F[i][1] == "H":
                f[key].append(["+", D[F[i][0].lower()]])
            else:
                f[key].append(["-", D[F[i][0].lower()]])
            Seq = fasta[key].seq[int(F[i][2]) - 1:int(F[i][3])].upper()  # then we extract the feature sequence
            for j in range(3, len(Seq), 3):
                if Seq[j:j + 3] == "T":  # these exceptions are for the truncated mitochondrial stop codons
                    pass
                # CodonsDict_2["TAA"] = CodonsDict_2["TAA"] +1
                elif Seq[j:j + 3] == "TA":
                    pass
                # CodonsDict_2["TAA"] = CodonsDict_2["TAA"] +1
                elif Seq[j:j + 3] not in CodonsDict_2:  # if a codon is not in our table we skip it
                    continue
                else:
                    CodonsDict_2[Seq[j:j + 3]] = CodonsDict_2[Seq[j:j + 3]] + 1
    # for % AT in URs
    aa = 0  # number of As in URs
    tt = 0  # number of Ts in URs
    for i in range(len(L)):
        if L[i] == 1:
            if fasta[key].seq.upper()[i] == "A":
                aa = aa + 1
            elif fasta[key].seq.upper()[i] == "T":
                tt = tt + 1
    # and Median length of URs
    list_string = map(str, L)  # first we remap the integers to strings
    list_string_2 = "".join(list_string).split("0")  # then i split the joined string to obtain empty strings instead of 0s
    list_string = list(filter(None, list_string_2))  # and finally i filter out the empty strings
    urs_len = []
    for i in list_string:  # to evaluate median length of URs i use split to remove the 0s (coding portions)
        urs_len.append(len(i))
    # print(gb_file[gb_record].id)
    # print(L)
    # print(urs_len)
    # print(median(urs_len))
    # print(sum(urs_len) / len(urs_len))
    # print(list_string_2)  #questo comando è bellissimo. Fa vedere proprio le regioni non codificanti come 1 e le codificanti come stringhe vuote!
    l[key] = ((aa + tt) / float(sum(L))) * 100
    q[key] = median(urs_len)
    # OLD
    # now for CAI we find the frequencies of each codon for a residue
    # for i in SynonymousCodons.keys():
    # 	gruppo=0
    # 	Lista_gruppo= []
    # 	for j in SynonymousCodons[i]:
    # 		gruppo = gruppo+CodonsDict_2[j]
    # 	for j in SynonymousCodons[i]:
    # 		Lista_gruppo.append(CodonsDict_2[j]/float(gruppo))
    # 	max_lista = max(Lista_gruppo)
    # 	for j in range(len(Lista_gruppo)):
    # 		Lista_gruppo[j]=Lista_gruppo[j]/max_lista
    # 		Cai_list.append(Lista_gruppo[j])
    # prodotto=1
    # for i in Cai_list:
    # 	if i == 0:
    # 		continue
    # 	else:
    # 		prodotto = prodotto * i
    # Cai_index = prodotto**(float(1)/len(Cai_list))
    # NEW
    CAI_dict = {}
    for i in list(
            SynonymousCodons.keys()):  # HERMES1 version (commented) had three steps using gruppo but now we just divide the codon count by the max count in the synonymouse group
        # gruppo=0
        # for j in SynonymousCodons[i]:
        #	gruppo = gruppo+CodonsDict_2[j]
        Lista_gruppo = []
        for j in SynonymousCodons[i]:
            Lista_gruppo.append(CodonsDict_2[j])  # /float(gruppo))
        print(Lista_gruppo)
        max_lista = max(Lista_gruppo)  # this is the highest count in the codon count table
        for j in SynonymousCodons[i]:  # here j is the codon and i the synonymous group
            if max_lista == 0:
                CAI_dict[j] = [CodonsDict_2[j], 0]	#the two elements are the codon count and the count divided by the max
            else:
                CAI_dict[j] = [CodonsDict_2[j], CodonsDict_2[j]/float(max_lista)]	#the two elements are the codon count and the count divided by the max
                #Cai_list.append(Lista_gruppo[j])

    somma_num = 0
    somma_den = 0
    for i in list(CAI_dict.keys()):
        if CAI_dict[i][0] == 0:  # skip codons with 0 frequencies
            continue
        else:
            somma_num = somma_num + (CAI_dict[i][0] * log(CAI_dict[i][1]))  # here we sum all te (Fij * log(Wij))
            somma_den = somma_den + CAI_dict[i][0]  # here all the Fij are summed up
    Cai_index = somma_num / somma_den  # the new CAI index is now the division

    L = sum(L)
    g[key] = (float(L) / LL) * 100
    s[key] = [str(AT), str(ATskew), str(GCskew), str(NG), str(Glength), str(Cai_index)]

    print("Genome: " + str(key) + "\nCDS: " + str(CDS) + " || Genes: " + str(NG))

first = []

# to obtain the correct gene order
for key in f:
    for j in range(len(f[key])):
        if f[key][j][1] == 'CO1' and f[key][j][0] == '+':
            if [y[1] for y in f[key]].index('CO1') == 0:  # if cox1 is the first gene and it's on the +
                c = 0
            else:  # if not, we move all the genes before cox1 to the end
                for k in range(0, [y[1] for y in f[key]].index('CO1')):
                    f[key].insert(len(f[key]), f[key].pop(0))
        elif f[key][j][1] == 'CO1' and f[key][j][0] == '-':
            if [y[1] for y in f[key]].index('CO1') == len(f[key]):  # if cox1 is the last gene on the - strand
                f[key].reverse()  # we just need to reverse the order
                for k in f[key]:
                    change_sign(k)  # and change the signs
            else:
                for k in reversed(list(range([y[1] for y in f[key]].index('CO1') + 1, len(f[key])))):
                    f[key].insert(0, f[key].pop(len(f[key]) - 1))
                f[key].reverse()
                for k in f[key]:
                    change_sign(k)
    first.append("_".join([x[1] for x in f[key]]))

# to evaluate AMIGA score
second = set(first)
third = {}
fourth = {}
for i in second:
    third[i] = 0
for key in f:
    j = "_".join([x[1] for x in f[key]])
    third[j] += 1
for key in f:
    j = "_".join([x[1] for x in f[key]])
    fourth[key] = (third[j] - 1) / float(len(f) - 1)

# to evaluate the absolute value of SU skew
for key in h:
    m, p = 0, 0
    for i in h[key][1]:
        if i == 1:
            p = p + 1
        else:
            m = m + 1
    h[key][0] = abs(float(p - m) / (p + m))

# to find Root-to-tip distance

# First we rename the root with a non ambigous name (to avoid the 100 bootstrap score for the root)
args.tree = open(args.tree).readline().split(")")
args.tree[len(args.tree) - 1] = "root;\n"
args.tree = ")".join(args.tree)

# then we use the tree to evaluate the root-to-tip distance
t = Tree(args.tree, format=1)
T = {}

for el in t.get_leaves():
        T[N[el.name]] = round(el.get_distance(t.name), 3)
# ML_distance
W = {}
W[N[args.outgroup]] = "NA"  # the distance between the outgroup and itself is set here.
optionals = [args.threads, ]
if args.partition:
    pa_bool = True
else:
    pa_bool = False

subprocess.call(["./raxml", "-f", "x", "-s", args.alignment, "-n", output, "-m", args.model, "-q" * pa_bool,
                 str(args.partition) * pa_bool, "-T", args.threads, "-p",
                 "123456"])  # ,stdout=open(os.devnull, "wb"))	#questo l'ho commentato per vedere l'errore
distances = open("RAxML_distances." + output).readlines()

for i in range(len(distances)):
    distances[i] = distances[i].split()
    b = [distances[i][0], distances[i][1]]
    if args.outgroup in b:
        b.remove(args.outgroup)
        W[N[b[0]]] = distances[i][2]


# save all variables in one file, it is repeated: the first block if colors are indicated, second block if no color is indicated
#block1
if len(O) > 0:  # if len(O) > 0 means there is a third column in the species file, and we have to add the Shade column to HERMES_variables.txt file.
    file_for_R.write("ID\tURs\tAMIGA\tSUskew\tRtoTdist\tMLdist\tAT\tATskew\tGCskew\tGenes\tLength\tCAI\tCircularity\tURs_AT\tURs_MedLen\tShades\n")
    for key in sorted(O, key=O.get, reverse=True):
        if key == N[args.outgroup]:
            continue
        else:
            if key in list(T.keys()):
                file_for_R.write(
                    M[key] + "\t" + str(g[key]) + "\t" + str(fourth[key]) + "\t" + str(h[key][0]) + "\t" + str(
                        T[key]) + "\t" + str(W[key]) + "\t" + "\t".join(s[key]) + "\t" + str(e[key])+ "\t" + str(
                        l[key]) + "\t" + str(q[key]) +"\t" + str(O[key]) + "\n")
            else:  # if the species isn't in the tree
                file_for_R.write(M[key] + "\t" + str(g[key]) + "\t" + str(fourth[key]) + "\t" + str(
                    h[key][0]) + "\t" + "NA" + "\t" + str(W[key]) + "\t" + "\t".join(s[key]) + "\t" + str(e[key]) + "\t" + str(
                    l[key]) + "\t" + str(q[key]) +"\t" + str(O[key]) + "\n")
#block2
else:  # if O is empty, means that there isn't a third column in the species file and that there is no column O[key] to be added (shades).
    file_for_R.write("ID\tURs\tAMIGA\tSUskew\tRtoTdist\tMLdist\tAT\tATskew\tGCskew\tGenes\tLength\tCAI\tCircularity\tURs_AT\tURs_MedLen\n")
    for key in f:   #these are the species and for each species we obtain all the values in each dictionary (variables)
        if key == N[args.outgroup]:
            continue
        else:
            if key in list(T.keys()):   #here we use the species (keys) to obtain the values in each dictionary (variable)
                file_for_R.write(
                    M[key] + "\t" + str(g[key]) + "\t" + str(fourth[key]) + "\t" + str(h[key][0]) + "\t" + str(
                        T[key]) + "\t" + str(W[key]) + "\t" + "\t".join(s[key]) + "\t" + str(e[key]) + "\t" + str(
                        l[key]) + "\t" + str(q[key]) + "\n")
            else:
                file_for_R.write(M[key] + "\t" + str(g[key]) + "\t" + str(fourth[key]) + "\t" + str(
                    h[key][0]) + "\t" + "NA" + "\t" + str(W[key]) + "\t" + "\t".join(s[key]) + "\t" + str(e[key]) + "\t" + str(
                    l[key]) + "\t" + str(q[key]) + "\n")
file_for_R.close()

if type(args.alpha) == float:  # if alpha isn't specified in the command line, it is equal to 'None' and so it isn't a float.
    subprocess.call(["Rscript", "--vanilla", "HERMES-v" + version + ".R", str(args.alpha)])
else:
    subprocess.call(["Rscript", "--vanilla", "HERMES-v" + version + ".R"])

# os.rename("./Rplots.pdf","./Results/HERMES.pdf") old version
os.remove("RAxML_distances." + output)
os.remove("RAxML_info." + output)
os.remove("RAxML_parsimonyTree." + output)
