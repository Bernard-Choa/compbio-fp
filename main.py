from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import draw




# !! CHANGE PATH IF THIS CODE IS RUN ON OTHER COMPUTERS !!
MY_SEQDATA_DIR_PATH = "C:/Users/berna/VSCodeProjects/compbio-handson/msaTest/seqData/"

TXSeq = ''
with open(MY_SEQDATA_DIR_PATH+"TX.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        TXSeq += record.seq
TXSeqObj = SeqRecord(Seq(TXSeq), id="xantha")

TSSeq = ''
with open(MY_SEQDATA_DIR_PATH+"TS.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        TSSeq += record.seq
TSSeqObj = SeqRecord(Seq(TSSeq), id="swinhoei")

TMSeq = ''
with open(MY_SEQDATA_DIR_PATH+"TM.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        TMSeq += record.seq
TMSeqObj = SeqRecord(Seq(TMSeq), id="mirabilis")

TCuSeq = ''
with open(MY_SEQDATA_DIR_PATH+"TCu.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        TCuSeq += record.seq
TCuSeqObj = SeqRecord(Seq(TCuSeq), id="cupola")

TCySeq = ''
with open(MY_SEQDATA_DIR_PATH+"TCy.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        TCySeq += record.seq
TCySeqObj = SeqRecord(Seq(TCySeq), id="cylindrica")

TDSeq = ''
with open(MY_SEQDATA_DIR_PATH+"TD.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        TDSeq += record.seq
TDSeqObj = SeqRecord(Seq(TDSeq), id="deliqua")

TMaSeq = ''
with open(MY_SEQDATA_DIR_PATH+"TMa.fasta") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        TMaSeq += record.seq
TMaSeqObj = SeqRecord(Seq(TMaSeq), id="maricae")



MSA = MultipleSeqAlignment([TXSeqObj, TSSeqObj, TMSeqObj, TCuSeqObj, TMaSeqObj, TCySeqObj, TDSeqObj])

calc = DistanceCalculator('identity')
builder = DistanceTreeConstructor(calc)
pohon = builder.build_tree(MSA)
pohon.rooted = True

fig = draw(pohon)

# for record in SeqIO.parse("TX_1.fasta", "fasta"):
#     print(record.id)