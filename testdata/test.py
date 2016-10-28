from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
sequenceset = ['AGATAAGTTCACGTTACCAT', 'GTATT', 'TAGATGAAGCGGGATAGTCTTTTTCTGATATGCACTTATCAGTTCACTAGCAGT', 'ACTGAACGTGATTGATGAAGCT', 'ATCTA']
record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"))
records = [record]
SeqIO.write(records, "aa.txt", "fasta")