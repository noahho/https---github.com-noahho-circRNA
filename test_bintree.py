from bintrees import AVLTree

from Bio.Seq import Seq
output_handle = open("mitofrags.fasta", "w")
SeqIO.write(mito_frags, output_handle, "fasta")
output_handle.close()