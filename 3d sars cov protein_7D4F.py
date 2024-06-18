import Bio
import pylab
import urllib
import pandas as pd
import nglview as nv
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Alphabet import IUPAC
from collections import Counter
from Bio.Data import CodonTable
from Bio import SeqIO, SearchIO
from Bio.PDB import PDBParser,MMCIFParser
from Bio.SeqUtils import GC,molecular_weight
from Bio.Alphabet import generic_dna,generic_rna,generic_protein

dir(Bio)
seq_file_read = SeqIO.read("Sequence_data/sequence.fasta","fasta")
seq_file_read
type(seq_file_read)

mRNA = seqfromfile.transcribe()
mRNA.back_transcribe()

print(CodonTable.unambiguous_dna_by_id[1])
common_amino = Counter(protein_seq)
common_amino.most_common(10)

del common_amino['*']
pylab.bar(common_amino.keys(),common_amino.values())
pylab.title('protein sequence frequency')
pylab.xlabel('amino acids')
pylab.ylabel('frequency')

protein_list = [str(i) for i in protein_seq.split('*')]
protein_list[:10]

large_proteins = [x for x in protein_list if len(x)>10]
df = pd.DataFrame({'protein_seq':large_proteins})

df['length'] = df['protein_seq'].apply(len)
df.head()

df.sort_values(by = ['length'],ascending = False)[:10]

one_large_protein = df.nlargest(1,'length')
single_prot = one_large_protein.iloc[0,0]
single_prot

with open ('Sequence_data/single_prot.fasta','w') as file:
    file.write('>large protein\n'+single_prot)

read = SeqIO.read("Sequence_data/single_prot.fasta","fasta")

%%time

result_handle = NCBIWWW.qblast("blastp","pdb",read.seq)
blast_qresult = SearchIO.read(result_handle,"blast-xml")

seqid = blast_qresult[0]

details = seqid[0]

print(f"\
Sequence ID:{seqid.id}\n\
description:{seqid.description}\n\
E value:    {details.evalue} \n\
Bit Score:  {details.bitscore}\n\
")

print(f"alignment:\n{details.aln}")

seqid.id
seqid.id.split('|')[1]
urllib.request.urlretrieve('https://files.rcsb.org/download/7D4F.pdb','Sequence_data/7D4F.pdb')
parser = PDBParser()
structure = parser.get_structure("7D4F",'Sequence_data/7D4F.pdb')
view = nv.show_biopython(structure)
view.render_image()

nv.show_biopython(structure,gui = True)