from Bio import SeqIO
import os


home_folder=os.path.expanduser('~')
gb_file = f"{home_folder}/ecoli-genomes/sequence_ecoli.gb"
for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
      print(gb_record.seq[270575:270608])