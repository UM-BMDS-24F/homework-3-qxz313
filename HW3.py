from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from io import StringIO

#On the Windows command line, enter makeblastdb -in mouse.fa -out mouse -dbtype prot to create the mouse database.

human_fa = "./human.fa"
mouse_db = "./mouse"
output_file = "blast_results.txt"

human_sequences = list(SeqIO.parse(human_fa, "fasta"))

for nt in human_sequences:
    print(nt)
    
with open(output_file, "w") as f:
    for nt in human_sequences:
        temp_query = "temp_query.fa"
        with open(temp_query, "w") as temp_f:
            SeqIO.write(nt, temp_f, "fasta")
        
        blastp_cline = NcbiblastpCommandline(cmd="C:/Program Files/NCBI/blast-2.16.0+/bin/blastp",
                                             query=temp_query, db=mouse_db, evalue=0.001, outfmt=5,out="out.xml")
        stdout, stderr = blastp_cline()
        result_handle = open("out.xml")
        blast_records =NCBIXML.parse(result_handle)
        best_alignment = None
        best_bitscore = 0
        i=0
        for blast_record in blast_records:
            i+=1
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.bits > best_bitscore:
                        best_alignment = alignment
                        best_bitscore = hsp.bits
 
        f.write(f"Human Sequence ID: {nt.id}\n")
        if best_alignment:
            best_hsp = best_alignment.hsps[0]

            f.write(f"Mouse Sequence ID: {best_alignment.hit_def}\n")
            f.write(f"Alignment:\n{best_hsp.sbjct}\n")
            f.write(f"E-value: {best_hsp.expect}\n")
            f.write(f"Bitscore: {best_hsp.bits}\n")
        else:
            f.write("No significant alignments found.\n")

        f.write("\n" + "="*40 + "\n\n")
