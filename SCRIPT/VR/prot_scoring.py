## read to prot fasta file
import os
import subprocess

class HSP:
    def __init__(self, query_start, query_end, subject_start, subject_end, bitscore, positives, identity, length, pos_query=None, pos_subject=None):
        self.query_start = query_start
        self.query_end = query_end
        self.subject_start = subject_start
        self.subject_end = subject_end
        self.bitscore = bitscore
        self.positives = positives
        self.identity = identity
        self.length = length
        self.pos_query = pos_query
        self.pos_subject = pos_subject
        

    def __repr__(self):
        return (f"HSP(query_start={self.query_start}, query_end={self.query_end}, "
                f"subject_start={self.subject_start}, subject_end={self.subject_end}, "
                f"bitscore={self.bitscore}, positives={self.positives}, identity={self.identity}, "
                f"length={self.length}, pos_query={self.pos_query}, pos_subject={self.pos_subject})")


def run_blast(fasta_fileQ, fasta_fileS, output_file):
    subprocess.run([
        "blastp", "-query", fasta_fileQ, "-subject", fasta_fileS,
        "-outfmt", "6 qstart qend sstart send bitscore positive nident length",
        "-out", output_file
    ], check=True)

def parse_blast_output(filename):
    hsps = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            hsps.append(HSP(
                int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]),
                float(parts[4]), int(parts[5]), int(parts[6]), int (parts[7])
            ))
    return hsps

def add_position_info(hsps):
    hsps.sort(key=lambda h: h.query_start)
    for i, h in enumerate(hsps):
        h.pos_query = i

    hsps.sort(key=lambda h: h.subject_start)
    for i, h in enumerate(hsps):
        h.pos_subject = i

def is_compatible_hsp_list(comp_hsps):
    for i in range(1, len(comp_hsps)):
        prev_hsp = comp_hsps[i - 1]
        curr_hsp = comp_hsps[i]
        if not (curr_hsp.query_start >= prev_hsp.query_end and
                curr_hsp.subject_start >= prev_hsp.subject_end and
                curr_hsp.pos_subject > prev_hsp.pos_subject):
            return False
    return True

def construct_compatible_hsp_list(hsps):
    hsps.sort(key=lambda h: h.bitscore, reverse=True)
    comp_hsps = []

    for hsp in hsps:
        temp_hsps = comp_hsps + [hsp]
        temp_hsps.sort(key=lambda h: h.subject_start)

        if is_compatible_hsp_list(temp_hsps):
            comp_hsps = temp_hsps

    return comp_hsps

def measure_prot_similarity(fasta_file1, fasta_file2):
    run_blast(fasta_file1, fasta_file2, "blast_output.tsv")

    hsps = parse_blast_output("blast_output.tsv")
    add_position_info(hsps)
    comp_hsps = construct_compatible_hsp_list(hsps)

    return comp_hsps

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python script.py <fasta_file1> <fasta_file2>")
        sys.exit(1)

    fastaQuery = sys.argv[1]
    fastaSubject = sys.argv[2]

    comp_hsps = measure_prot_similarity(fastaQuery, fastaSubject)
    sumBitscore=0
    sumPositives=0
    sumIdentity=0
    sumLength=0
    for hsp in comp_hsps:
        #print(hsp)
        sumBitscore += hsp.bitscore;
        sumPositives += hsp.positives;
        sumIdentity += hsp.identity;
        sumLength += hsp.length;
    print(f"{sumPositives} {sumIdentity} {sumLength} {sumBitscore}");