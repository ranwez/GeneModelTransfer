import sys
import subprocess

def run_blastp(predicted_fasta, template_fasta, output_file):
    """Run BLASTP and save the output to a file."""
    blastp_cmd = [
        'blastp',
        '-query', predicted_fasta,
        '-subject', template_fasta,
        '-out', output_file,
        '-outfmt', '0'
    ]
    subprocess.run(blastp_cmd, check=True)

def parse_blast_output(blast_output):
    """Parse the BLASTP output and construct the similarity table."""
    in_alignment_section = False
    alignment_started = False
    pred_prod_length = 0
    current_hsp_start = 0
    current_hsp_end = 0
    table = []

    with open(blast_output, 'r') as f:
        for line in f:
            if line.startswith("Query= "):
                next(f)  # Skip the blank line
                pred_prod_length = int(next(f).strip().split('=')[1])
                table = [0] * pred_prod_length
                continue
            
            if line.startswith(">"):
                in_alignment_section = True
                continue
            
            if in_alignment_section and line.startswith(" Score ="):
                alignment_started = True
                continue
            
            if alignment_started:
                if line.startswith("Query  "):
                    parts = line.split()
                    current_hsp_start = int(parts[1]) - 1
                    current_hsp_end = int(parts[3])
                    query_seq = parts[2]
                    
                    # Calculate indices for extraction
                    seq_start = line.find(query_seq)
                    seq_end = seq_start + len(query_seq)

                    # Extract the subsequences using the indices
                    alignment_line = next(f)
                    subject_line = next(f)
                    alignment_seq = alignment_line[seq_start:seq_end]
                    subject_seq = subject_line[seq_start:seq_end]
                
                    for i in range(current_hsp_start, current_hsp_end):
                        query_index = i - current_hsp_start
                        if query_seq[query_index] == subject_seq[query_index]:
                            table[i] = 2
                        elif alignment_seq[query_index] == '+' and table [i] == 0:
                            table[i] = 1
    num_identical = table.count(2)
    num_positive = sum(1 for x in table if x >= 1)
    
    return num_identical, num_positive

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <predicted_prot.fasta> <template_prot.fasta>")
        sys.exit(1)
    
    predicted_fasta = sys.argv[1]
    template_fasta = sys.argv[2]
    output_file = "blast_output.txt"
    
    run_blastp(predicted_fasta, template_fasta, output_file)
    num_identical, num_positive = parse_blast_output(output_file)
    #print(f"Identical Matches: {num_identical}")
    #print(f"Positive Matches: {num_positive}")
    print(f"{num_positive} {num_identical} 0 0");

if __name__ == "__main__":
    main()