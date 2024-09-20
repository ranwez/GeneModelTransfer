import csv
import sys
import argparse
"""
@author: Vincent Ranwez
@description: Fix LRR annotation files (assuming each gene as a single mrna, no alternative splicing)
"""

def modify_feature_ids(gff_file, ids_prefix, remove_comments):
    # Read the GFF file and store features
    gene_features = []
    cds_features = []
    mrna_features = []
    with open(gff_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if not row[0].startswith('#'):
                if row[2]=="gene":
                    gene_features.append(row)
                elif row[2]=="mRNA":
                    mrna_features.append(row)
                elif row[2]=="CDS": 
                    cds_features.append(row)

    # Create dictionaries to hold feature information and relationships
    gene_id_mapping = {}
    mrna_id_mapping = {}
    cds_count = {}
    mrna_count = {}
    # Create dictionaries to hold hierarchical relationships
    gene_to_mrna = {}
    mrna_to_cds = {}
    orphan_mrna=[]
    orphan_cds=[]


    # Process genes, mRNAs, and CDS features. Since we need parent information we process all genes first, then all mrna etc.
    print ("nb gene: "+ str(len(gene_features)))
    print ("nb mrna: " + str(len(mrna_features)))
    print ("nb cds: " + str(len (cds_features)))
    
    # Sort features by start position to get correct cds/exon indices
    gene_features.sort(key=lambda x: int(x[3]))
    mrna_features.sort(key=lambda x: int(x[3]))
    cds_features.sort(key=lambda x: int(x[3]))

    for row in gene_features:
        process_gene(row, gene_id_mapping)
    for row_id, row in enumerate(mrna_features):
        process_mrna(row, row_id, gene_id_mapping, mrna_id_mapping, gene_to_mrna, orphan_mrna,mrna_count)
    for row_id, row in enumerate(cds_features):
        process_cds(row, row_id, mrna_id_mapping, cds_count, mrna_to_cds, orphan_cds)

    # Organize and sort mRNAs and CDS for each gene
    sorted_features = []
    for gene in gene_features:
        gene_id = gene[8].split(';')[0].split('=')[1]
        sorted_features.append(gene)
        if gene_id in gene_to_mrna:
            for mrna_row_id in gene_to_mrna[gene_id]:
                mrna = mrna_features[mrna_row_id]
                if(mrna[3]<gene[3] or mrna[4]>gene[4]):
                    print (f"ERROR: incompatible mRNA gene bounds for {mrna} \n")
                append_mrna(mrna, sorted_features,mrna_to_cds, cds_features)
    
    # Add remaining mRNAs (those without parents)
    for mrna_row_id in orphan_mrna:
        mrna = mrna_features[mrna_row_id]
        append_mrna(mrna, sorted_features,mrna_to_cds, cds_features)

    # Add remaining CDS (those without parents)
    for cds_row_id in orphan_cds:
        cds_row=cds_features[cds_row_id]
        sorted_features.append(cds_row)
        sorted_features.append(to_exon(cds_row))
        
    # optionnally add prefix to IDs and remove non essential comment
    if ids_prefix or remove_comments:
        for row in sorted_features:
            (id, parent_id, others_attr)=extract_id_and_pid(row[8])
            real_ids_prefix = f"{ids_prefix}_" if ids_prefix else ''

            new_id=f"{real_ids_prefix}{id}"
            new_parent_id = f"{real_ids_prefix}{parent_id}" if parent_id else ''
            new_others_attr = '' if remove_comments else others_attr
            row[8] = update_id_and_pid(new_id, new_parent_id, new_others_attr)
    return (sorted_features)
   
# using csv writer, lead to new line issues even when replacing \r in last field ..
def write_gff(gff_rows, output_file_name):
    with open(output_file_name, 'w', newline='\n') as file:
        for row in gff_rows:
            line = '\t'.join(row) + '\n'  # Construct each line manually
            line.replace('\r\n', '\n').replace('\r', '')
            file.write(line)

def to_exon(cds):
    # create a copy of the CDS and turn it to an exon feature
    exon = list(cds)
    exon[2] = 'exon'
    (cds_id, parent_mrna_id, others_attr)=extract_id_and_pid(cds[8])
    exon_id = cds_id.replace('CDS', 'exon')
    exon[8] = update_id_and_pid(exon_id, parent_mrna_id, others_attr)
    return (exon)
    
def append_mrna(mrna, sorted_features,mrna_to_cds, in_cds_features):
    sorted_features.append(mrna)
    mrna_id = mrna[8].split(';')[0].split('=')[1]
    last_bounds=0
    if mrna_id in mrna_to_cds:
        for cds_row_id in mrna_to_cds[mrna_id]:
            cds = in_cds_features[cds_row_id]
            if(int(cds[3])<int(mrna[3]) or int(cds[4])>int(mrna[4])):
                print (f"ERROR: incompatible cds / mRNA bounds for \n {cds} \n")
            if (int(cds[4])< last_bounds):
                print (f"ERROR: overlapping cds \n {cds} \n")
            last_bounds = int(cds[4])
            sorted_features.append(cds)
            sorted_features.append(to_exon(cds))

def process_gene(row, gene_id_mapping):
    seq_id, start, attributes = row[0], row[3], row[8]
    (former_gene_id, parent_gene_id, others_attr)=extract_id_and_pid(row[8])
    start10=start.zfill(10)
    gene_id = f"{seq_id}_{start10}"
    gene_id_mapping[former_gene_id] = gene_id
    row[8] = update_id_and_pid(gene_id, parent_gene_id, others_attr)

def extract_id_and_pid(infos):
    id_attr = [attr for attr in infos.split(';') if attr.startswith('ID=')]
    id_value = id_attr[0].split('=')[1] if len(id_attr) == 1 else ''

    pid_attr = [attr for attr in infos.split(';') if attr.startswith('Parent=')]    
    pid_value = pid_attr[0].split('=')[1] if len(pid_attr) == 1 else ''

    others_attr=[attr for attr in infos.split(';') if not (attr.startswith('ID=') or attr.startswith('Parent='))]
    
    return (id_value, pid_value, others_attr)

def update_id_and_pid(new_id, new_parent_id, others_attr):
    attributes_parts = [f"ID={new_id}"]
    if new_parent_id !='':
        attributes_parts.append(f"Parent={new_parent_id}")
    others_attr_str = ';'.join(others_attr)
    if others_attr_str:  # Add the joined string only if it's not empty
        attributes_parts.append(others_attr_str)
    return ';'.join(attributes_parts)

def process_mrna(row, row_id, gene_id_mapping, mrna_id_mapping, gene_to_mrna, orphan_mrna, mrna_count):
    start, end = row[3], row[4]
    # Extract the current parent and mRNA IDs
    (former_mrna_id, parent_gene_id, others_attr)=extract_id_and_pid(row[8])
    new_parent_id = ''
    if (parent_gene_id !='') and (parent_gene_id in gene_id_mapping):
        new_parent_id = gene_id_mapping.get(parent_gene_id)
        gene_to_mrna.setdefault(new_parent_id, []).append(row_id)
    else:
        orphan_mrna.append(row_id)
        print("\nWarning: No valid parent gene ID found for mRNA. Removing parent attribute.")
        print(row)

   # Count the position for mRNA within its parent gene
    if new_parent_id != '':
        mrna_count[new_parent_id] = mrna_count.get(new_parent_id, 0) + 1
        position = mrna_count[new_parent_id]
    else:
        position = 1  # Default position if no parent is found

    # Generate the new mRNA ID
    mrna_id = f"{new_parent_id}_mrna_{position}" if new_parent_id else f"mrna_{start}_{end}"
    mrna_id_mapping [former_mrna_id] = mrna_id
    row[8]=update_id_and_pid(mrna_id, new_parent_id, others_attr)
    

def process_cds(row, row_id, mrna_id_mapping, cds_count, mrna_to_cds, orphan_cds):
    (former_cds_id, parent_mrna_id, others_attr)=extract_id_and_pid(row[8])
    
    # Check if a valid parent mRNA attribute is available
    new_parent_id = ''
    if (parent_mrna_id != '') and (parent_mrna_id in mrna_id_mapping):
        new_parent_id = mrna_id_mapping.get(parent_mrna_id)
        mrna_to_cds.setdefault(new_parent_id, []).append(row_id)
    else:
        orphan_cds.append(row_id)
        print("\nWarning: No valid parent mRNA ID found for CDS. Removing parent attribute.")
        print(row)

    # Count the position for CDS within its parent mRNA
    if new_parent_id != '':
        cds_count[new_parent_id] = cds_count.get(new_parent_id, 0) + 1
        position = cds_count[new_parent_id]
    else:
        position = 1  # Default position if no parent is found

    # Generate the new CDS ID
    cds_id = f"{new_parent_id}_CDS_{position}" if new_parent_id else f"CDS_{position}"
    row[8]=update_id_and_pid(cds_id, new_parent_id, others_attr)

#-g /Users/ranwez/Desktop/BUG_CLEANER/bug_cleaner.gff -o /Users/ranwez/Desktop/BUG_CLEANER/bug_cleaner_out.gff
def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Modify feature IDs in a GFF file. \
        This script updates gene, mRNA, and CDS feature IDs in a GFF file based on their sequence ID and start positions. \
        It then sorts and organizes the features such that genes are sorted by start position, \
        followed by their corresponding mRNA and CDS entries. \
        An optionnal prefix (e.g. genome name) could be provided and added at the begining of each ID. \
        The script expect UTF8 encoded file you can convert using iconv bash command e.g. \
        iconv -f ISO-8859-1 -t UTF-8 input_latin1.gff -o input_utf8.gff')
    
    parser.add_argument("-g", "--gff", metavar='input_ut8.gff', type=str, required=True, help="path to the input GFF file.")
    parser.add_argument("-o", "--output", metavar='output_cleaned.gff', type=str, required=True, help="name of the reformated GFF file.")
    parser.add_argument("-p", "--prefix", metavar='DWSvevo1', type=str,required=False, help="a prefix to add at the begining of each feature ID.")
    parser.add_argument("-r", "--removeComments", action='store_true', help="If used, comments other than ID and Parent are removed.")

    # Check if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    updated_gff_rows=modify_feature_ids(args.gff, args.prefix or '', args.removeComments)
    

    write_gff(updated_gff_rows,args.output )    

if __name__ == "__main__":
    main()
