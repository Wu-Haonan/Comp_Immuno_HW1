import csv
from collections import Counter
from collections import defaultdict


def count_vj_pairs(file_path):
    vj_counter = Counter()

    with open(file_path, newline='') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')

        for row in reader:
            v_hit = row['V_hit'].split('*')[0]
            j_hit = row['J_hit'].split('*')[0]
            
            vj_pair = (v_hit, j_hit) 
            vj_counter[vj_pair] += 1

    return vj_counter

def print_vj_frequencies(vj_counter):
    print("VJ Pair Frequency:")
    for vj_pair, count in vj_counter.items():
        print(f"V_hit: {vj_pair[0]}, J_hit: {vj_pair[1]}, Frequency: {count}")

def get_top_v_genes(file_path, top_n=10):
    gene_data = []

    with open(file_path, newline='') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')

        for row in reader:
            gene_name = row['Gene_name']
            frequency = int(row['Multiplicity'])  
            gene_data.append((gene_name, frequency))
    sorted_genes = sorted(gene_data, key=lambda x: x[1], reverse=True)[:top_n]

    return sorted_genes

def parse_gene_name(gene_name):
    return gene_name.split('*')[0]

def count_mutations_per_alignment(shm_file_path, genes_of_interest):
    v_genes = {gene[0] for gene in genes_of_interest}
    
    mutation_counts = defaultdict(list)
    current_gene = None
    current_count = 0
    
    with open(shm_file_path, 'r') as file:
        next(file)
        
        for line in file:
            if line.startswith('Read_name'):
                # prase last gene
                if current_gene in v_genes:
                    mutation_counts[current_gene].append(current_count)
            
                gene_part = [part for part in line.split() if part.startswith('Gene_name:')][0]
                current_gene = gene_part.split(':')[1].split('*')[0]
                current_count = 0
            
            else:
                if current_gene in v_genes:
                    current_count += 1
        
        # last alignment
        if current_gene in v_genes:
            mutation_counts[current_gene].append(current_count)
    
    return mutation_counts

def print_mutation_stats(mutation_counts):
    print("\nMutation Statistics per Gene:")
    print("============================")
    
    for gene in mutation_counts:
        counts = mutation_counts[gene]
        print(f"\nGene: {gene}")
        print(f"Number of alignments with mutations: {len(counts)}")
        print(f"Mutation counts per alignment: {counts}")

if __name__ == "__main__":
    default = './results/plasma'
    cdr_details_path = default + '/cdr_details.txt' 
    vj_counter = count_vj_pairs(cdr_details_path)
    print_vj_frequencies(vj_counter)

    v_usage_path = default + '/visualizer/plots/gene_usage_plots/v_usage.txt' 
    top_n = 10 
    sorted_genes = get_top_v_genes(v_usage_path, top_n)
    print(sorted_genes)

    shm_details_path = default + '/shm_details.txt' 
    mutation_counts = count_mutations_per_alignment(shm_details_path, sorted_genes)
    print_mutation_stats(mutation_counts)