import csv
from collections import Counter
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np


def count_vj_pairs(file_path):
    vj_counter = Counter()
    count = 0
    non_productive = 0
    with open(file_path, newline='') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')

        for row in reader:
            v_hit = row['V_hit'].split('*')[0]
            j_hit = row['J_hit'].split('*')[0]
            if row['Productive'] == '0':
                non_productive += 1
            vj_pair = (v_hit, j_hit) 
            vj_counter[vj_pair] += 1
            count += 1
    print(count)
    return vj_counter, (non_productive/count)

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
    SHM = defaultdict(list)
    current_gene = None
    current_count = 0
    
    with open(shm_file_path, 'r') as file:
        next(file)
        
        for line in file:
            if line.startswith('Read_name'):
                # prase last gene
                if current_gene in v_genes:
                    mutation_counts[current_gene].append(current_count)
                    SHM[current_gene].append(current_count/gene_length)
            
                gene_part = [part for part in line.split() if part.startswith('Gene_name:')][0]
                current_gene = gene_part.split(':')[1].split('*')[0]
                length_part = [part for part in line.split() if part.startswith('Gene_length:')][0]
                gene_length = int(length_part.split(':')[1])
                current_count = 0
            
            else:
                if current_gene in v_genes:
                    current_count += 1
        
        # last alignment
        if current_gene in v_genes:
            mutation_counts[current_gene].append(current_count)
            SHM[current_gene].append(current_count/gene_length)
    
    return mutation_counts, SHM

def print_mutation_stats(mutation_counts):
    print("\nMutation Statistics per Gene:")
    print("============================")
    
    for gene in mutation_counts:
        counts = mutation_counts[gene]
        print(f"\nGene: {gene}")
        print(f"Number of alignments with mutations: {len(counts)}")
        print(f"Mutation counts per alignment: {counts}")

def plot_shm_rates(shm_rates, output_path):
    genes = list(shm_rates.keys())
    data = [shm_rates[gene] for gene in genes]
    

    plt.figure(figsize=(12, 6))
    
    bp = plt.boxplot(data, tick_labels=genes)
    
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('SHM Rate')
    plt.xlabel('V Genes')
    plt.title('SHM Rates Distribution for Top 10 V Genes')
    
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def get_intersection_genes(sorted_genes_dict):
    gene_sets = []
    for sample_type, genes in sorted_genes_dict.items():
        gene_names = {gene[0] for gene in genes}
        gene_sets.append(gene_names)
    
    intersection = set.intersection(*gene_sets)
    return sorted(list(intersection))

def plot_intersection_shm(gene, shm_data_dict, output_path):
    plt.figure(figsize=(8, 6))

    data = []
    labels = []
    for sample_type in ['naive', 'memory', 'plasma']:
        if gene in shm_data_dict[sample_type]:
            data.append(shm_data_dict[sample_type][gene])
            labels.append(sample_type)

    bp = plt.boxplot(data, tick_labels=labels)
    
    plt.title(f'SHM Rates Distribution for {gene}')
    plt.ylabel('SHM Rate')
    plt.xlabel('Sample Type')
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

if __name__ == "__main__":
    all_sorted_genes = {}
    all_shm_data = {}
    
    for sample_type in ['plasma', 'naive', 'memory']:
        print(sample_type)
        default = './results/' + sample_type
        
        cdr_details_path = default + '/cdr_details.txt' 
        vj_counter, non_productive_rate = count_vj_pairs(cdr_details_path)
        print_vj_frequencies(vj_counter)
        print("fraction of non-productive: ", non_productive_rate)

        v_usage_path = default + '/visualizer/plots/gene_usage_plots/v_usage.txt' 
        sorted_genes = get_top_v_genes(v_usage_path, 10)
        all_sorted_genes[sample_type] = sorted_genes
        print(sorted_genes)

        shm_details_path = default + '/shm_details.txt' 
        mutation_counts, SHM = count_mutations_per_alignment(shm_details_path, sorted_genes)
        all_shm_data[sample_type] = SHM
        print_mutation_stats(mutation_counts)
        print_mutation_stats(SHM)

        output_file = sample_type + '_shm_rates_boxplot.png'
        plot_shm_rates(SHM, output_file)
    
    intersection_genes = get_intersection_genes(all_sorted_genes)
    print("\nIntersection genes:", intersection_genes)

    for gene in intersection_genes:
        output_file = f'intersection_{gene}_comparison.png'
        plot_intersection_shm(gene, all_shm_data, output_file)