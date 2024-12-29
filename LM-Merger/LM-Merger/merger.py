import re
import boolean
from parse_utils import parse_expression, simplify_expression

def merge_networks(networks, method="OR", descriptive=True):
    """
    Merges multiple Boolean network models using the specified method (OR, AND, or Inhibitor Wins).
    Outputs a description of the merged model if descriptive is True.
    """
    if method not in ["OR", "AND", "Inhibitor Wins"]:
        raise ValueError("Invalid method. Use 'OR', 'AND', or 'Inhibitor Wins'.")

    merged_network = {}
    overlap_genes = set()
    descriptions = {}
    warning_gene = {}
    algebra = boolean.BooleanAlgebra()

    # Track the number of genes in each individual model
    individual_gene_counts = [len(network) for network in networks]

    # Iterate through each network
    for idx, network in enumerate(networks):
        for gene, expression in network.items():
            if gene in merged_network:
                overlap_genes.add(gene)
                # Ignore if the expression is the same
                if algebra.parse(merged_network[gene]) == algebra.parse(expression):
                    if descriptive:
                        descriptions.setdefault(gene, []).append((f"Model {idx + 1}", expression))
                        descriptions[gene].append(("Merged", merged_network[gene]))
                    continue
                # Determine how to merge the expressions based on the method
                if method == "OR":
                    merged_expression = f"({merged_network[gene]})|({expression})"
                elif method == "AND":
                    merged_expression = f"({merged_network[gene]})&({expression})"
                elif method == "Inhibitor Wins":
                    activators, inhibitors = parse_expression(expression)
                    existing_activators, existing_inhibitors = parse_expression(merged_network[gene])
                    combined_activators = list(set(existing_activators + activators))
                    combined_inhibitors = list(set(existing_inhibitors + inhibitors))
                    
                    activator_expr = '|'.join(filter(None, combined_activators)) if combined_activators else ""
                    inhibitor_expr = '|'.join(filter(None, combined_inhibitors)) if combined_inhibitors else ""
                    
                    activator_genes = set(re.split(r'\|', activator_expr))
                    inhibitor_genes = set(re.split(r'\|', inhibitor_expr))
                    if activator_genes & inhibitor_genes:
                        overlappings = list(activator_genes & inhibitor_genes)
                        for overlapping in overlappings:
                            activator_expr = re.sub(rf'(\||&)?{overlapping}(\||&)?', '', activator_expr).strip('|&')
                        warning_gene[gene] = overlappings
                    
                    if inhibitor_expr:
                        if '&' in inhibitor_expr or '|' in inhibitor_expr:
                            merged_expression = f"({activator_expr})&!({inhibitor_expr})"
                        else:
                            merged_expression = f"({activator_expr})&!{inhibitor_expr}"
                    else:
                        merged_expression = activator_expr
                else:
                    merged_expression = expression
                
                # Simplify the merged expression
                # print(gene)
                # print(merged_expression)
                merged_network[gene] = simplify_expression(merged_expression)

                if descriptive:
                    descriptions.setdefault(gene, []).append((f"Model {idx + 1}", expression))
                    descriptions[gene].append(("Merged", merged_network[gene]))
                
            else:
                merged_network[gene] = expression
                if descriptive:
                    descriptions[gene] = [(f"Model {idx + 1}", expression)]

    if descriptive:
        print(f"Merging Method: {method}")
        print(f"Total Genes in Merged Network: {len(merged_network)}")
        print(f"Number of Genes in Each Individual Model:")
        for i, count in enumerate(individual_gene_counts, 1):
            print(f"  Model {i}: {count} genes")
        print(f"Overlapping Genes: {len(overlap_genes)}")
        if overlap_genes:
            print(f"Overlapping Genes List: {', '.join(overlap_genes)}")
            for gene in overlap_genes:
                print(f"\nGene: {gene}")
                if gene in warning_gene.keys():
                    print(f"Warning: possible conflicts for {warning_gene[gene]}, keeping only inhibitor.")
                for desc in descriptions[gene]:
                    print(f"  {desc[0]} Function: {desc[1]}")
        else:
            print("No overlapping genes found.")
     
    return merged_network