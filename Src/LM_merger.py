import re
import libsbml
import os
import boolean
# from rpy2.robjects.packages import importr
# boolnet = importr("BoolNet")

def read_network(file_path, body_separator=","):
    """
    Reads a Boolean network model from a text file or an SBML-qual file and returns a dictionary representation.
    The Boolean expressions are cleaned of unnecessary parentheses.
    
    :param file_path: Path to the model file (text file or SBML-qual file).
    :param body_separator: Character used to separate targets and factors in the text file.
    :return: Dictionary with gene names as keys and their cleaned Boolean expressions as values.
    """
    # Determine the file type based on the extension
    file_extension = os.path.splitext(file_path)[1].lower()
    
    if file_extension in [".txt", ".csv"]:
        return read_text_network(file_path, body_separator)
    elif file_extension in [".xml", ".sbml"]:
        return read_sbml_qual_network(file_path)
    else:
        raise ValueError("Unsupported file format. Please provide a .txt, .csv, .xml, or .sbml file.")

def read_text_network(file_path, body_separator=","):
    """
    Reads a Boolean network model from a text file and returns a dictionary representation.
    The Boolean expressions are cleaned of unnecessary parentheses.
    
    :param file_path: Path to the model file.
    :param body_separator: Character used to separate targets and factors in the file.
    :return: Dictionary with gene names as keys and their cleaned Boolean expressions as values.
    """
    network = {}

    with open(file_path, 'r') as file:
        lines = file.readlines()
        
    # Strip comments and empty lines
    lines = [re.sub(r'#.*', '', line).strip() for line in lines if line.strip() and not line.startswith("#")]
    
    # Check if the file is empty after stripping comments and empty lines
    if not lines:
        raise ValueError("The file is empty or contains only comments.")
    
    # Check and parse the header
    header = [part.strip() for part in lines[0].split(body_separator)]
    if len(header) < 2 or header[0].lower() != "targets" or header[1].lower() not in ["functions", "factors"]:
        raise ValueError(f"Invalid header: {lines[0]}")
    
    # Process each rule
    for line in lines[1:]:
        # Split the line based on the body_separator, but only at the first occurrence
        parts = [part.strip() for part in line.split(body_separator, 1)]
        if len(parts) < 2:
            raise ValueError(f"Invalid rule format: {line}")
        
        target = parts[0].strip()
        expression = parts[1].strip()
        
        # Validate gene names
        if not re.match(r'^[a-zA-Z_][a-zA-Z0-9_]*$', target):
            raise ValueError(f"Invalid gene name: {target}")
        
        network[target] = simplify_expression(expression)
    
    # Capitalize all the genes and expressions to make the networks case-insensitive
    network = {gene.upper(): expression.upper() for gene, expression in network.items()}
    
    return network

def read_sbml_qual_network(file_path):
    """
    Reads a Boolean network model from an SBML-qual file and returns a dictionary representation.
    The Boolean expressions are cleaned of unnecessary parentheses.
    
    :param file_path: Path to the SBML-qual file.
    :return: Dictionary with gene names as keys and their cleaned Boolean expressions as values.
    """
    reader = libsbml.SBMLReader()
    document = reader.readSBML(file_path)
    model = document.getModel()

    # Parse qualitative species (genes)
    genes = {}
    listOfSpecies = model.getPlugin("qual").getListOfQualitativeSpecies()

    for i in range(listOfSpecies.size()):
        species = listOfSpecies.get(i)
        gene_name = species.getId()
        genes[gene_name] = species.getName() or gene_name

    # Parse transitions (interactions)
    network = {}
    listOfTransitions = model.getPlugin("qual").getListOfTransitions()

    for i in range(listOfTransitions.size()):
        transition = listOfTransitions.get(i)
        output_list = transition.getListOfOutputs()
        function_terms = transition.getListOfFunctionTerms()

        # Get function terms and convert to Boolean expressions
        expression = []
        for j in range(function_terms.size()):
            function_term = function_terms.get(j)
            math = function_term.getMath()
            expression.append(parse_mathml(math, genes))
        # print(expression)
        if len(expression) > 1:
            expression = " | ".join(expression)
        else:
            expression = expression[0]

        for output in output_list:
            gene = genes[output.getQualitativeSpecies()]
            network[gene] = simplify_expression(expression)

    # Capitalize all the genes and expressions to make the networks case-insensitive
    network = {gene.upper(): expression.upper() for gene, expression in network.items()}
    return network


def parse_mathml(math, genes):
    """
    Parses MathML expressions to Boolean expressions.
    
    :param math: MathML object from libSBML.
    :param genes: Dictionary of gene names.
    :return: Boolean expression as a string.
    """
    name = math.getType()
    
    if name == libsbml.AST_CONSTANT_TRUE:
        return "1"
    elif name == libsbml.AST_CONSTANT_FALSE:
        return "0"
    elif name == libsbml.AST_NAME:
        return genes.get(math.getName(), math.getName())
    elif name == libsbml.AST_REAL:
        return str(math.getReal())
    elif name == libsbml.AST_INTEGER:
        return str(math.getInteger())
    elif name == libsbml.AST_RELATIONAL_EQ:
        left_child = parse_mathml(math.getChild(0), genes)
        right_child = parse_mathml(math.getChild(1), genes)
        if right_child == "1":
            return left_child  # GATA1 == 1 should return GATA1
        elif right_child == "0":
            return f"!{left_child}"  # GATA1 == 0 should return !GATA1
        else:
            return f"({left_child} == {right_child})"
    elif name == libsbml.AST_LOGICAL_AND:
        children = [parse_mathml(math.getChild(i), genes) for i in range(math.getNumChildren())]
        return f"({' & '.join(children)})"
    elif name == libsbml.AST_LOGICAL_OR:
        children = [parse_mathml(math.getChild(i), genes) for i in range(math.getNumChildren())]
        return f"({' | '.join(children)})"
    elif name == libsbml.AST_LOGICAL_NOT:
        child = parse_mathml(math.getChild(0), genes)
        if math.getChild(0).getType() in (libsbml.AST_LOGICAL_AND, libsbml.AST_LOGICAL_OR):
            return f"!({child})"
        else:
            return f"!{child}"
    else:
        raise ValueError(f"Unsupported MathML element: {libsbml.formulaToString(math)}")


def simplify_expression(expression):
    """
    Simplifies the Boolean expression by removing redundant parentheses using boolean.py package.
    
    :param expression: The Boolean expression as a string.
    :return: The simplified expression.
    """
    algebra = boolean.BooleanAlgebra()
    algebra.parse(expression)
    simplified_expression = str(algebra.parse(expression).simplify())
    # Replace '~' with '!'
    simplified_expression = simplified_expression.replace('~', '!')

    return simplified_expression


def parse_expression(expression):
    """
    Parses a Boolean expression to separate activators and inhibitors.
    """
    activators = []
    inhibitors = []

    # Detect complex inhibitors like !(GATA1 & GATA2)
    complex_inhibitors = re.findall(r'!\(([^)]+)\)', expression)
    for complex_inhibitor in complex_inhibitors:
        inhibitors.append(complex_inhibitor.strip())
        expression = expression.replace(f"!({complex_inhibitor.strip()})", '')

    # Detect simple inhibitors like !GATA1
    simple_inhibitors = re.findall(r'![a-zA-Z_][a-zA-Z0-9_]*', expression)
    inhibitors.extend([inhibitor[1:] for inhibitor in simple_inhibitors])

    # Remove all inhibitors from the expression to identify activators
    for inhibitor in inhibitors:
        expression = expression.replace(f"&!{inhibitor}", '').replace(f"|!{inhibitor}", '').replace(f"!{inhibitor}", '')

    # Detect activator groups connected by '&' or within parentheses
    activator_groups = re.findall(r'\([^()]*\)|[^()]*&[^()]*', expression)
    # Remove the symbol & behind the last activator name and filter out empty strings and exactly '(|)'
    activator_groups = [group.rstrip('&').strip() for group in activator_groups if group.strip() != '(|)']
    
    for group in activator_groups:
        if group.strip():
            activators.append(group.strip())
        expression = expression.replace(group, '')

    # Detect remaining simple activators
    activators.extend(re.findall(r'[a-zA-Z_][a-zA-Z0-9_]*', expression))

    # Remove empty parentheses from the activators list
    activators = [activator for activator in activators if activator != '()']
    # Remove symbol if it is the first character
    activators = [activator[1:] if activator[0] in ['&', '|'] else activator for activator in activators]
    # Remove symbol and the parentheses if it is the first character within parentheses
    activators = [activator[2:-1] if activator.startswith('(&') or activator.startswith('(|') else activator for activator in activators]
    # Remove duplicates
    activators = list(set(activators))
    inhibitors = list(set(inhibitors))
    return activators, inhibitors

def check_inhibitor_wins_rule(expression):
    """
    Ensures that the "Inhibitor Wins" rule is respected within a single expression.
    """
    components = re.split(r'\|', expression)
    
    for component in components:
        component = component.strip()
        if '&' in component:
            activators, inhibitors = parse_expression(component)
            if activators and inhibitors:
                continue
        if "!" not in component:
            continue
        activators, inhibitors = parse_expression(component)
        if activators and inhibitors:
            return False
    
    return True

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

def write_network_to_file(network, filename, format="text"):
    """
    Writes the merged Boolean network to a file in the specified format.
    
    :param network: Dictionary representing the merged Boolean network.
    :param filename: The name of the file to write the network to.
    :param format: The format of the output file, either "text" or "sbml".
    """
    if format == "text":
        with open(filename + ".txt", 'w') as file:
            file.write("targets, factors\n")
            for gene, expression in network.items():
                file.write(f"{gene}, {expression}\n")
        print(f"Network successfully written to {filename}.txt")
    
    elif format == "sbml":
        # Write to a temporary text file
        temp_text_file = "temp_network.txt"
        with open(temp_text_file, 'w') as file:
            file.write("targets, factors\n")
            for gene, expression in network.items():
                file.write(f"{gene}, {expression}\n")
        
        # Load the network from the text file using boolnet
        net = boolnet.loadNetwork(temp_text_file)
        
        # Export to SBML format
        boolnet.toSBML(net, filename + ".sbml")
        
        # Delete the temporary text file
        os.remove(temp_text_file)
        
        print(f"Network successfully written to {filename}.sbml")
    
    else:
        raise ValueError("Invalid format specified. Use 'text' or 'sbml'.")

def customize_node(merged_network, networks, nodes, method):
    """
    Customize the merge method for specific nodes after merge_networks.
    For the specified nodes, re-merge their rules from the original networks using the given method.
    The rest of the merged_network remains unchanged.

    :param merged_network: The merged network dictionary (output of merge_networks).
    :param networks: List of original network dictionaries.
    :param nodes: List of node names (or a single node name) to customize.
    :param method: Merge method to use for these nodes ("OR", "AND", or "Inhibitor Wins").
    :return: The updated merged_network dictionary.
    """
    if isinstance(nodes, str):
        nodes = [nodes]
    nodes = [n.upper() for n in nodes]
    if method not in ["OR", "AND", "Inhibitor Wins"]:
        raise ValueError("Invalid method. Use 'OR', 'AND', or 'Inhibitor Wins'.")

    algebra = boolean.BooleanAlgebra()

    for node in nodes:
        print(f"Customizing node: {node}")
        # Collect all expressions for this node from all networks
        node_expressions = [network[node] for network in networks if node in network]
        if not node_expressions:
            continue  # Node not present in any network
        # Merge expressions for this node using the specified method
        merged_expr = node_expressions[0]
        for expr in node_expressions[1:]:
            if algebra.parse(merged_expr) == algebra.parse(expr):
                continue
            if method == "OR":
                merged_expr = f"({merged_expr})|({expr})"
            elif method == "AND":
                merged_expr = f"({merged_expr})&({expr})"
            elif method == "Inhibitor Wins":
                activators, inhibitors = parse_expression(expr)
                existing_activators, existing_inhibitors = parse_expression(merged_expr)
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
                if inhibitor_expr:
                    if '&' in inhibitor_expr or '|' in inhibitor_expr:
                        merged_expr = f"({activator_expr})&!({inhibitor_expr})"
                    else:
                        merged_expr = f"({activator_expr})&!{inhibitor_expr}"
                else:
                    merged_expr = activator_expr
        # Simplify and update the merged_network for this node
        print(f"{node} function before customization: {merged_network[node]}")
        merged_expr = simplify_expression(merged_expr)
        print(f"{node} function after customization: {merged_expr}")
        merged_network[node] = simplify_expression(merged_expr)
        print('\n')
    return merged_network
