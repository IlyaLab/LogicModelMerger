import os
import re
import libsbml
import boolean
from rpy2.robjects.packages import importr
boolnet = importr("BoolNet")
from parse_utils import parse_mathml, simplify_expression

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