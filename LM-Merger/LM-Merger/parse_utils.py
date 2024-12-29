import libsbml
import boolean

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