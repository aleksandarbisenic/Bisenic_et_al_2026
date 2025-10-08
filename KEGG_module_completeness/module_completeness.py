import re
import os
import requests
import csv
import pandas as pd

# Pre-fetch KEGG modules cache
module_cache = {}

# Function to validate and fix parentheses structure in a step
def validate_and_fix_parentheses(step):
    """Validates and ensures parentheses in the step are balanced."""
    stack = []
    for char in step:
        if char == '(':
            stack.append('(')
        elif char == ')':
            if not stack:
                return False
            stack.pop()
    return not stack

# Function to balance parentheses for module steps
def balance_parentheses(steps):
    """Balances parentheses for module steps, ensuring that steps are not split across multiple lines."""
    balanced_steps = []
    temp_step = ""

    for step in steps:
        temp_step += step.strip()
        if validate_and_fix_parentheses(temp_step):
            balanced_steps.append(temp_step)
            temp_step = ""
        else:
            temp_step += " "

    return balanced_steps

# Function to convert KO terms in the expression to True/False based on the KO terms present in the bacterial genome
def convert_ko_to_boolean(expression, ko_terms_in_genome):
    """Converts KO terms in the expression to their corresponding True/False values."""
    for ko in re.findall(r'K\d+', expression):  # Find all KO terms like K00036, K19243 etc.
        expression = expression.replace(ko, str(ko_terms_in_genome.get(ko, False)))
    return expression

# Function to resolve OR (comma-separated) and AND (space-separated) logic
def resolve_or_stretch(expression):
    """Resolve OR stretches, where values are separated by commas."""
    def eval_or(match):
        terms = match.group(0).split(',')
        return 'True' if any(term.strip() == 'True' for term in terms) else 'False'
    
    return re.sub(r'(True|False)(,\s*(True|False))+', eval_or, expression)

def resolve_and_stretch(expression):
    """Resolve AND stretches, where values are separated by spaces."""
    def eval_and(match):
        terms = match.group(0).split()
        return 'True' if all(term.strip() == 'True' for term in terms) else 'False'
    
    return re.sub(r'(True|False)(\s+(True|False))+', eval_and, expression)

# Function to eliminate parentheses around single True/False terms
def eliminate_single_parentheses(expression):
    """Eliminates parentheses that surround only a single word (True or False)."""
    return re.sub(r'\((True|False)\)', r'\1', expression)

# Function to handle module steps with + or - symbols
def process_plus_minus_step(step, ko_terms_in_genome):
    """Processes module steps with + or - symbols according to specified logic."""
    
    # Convert KO terms to True/False based on the presence in ko_terms_in_genome
    step = convert_ko_to_boolean(step, ko_terms_in_genome)

    if step.strip() in ['-True', '-False']:
        step = step.replace('-', '')

    while True:
        # Remove any -True or -False parts
        step = re.sub(r'-True|-False', '', step)

        # Resolve any '+' separated stretches
        def resolve_plus(match):
            terms = match.group(0).split('+')
            return '(' + ' '.join(terms) + ')'

        step = re.sub(r'(True|False)(\+True|\+False)+', resolve_plus, step)

        # Apply OR and AND logic and remove parentheses around single True/False terms
        step = resolve_or_stretch(step)
        step = resolve_and_stretch(step)
        step = eliminate_single_parentheses(step)

        # Exit only if step is exactly "True" or "False" with no other characters
        if step.strip() in ["True", "False"]:
            break

    return step

# Function to resolve module steps
def resolve_module_steps(module_ko_list, ko_terms_in_genome):
    total_steps = 0
    fulfilled_steps = 0

    for step in module_ko_list:
        # Skip step if it contains "--" and do not count it in total_steps
        if '--' in step:
            continue  # Skip this step

        # Check if the step contains "M" followed by 5 digits; if so, consider it False and skip further processing
        if re.search(r'M\d{5}', step):
            total_steps += 1  # Count the step but consider it as False
            continue  # Move to the next step

        total_steps += 1  # Count the step only if it does not contain "--" or "MXXXXX"

        # If the step contains + or -, use the special logic
        if '+' in step or '-' in step:
            step = process_plus_minus_step(step, ko_terms_in_genome)
        else:
            # Use the regular logic for steps without + or -
            step = convert_ko_to_boolean(step, ko_terms_in_genome)

            while True:
                # Resolve OR and AND logic
                step = resolve_or_stretch(step)
                step = resolve_and_stretch(step)
                # Eliminate parentheses around single True/False terms
                step = eliminate_single_parentheses(step)

                # Exit only if step is exactly "True" or "False" with no other characters
                if step.strip() in ["True", "False"]:
                    break

        # Count fulfilled steps
        if step.strip("()") == "True":
            fulfilled_steps += 1

    return total_steps, fulfilled_steps

# Function to fetch and cache KEGG module information via the KEGG REST API
def fetch_kegg_module_info_api(module_id):
    if module_id in module_cache:
        return module_cache[module_id]

    base_url = f"https://rest.kegg.jp/get/{module_id}"
    response = requests.get(base_url)

    if response.status_code != 200:
        return None, None

    lines = response.text.splitlines()

    module_name = None
    steps = []
    collecting_steps = False

    for line in lines:
        if line.startswith("NAME"):
            module_name = " ".join(line.split()[1:])
        elif line.startswith("ORTHOLOGY"):
            break
        elif line.startswith("DEFINITION"):
            definition = line.split("DEFINITION")[1].strip()
            raw_steps = [definition]
            collecting_steps = True
        elif collecting_steps:
            raw_steps.append(line.strip())

    if raw_steps:
        steps = balance_parentheses(re.split(r'\s+', " ".join(raw_steps).strip()))

    module_cache[module_id] = (module_name, steps)
    return module_name, steps

# Function to fetch modules in batches of 10
def fetch_modules_in_batches(start_id, batch_size=10):
    fetched_modules = []
    end_id = start_id + batch_size

    for module_id in range(start_id, end_id):
        module_id_str = f"M{module_id:05d}"
        module_name, steps = fetch_kegg_module_info_api(module_id_str)

        if module_name and steps:
            fetched_modules.append((module_name, steps))
    
    print(f"Fetched {len(fetched_modules)} modules in batch starting with ID {start_id}")

    return fetched_modules

# Function to append modules to CSV in a module-per-column manner with a fixed maximum number of steps
def append_modules_to_csv(fetched_modules):
    max_steps = 30

    if os.path.exists('kegg_module_ko_terms.csv'):
        df = pd.read_csv('kegg_module_ko_terms.csv')
    else:
        df = pd.DataFrame()

    current_columns = df.shape[1]

    for module_name, steps in fetched_modules:
        if len(steps) > max_steps:
            print(f"Module '{module_name}' has more than {max_steps} steps ({len(steps)} steps). Consider increasing the limit.")
        
        truncated_steps = steps[:max_steps]
        padded_steps = [module_name] + truncated_steps + [''] * (max_steps - len(truncated_steps))
        df[f'Module_{current_columns}'] = padded_steps
        current_columns += 1

    df.to_csv('kegg_module_ko_terms.csv', index=False)
    print(f"Appended {len(fetched_modules)} modules to kegg_module_ko_terms.csv")
    del df

# Function to process KO files and append results to CSV with module names as headers
def append_ko_completeness_to_csv(ko_files, fetched_modules):
    binary_file = 'kegg_module_completeness_binary.csv'
    percentage_file = 'kegg_module_completeness_percentage.csv'

    # Load existing data if files exist, otherwise create new DataFrames
    if os.path.exists(binary_file):
        binary_df = pd.read_csv(binary_file)
        percentage_df = pd.read_csv(percentage_file)
    else:
        binary_df = pd.DataFrame(columns=["Bacterium"])
        percentage_df = pd.DataFrame(columns=["Bacterium"])

    # Ensure a row for each KO file in both DataFrames
    for ko_file in ko_files:
        bacterium_name = os.path.splitext(os.path.basename(ko_file))[0]
        if bacterium_name not in binary_df['Bacterium'].values:
            # Add new rows for the bacterium using pd.concat()
            new_row = pd.DataFrame({"Bacterium": [bacterium_name]})
            binary_df = pd.concat([binary_df, new_row], ignore_index=True)
            percentage_df = pd.concat([percentage_df, new_row], ignore_index=True)

    # Add new columns for each module in this batch, using module names as column headers
    for module_name, module_ko_list in fetched_modules:
        if module_name not in binary_df.columns:
            binary_df[module_name] = ''
            percentage_df[module_name] = ''

    # Process each KO file
    for ko_file in ko_files:
        bacterium_name = os.path.splitext(os.path.basename(ko_file))[0]

        # Load KO terms for the current bacterium
        with open(ko_file, 'r') as f:
            ko_terms = {line.strip(): True for line in f}

        # Calculate completeness for each fetched module
        for module_name, module_ko_list in fetched_modules:
            total_steps, fulfilled_steps = resolve_module_steps(module_ko_list, ko_terms)

            binary_value = 1 if (fulfilled_steps == total_steps or (total_steps >= 3 and fulfilled_steps == total_steps - 1)) else 0
            percentage_value = round((fulfilled_steps / total_steps) * 100, 2) if total_steps > 0 else 0

            # Set the completeness values in the appropriate row and column for this bacterium and module
            row_index = binary_df.index[binary_df['Bacterium'] == bacterium_name].tolist()[0]
            binary_df.at[row_index, module_name] = binary_value
            percentage_df.at[row_index, module_name] = percentage_value

    binary_df.to_csv(binary_file, index=False)
    percentage_df.to_csv(percentage_file, index=False)
    print(f"Appended completeness results for {len(fetched_modules)} modules to binary and percentage CSVs")

# Main function to run the batch processing
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Assess KEGG Module completeness based on KO terms.")
    parser.add_argument('ko_files', metavar='F', type=str, nargs='+', help='KO term files for bacterial genomes')

    args = parser.parse_args()

    batch_size = 10
    module_id = 1

    while module_id <= 1000:
    	fetched_modules = fetch_modules_in_batches(module_id, batch_size)

    	if fetched_modules:  # Only process if there are fetched modules
        	append_modules_to_csv(fetched_modules)
        	append_ko_completeness_to_csv(args.ko_files, fetched_modules)

    	module_id += batch_size
