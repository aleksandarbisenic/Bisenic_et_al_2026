import re
import os
import requests
import pandas as pd

# Cache for fetched modules
module_cache = {}


def validate_parentheses_balanced(s: str) -> bool:
    """Return True if parentheses in s are balanced."""
    stack = 0
    for ch in s:
        if ch == "(":
            stack += 1
        elif ch == ")":
            stack -= 1
            if stack < 0:
                return False
    return stack == 0


def balance_parentheses(tokens):
    """
    Takes a list of whitespace-split tokens and merges adjacent tokens until parentheses are balanced.
    This reproduces your old behavior: one "step" per merged token group.
    """
    steps = []
    buf = ""

    for tok in tokens:
        tok = str(tok).strip()
        if not tok:
            continue

        buf += tok
        if validate_parentheses_balanced(buf):
            steps.append(buf.strip())
            buf = ""
        else:
            buf += " "  # keep merging (this reproduces inner spaces like "((...) K00927,...)")

    if buf.strip():
        # Keep remainder rather than discarding
        steps.append(buf.strip())

    return steps


def convert_ko_to_boolean(expression: str, ko_terms_in_genome: dict) -> str:
    """Replace each K##### with True/False depending on presence in KO dictionary."""
    for ko in re.findall(r"K\d{5}", expression):
        expression = expression.replace(ko, str(ko_terms_in_genome.get(ko, False)))
    return expression


def resolve_or_stretch(expression: str) -> str:
    """Resolve OR stretches: comma-separated True/False values."""
    def eval_or(match):
        terms = match.group(0).split(",")
        return "True" if any(t.strip() == "True" for t in terms) else "False"

    return re.sub(r"(True|False)(,\s*(True|False))+", eval_or, expression)


def resolve_and_stretch(expression: str) -> str:
    """Resolve AND stretches: space-separated True/False values."""
    def eval_and(match):
        terms = match.group(0).split()
        return "True" if all(t.strip() == "True" for t in terms) else "False"

    return re.sub(r"(True|False)(\s+(True|False))+", eval_and, expression)


def eliminate_single_parentheses(expression: str) -> str:
    """Remove parentheses around a single True/False."""
    return re.sub(r"\((True|False)\)", r"\1", expression)


def process_plus_minus_step(step: str, ko_terms_in_genome: dict, max_iter: int = 2000) -> str:
    """
    Special handling for steps containing '+' or '-'.
    Includes a safety guard to prevent infinite loops.
    """
    step = convert_ko_to_boolean(step, ko_terms_in_genome)

    if step.strip() in ["-True", "-False"]:
        step = step.replace("-", "")

    prev = None
    for _ in range(max_iter):
        if step == prev:
            return "False"  # no progress possible -> fail safely
        prev = step

        # Remove any -True/-False fragments
        step = re.sub(r"-True|-False", "", step)

        # Turn True+False+True ... into (True False True) so it becomes an AND-stretch
        def resolve_plus(match):
            terms = match.group(0).split("+")
            return "(" + " ".join(terms) + ")"

        step = re.sub(r"(True|False)(\+(True|False))+", resolve_plus, step)

        step = resolve_or_stretch(step)
        step = resolve_and_stretch(step)
        step = eliminate_single_parentheses(step)

        if step.strip() in ["True", "False"]:
            return step.strip()

    return "False"


def resolve_module_steps(module_steps, ko_terms_in_genome: dict, max_iter: int = 2000):
    total_steps = 0
    fulfilled_steps = 0

    for step in module_steps:
        if not isinstance(step, str):
            continue
        step = step.strip()
        if not step:
            continue

        # Skip these lines entirely (don’t count them)
        if "--" in step:
            continue

        # If the step contains M##### (module reference), count as a step but treat as False
        if re.search(r"M\d{5}", step):
            total_steps += 1
            continue

        total_steps += 1

        if "+" in step or "-" in step:
            step_eval = process_plus_minus_step(step, ko_terms_in_genome, max_iter=max_iter)
        else:
            step_eval = convert_ko_to_boolean(step, ko_terms_in_genome)

            prev = None
            for _ in range(max_iter):
                if step_eval == prev:
                    step_eval = "False"
                    break
                prev = step_eval

                step_eval = resolve_or_stretch(step_eval)
                step_eval = resolve_and_stretch(step_eval)
                step_eval = eliminate_single_parentheses(step_eval)

                if step_eval.strip() in ["True", "False"]:
                    break

        if str(step_eval).strip("()") == "True":
            fulfilled_steps += 1

    return total_steps, fulfilled_steps


def fetch_kegg_module_info_api(module_id: str):
    """
    Correct KEGG flat-file parsing:
    - Read NAME
    - Read only DEFINITION and its continuation lines (field area empty)
    - Stop DEFINITION when the next field begins (e.g., DIAGRAM/CLASS/ORTHOLOGY/...)
    - Then split DEFINITION into whitespace tokens and balance parentheses to recover steps
    """
    if module_id in module_cache:
        return module_cache[module_id]

    url = f"https://rest.kegg.jp/get/{module_id}"

    try:
        resp = requests.get(url, timeout=60)
    except requests.RequestException:
        return None, None

    if resp.status_code != 200:
        return None, None

    module_name = None
    definition_lines = []
    in_definition = False

    for line in resp.text.splitlines():
        field = line[:12].strip() if len(line) >= 12 else line.strip()
        value = line[12:].rstrip() if len(line) >= 12 else ""

        if field == "NAME":
            module_name = value.strip()

        elif field == "DEFINITION":
            in_definition = True
            definition_lines.append(value.strip())

        elif in_definition:
            # Continuation lines have empty field region (first 12 chars spaces)
            if field == "":
                definition_lines.append(value.strip())
            else:
                # Next field begins (DIAGRAM/CLASS/ORTHOLOGY/...) -> stop
                break

    definition_text = " ".join(definition_lines).strip()
    if not module_name or not definition_text:
        module_cache[module_id] = (module_name, [])
        return module_name, []

    # IMPORTANT: Restore original behavior: split DEFINITION into tokens -> steps
    tokens = re.split(r"\s+", definition_text)
    steps = balance_parentheses(tokens)

    module_cache[module_id] = (module_name, steps)
    return module_name, steps


def fetch_modules_in_batches(start_id: int, batch_size: int = 10):
    fetched_modules = []
    end_id = start_id + batch_size

    for module_id in range(start_id, end_id):
        module_id_str = f"M{module_id:05d}"
        module_name, steps = fetch_kegg_module_info_api(module_id_str)

        if module_name and steps:
            fetched_modules.append((module_name, steps))

    print(f"Fetched {len(fetched_modules)} modules in batch starting with ID {start_id}")
    return fetched_modules


def append_modules_to_csv(fetched_modules, out_csv="kegg_module_ko_terms.csv", max_steps=30):
    if os.path.exists(out_csv):
        df = pd.read_csv(out_csv)
    else:
        df = pd.DataFrame()

    current_columns = df.shape[1]

    for module_name, steps in fetched_modules:
        if len(steps) > max_steps:
            print(
                f"Module '{module_name}' has more than {max_steps} steps ({len(steps)} steps). "
                f"Consider increasing the limit."
            )

        truncated_steps = steps[:max_steps]
        padded_col = [module_name] + truncated_steps + [""] * (max_steps - len(truncated_steps))
        df[f"Module_{current_columns}"] = padded_col
        current_columns += 1

    df.to_csv(out_csv, index=False)
    print(f"Appended {len(fetched_modules)} modules to {out_csv}")


def append_ko_completeness_to_csv(ko_files, fetched_modules):
    binary_file = "kegg_module_completeness_binary.csv"
    percentage_file = "kegg_module_completeness_percentage.csv"

    if os.path.exists(binary_file):
        binary_df = pd.read_csv(binary_file)
        percentage_df = pd.read_csv(percentage_file)
    else:
        binary_df = pd.DataFrame(columns=["Bacterium"])
        percentage_df = pd.DataFrame(columns=["Bacterium"])

    # Ensure rows exist for each bacterium
    for ko_file in ko_files:
        bacterium_name = os.path.splitext(os.path.basename(ko_file))[0]
        if bacterium_name not in set(binary_df["Bacterium"].values):
            binary_df = pd.concat([binary_df, pd.DataFrame({"Bacterium": [bacterium_name]})], ignore_index=True)
            percentage_df = pd.concat([percentage_df, pd.DataFrame({"Bacterium": [bacterium_name]})], ignore_index=True)

    # Ensure columns exist for each module (by module name)
    for module_name, _ in fetched_modules:
        if module_name not in binary_df.columns:
            binary_df[module_name] = ""
            percentage_df[module_name] = ""

    # Fill values
    for ko_file in ko_files:
        bacterium_name = os.path.splitext(os.path.basename(ko_file))[0]

        with open(ko_file, "r") as f:
            ko_terms = {line.strip(): True for line in f if line.strip()}

        row_index = binary_df.index[binary_df["Bacterium"] == bacterium_name].tolist()[0]

        for module_name, module_steps in fetched_modules:
            total_steps, fulfilled_steps = resolve_module_steps(module_steps, ko_terms)

            binary_value = 1 if (
                fulfilled_steps == total_steps or (total_steps >= 3 and fulfilled_steps == total_steps - 1)
            ) else 0

            percentage_value = round((fulfilled_steps / total_steps) * 100, 2) if total_steps > 0 else 0

            binary_df.at[row_index, module_name] = binary_value
            percentage_df.at[row_index, module_name] = percentage_value

    binary_df.to_csv(binary_file, index=False)
    percentage_df.to_csv(percentage_file, index=False)
    print(f"Appended completeness results for {len(fetched_modules)} modules to binary and percentage CSVs")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Assess KEGG Module completeness based on KO terms.")
    parser.add_argument("ko_files", metavar="F", type=str, nargs="+", help="KO term files for bacterial genomes")
    args = parser.parse_args()

    batch_size = 10
    module_id = 1

    # Adjust this upper bound if you want to scan more module IDs
    while module_id <= 1000:
        fetched_modules = fetch_modules_in_batches(module_id, batch_size)

        if fetched_modules:
            append_modules_to_csv(fetched_modules)
            append_ko_completeness_to_csv(args.ko_files, fetched_modules)

        module_id += batch_size
