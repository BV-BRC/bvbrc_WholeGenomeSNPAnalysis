import json
import os
import re
import sys
import shutil

# Ensure files adhere to the rules defined in the the chewbbaca allele call function

def chewbbaca_filename_format(filename, index):
    # Rule 1 Needs an unique prefix
    # Rule 2: Replace spaces and illegal characters with underscores
    name, ext = os.path.splitext(filename)
    name = re.sub(r'[^A-Za-z0-9_\-]', '_', name)
    new_filename = "{}_{}{}".format(index, name, ext)
    return new_filename

def copy_new_file(clean_fasta_dir, new_name, filename, original_path):
    # deal with moving the files 
    clean_path = os.path.join(clean_fasta_dir, new_name)
    # If the filename was changed, copy the renamed file to the output directory
    if filename != new_name:
        print("Renaming and copying: {} -> {}".format(filename, new_name))
        shutil.copy2(original_path, clean_path)
    else:
        print("Copying: {}".format(filename))
        shutil.copy2(original_path, clean_path)

def fasta_prep(raw_fasta_dir, clean_fasta_dir, config_file):
    # check analysis type
    analysis_type = config_file["params"]["analysis_type"]

    for index, filename in enumerate(sorted(os.listdir(raw_fasta_dir)), start=1):
        original_path = os.path.join(raw_fasta_dir, filename)
        
        # check analysis type
        if analysis_type == "Whole Genome SNP Analysis":
            new_name = ksnp4_filename_format(filename)
            
        elif analysis_type == "chewbbaca":
            new_name = chewbbaca_filename_format(filename, index)
        else:
            print("Invalid analysis type passed. {} \n Exiting Job.".format(config_file["params"]["analysis_type"]))
            sys.exit(1)
        # copy the file to the clean dir 
        copy_new_file(clean_fasta_dir, new_name, filename, original_path)


# Ensure files adhere to the rules defined in the kSNP4 documentation
# All patric ids contain a . that needs to be changed to a _
def ksnp4_filename_format(filename):
    # Update the filename according to kSNP4.1 rules.
    name, ext = os.path.splitext(filename)
    # Rule 1: Remove extra dots in the name (keep only one before the extension)
    name = name.replace(".", "_")

    # Rule 2: Replace spaces and illegal characters with underscores
    name = re.sub(r'[^A-Za-z0-9_\-]', '_', name)

    # Ensure we return a properly formatted filename
    return "{}{}".format(name, ext)

if __name__ == "__main__":
    # check inputs are correct
    if len(sys.argv) != 4:
        print("Usage: python fix_filenames.py <raw_fasta_dir> <clean_fasta_dir> <config_file>")
        sys.exit(1)

    raw_fasta_dir = sys.argv[1]
    clean_fasta_dir = sys.argv[2]
    config_file = sys.argv[3]

    # Not checking that directories exist because Snakemake handles
    # Load config
    try:
        with open(config_file, "r") as f:
            config_file = json.load(f)
    except FileNotFoundError:
        sys.stderr.write(f"Error: {config_file} not found.\n")
        sys.exit(1)
    except json.JSONDecodeError as e:
        sys.stderr.write(f"JSON parse error in {config_file}: {e}\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"Unexpected error: {e}\n")
        sys.exit(1)

    fasta_prep(raw_fasta_dir, clean_fasta_dir, config_file)
    msg = "Checkpoint 4: Filename review and file copy complete."
    sys.stderr.write(msg)

