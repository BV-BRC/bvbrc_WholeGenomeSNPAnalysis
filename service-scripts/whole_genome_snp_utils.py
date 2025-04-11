from Bio import Phylo

import click
import json
import os
import re
import shutil
import subprocess
import sys

def copy_new_file(clean_fasta_dir, new_name, filename, original_path):
    # Deal with moving the files 
    clean_path = os.path.join(clean_fasta_dir, new_name)
    # If the filename was changed, copy the renamed file to the output directory
    if filename != new_name:
        print("Renaming and copying: {} -> {}".format(filename, new_name))
        shutil.copy2(original_path, clean_path)
    else:
        print("Copying: {}".format(filename))
        shutil.copy2(original_path, clean_path) 

def edit_newick_genome_id(raw_nwk, clean_nwk):
    """Reverts genome ids in kSNP4 formatted Newick files to match genome ids for phyloxml"""
    click.echo("Reverting genome IDs {}".format(raw_nwk))
    if os.path.isfile(raw_nwk) == True and os.path.getsize(raw_nwk) > 0:
        # Using bio phylo package to ensure accuracy as some files have extra details
        fix_labels_with_phylo(raw_nwk, clean_nwk)
        return clean_nwk
    else:
        msg = "{} is either empty or not found... cannot edit genome ids for phyloxml".format(raw_nwk)
        sys.stderr.write(msg)

def fix_labels_with_phylo(raw_nwk, clean_nwk):
    tree = Phylo.read(raw_nwk, "newick")
    for clade in tree.find_clades():
        if clade.name:
            clade.name = clade.name.replace("_", ".")
    Phylo.write(tree, clean_nwk, "newick")

def ksnp4_filename_format(filename):
    # Update the filename according to kSNP4.1 rules.
    name, ext = os.path.splitext(filename)
    # Rule 1: Remove extra dots in the name (keep only one before the extension)
    name = name.replace(".", "_")

    # Rule 2: Replace spaces and illegal characters with underscores
    name = re.sub(r'[^A-Za-z0-9_\-]', '_', name)

    # Ensure we return a properly formatted filename
    return "{}{}".format(name, ext)

def organize_files_by_type(work_dir, destination_dir):
    if not os.path.exists(work_dir):
        sys.stderr.write("Work directory, {}, does not exist".format(work_dir))
        return
    for filename in os.listdir(work_dir):
        file_path = os.path.join(work_dir, filename)
        intermediate_dir = os.path.join(destination_dir, "Intermediate_Files")
        os.makedirs(intermediate_dir, exist_ok=True)
        firstword = filename.split("_")[0]

        # Skip directories
        if not os.path.isfile(file_path):
            continue
        
        # Sort files into directories according to the first word
        if firstword == "All" or firstword == "nonCore" or firstword == "SNPs":
            All_SNPs_dir = os.path.join(destination_dir, "All_SNPs")
            os.makedirs(All_SNPs_dir, exist_ok=True)
            shutil.copy(file_path, All_SNPs_dir)
        if firstword == "annotate" and os.path.getsize(file_path) > 0:
            shutil.copy(file_path, intermediate_dir)
        if firstword == "ClusterInfo.SNPs" or firstword == "ClusterInfo.core":
            cluster_dir = os.path.join(destination_dir, "Cluster_Information")
            os.makedirs(cluster_dir, exist_ok=True)
            shutil.copy(file_path, cluster_dir)
        if firstword == "core":
            core_snp_dir = os.path.join(destination_dir, "Core_SNPs")
            os.makedirs(core_snp_dir, exist_ok=True)
            shutil.copy(file_path, core_snp_dir)
            os.makedirs(core_snp_dir, exist_ok=True)
        if firstword == "COUNT" or firstword == "tip" or firstword == "Node" or firstword == "NJ.dist.matrix":
            shutil.copy(file_path, intermediate_dir)
        if firstword == "Homoplasy":
            homoplasy_dir = os.path.join(destination_dir, "Homoplasy")
            os.makedirs(homoplasy_dir, exist_ok=True)
            shutil.copy(file_path, homoplasy_dir)
        if firstword.startswith("VCF"):
            VCFs_dir = os.path.join(destination_dir, "VCFs")
            os.makedirs(VCFs_dir, exist_ok=True)
            shutil.copy(file_path, VCFs_dir)
    # trees get a seperate loop
    clean_tree_dir = os.path.join(work_dir,"clean_trees")
    tree_dir = os.path.join(destination_dir, "Trees")
    os.makedirs(tree_dir, exist_ok=True)
    for filename in os.listdir(clean_tree_dir):
        file_path = os.path.join(clean_tree_dir, filename)
        if os.path.isfile(file_path) == True and os.path.getsize(file_path) > 0:
            shutil.copy(file_path, tree_dir)

def parse_optimum_k(kchooser_report):
    with open(kchooser_report, 'r') as file:
        for line in file:
            match = re.search(r'The optimum value of k is (\d+)', line)
            if match:
                print(match.group(1))
                return
    print("Optimum value of k not found")

def run_newick_to_phyloxml(clean_nwk):
        # Run phyloxml command
        result = subprocess.run(["p3x-newick-to-phyloxml", "--verbose", "-l", "genome_id", "-g", "collection_year,host_common_name,isolation_country,strain,genome_name,genome_id,accession,subtype,lineage,host_group,collection_date,geographic_group,geographic_location", clean_nwk])
        msg = "{}".format(result)
        sys.stderr.write(msg)
        
@click.group()
def cli():
    """ This script supports the Whole Genome SNP service with multiple commands."""
    pass


@cli.command()
@click.argument("service_config")
def clean_fasta_filenames(service_config):
    """Ensure files adhere to the rules defined by kSNP4"""
    with open(service_config) as file:
        data = json.load(file)
        raw_fasta_dir = data["raw_fasta_dir"]
        clean_fasta_dir = data["clean_data_dir"]
        for index, filename in enumerate(sorted(os.listdir(raw_fasta_dir)), start=1):
            original_path = os.path.join(raw_fasta_dir, filename)
            new_name = ksnp4_filename_format(filename)
            copy_new_file(clean_fasta_dir, new_name, filename, original_path)

@cli.command()
@click.argument("service_config")
def convert_to_phyloxml_trees(service_config):
    """Use genome IDs in the tree files for phyloxml to connect the existing metadata. Iterate through each tree file to remove kSNP4 formating restrictions."""
    with open(service_config) as file:
        data = json.load(file)
    ### start organize output files ###
    work_dir = data["work_data_dir"]
    if not os.path.exists(work_dir):
        sys.stderr.write("Work directory, {}, does not exist".format(work_dir))
        return
    for filename in os.listdir(work_dir):
        file_path = os.path.join(work_dir, filename)
        # Skip directories
        if not os.path.isfile(file_path):
            continue
        clean_tree_dir = os.path.join(work_dir, "clean_trees")
        os.makedirs(clean_tree_dir, exist_ok=True)
        firstword = filename.split("_")[0]
        if firstword == "tree" or firstword == "tree.SNPs":
            clean_nwk_path = os.path.join(clean_tree_dir, filename)
            edit_newick_genome_id(file_path, clean_nwk_path)
            run_newick_to_phyloxml(clean_nwk_path)

@cli.command()
@click.argument("service_config")
def organize_output_files(service_config):
    """Organize files by type based on the first word. Trees are managed via convert_to_phyloxml_trees"""
    with open(service_config) as file:
        data = json.load(file)
    work_dir = data["work_data_dir"]
    destination_dir = data["output_data_dir"]
    organize_files_by_type(work_dir, destination_dir)

@cli.command()
@click.argument("kchooser_report")
def parse_kchooser_report(kchooser_report):
    parse_optimum_k(kchooser_report)

if __name__ == "__main__":
    cli()