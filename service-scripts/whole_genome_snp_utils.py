from Bio import Phylo

import click
import json
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.offline as offline
import re
import shutil
import subprocess
import sys

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, leaves_list

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
    # Files coming from the api do not end in fasta. If it does not end in .fasta add extension
    name, ext = os.path.splitext(filename)
    if ext != ".fasta":
        # add_extension = os.path.join(filename, ".fasta")
        add_extension = filename + ".fasta"
        name, ext = os.path.splitext(add_extension)
    # Rule 1: Remove extra dots in the name (keep only one before the extension)
    name = name.replace(".", "_")

    # Rule 2: Replace spaces and illegal characters with underscores
    name = re.sub(r'[^A-Za-z0-9_\-]', '_', name)

    # Ensure we return a properly formatted filename
    return "{}{}".format(name, ext)

def infer_output_type(filename):
    if "core_SNPs" in filename:
        return "Core_SNPs"
    elif "majority" in filename:
        return "Majority_SNPs"
    elif "SNPs_all" in filename:
        return "All_SNPs"

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
        if firstword == "All" or filename.startswith("SNPs_all"):
            All_SNPs_dir = os.path.join(destination_dir, "All_SNPs")
            os.makedirs(All_SNPs_dir, exist_ok=True)
            shutil.copy(file_path, All_SNPs_dir)
        if firstword == "annotate" and os.path.getsize(file_path) > 0:
            shutil.copy(file_path, intermediate_dir)
        if firstword == "ClusterInfo.SNPs" or firstword == "ClusterInfo.core":
            group = infer_output_type(filename)
            cluster_dir = os.path.join(destination_dir, group, "Cluster_Information")
            os.makedirs(cluster_dir, exist_ok=True)
            shutil.copy(file_path, cluster_dir)
        if firstword == "core" or firstword == "nonCore":
            core_snp_dir = os.path.join(destination_dir, "Core_SNPs")
            os.makedirs(core_snp_dir, exist_ok=True)
            shutil.copy(file_path, core_snp_dir)
        if firstword == "COUNT" or firstword == "tip" or firstword == "Node" or firstword == "NJ.dist.matrix":
            intermediate_dir = os.path.join(destination_dir, "Intermediate_Files")
            os.makedirs(intermediate_dir, exist_ok=True)
            shutil.copy(file_path, intermediate_dir)
        if firstword == "Homoplasy":
            group = infer_output_type(filename)
            homoplasy_dir = os.path.join(destination_dir, group, "Homoplasy")
            os.makedirs(homoplasy_dir, exist_ok=True)
            shutil.copy(file_path, homoplasy_dir)
        if "SNPs_in_majority" in filename and "matrix" in filename:
            group = infer_output_type(filename)
            majority_dir = os.path.join(destination_dir, group)
            os.makedirs(majority_dir, exist_ok=True)
            shutil.copy(file_path, majority_dir)
        # Capture specifically SNPs_in_majority0.5 where 5 could be any integer 0-9 without anything else
        if len(filename) == 19 and filename.startswith("SNPs_in_majority0."):
            group = infer_output_type(filename)
            majority_dir = os.path.join(destination_dir, group)
            os.makedirs(majority_dir, exist_ok=True)
            shutil.copy(file_path, majority_dir)
        if firstword.startswith("VCF"):
            VCFs_dir = os.path.join(destination_dir, "VCFs")
            os.makedirs(VCFs_dir, exist_ok=True)
            shutil.copy(file_path, VCFs_dir)
    # Trees get a seperate loop
    clean_tree_dir = os.path.join(work_dir,"clean_trees")
    for filename in os.listdir(clean_tree_dir):
        file_path = os.path.join(clean_tree_dir, filename)
        group = infer_output_type(filename)
        tree_dir = os.path.join(destination_dir, group)
        os.makedirs(tree_dir, exist_ok=True)

        name, ext = os.path.splitext(filename)

        if ext == ".phyloxml":
            shutil.copy(file_path, tree_dir)
        elif ext == ".tre":
            # Newick Trees going in their own subdirectory
            newick_tree_dir = os.path.join(tree_dir, "Newick_Files")
            os.makedirs(newick_tree_dir, exist_ok=True)
            shutil.copy(file_path, newick_tree_dir)
        else:
            msg = "UNKNOWN TREE not uploaded in output dir: {}".format(file_path)
            sys.stderr.write(msg)

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

def snp_distance_heatmap(ksnp_dist_report, output_html):
    # Set up data
    df = pd.read_csv(ksnp_dist_report, sep='\t', header=None)
    # df.columns = ['value', 'genome1', 'genome2']
    df.columns = ["value", "genome1", "genome2"]
    pivoted = df.pivot(index="genome1", columns="genome2", values="value")

   # Define zmin/zmax and map color stops using SNP value thresholds
    zmin, zmax = 0, df["value"].max()

    # Normalize real values to 0–1 range for color scale
    def normalize(v): return (v - zmin) / (zmax - zmin)

    # color_scale = [
    #     [normalize(0), "crimson"],
    #     [normalize(10), "orange"],
    #     [normalize(40), "yellow"],
    #     [1.0, "blue"]
    # ]

    color_scale = [
        [normalize(0), "crimson"],
        [normalize(10), "yellow"],
        [normalize(40), "mistyrose"],
        [1.0, "white"]
    ]

    # # color blind friendly colors
    # color_scale = [
    # [0.0, "#08306b"],    # dark blue (0)
    # [0.05, "#4292c6"],   # medium-light blue (~5)
    # [0.1, "#a6bddb"],    # light blue (~10)
    # [0.25, "#d9d9d9"],   # gray (~25)
    # [0.5, "#f0f0f0"],    # very light gray (~40)
    # [1.0, "#ffffff"]     # white (max)
# ]
    # Plot heatmap
    fig = px.imshow(pivoted,
                    labels=dict(x="Genome", y="Genome", color="SNP Distance"),
                    x=pivoted.columns,
                    y=pivoted.index,
                    color_continuous_scale=color_scale,
                    zmin=zmin,
                    zmax=zmax,
                    aspect='auto')
    # Custom ticks on colorbar legend - covering each other
    # fig.update_layout(
    #     title="SNP Count Heatmap",
    #     coloraxis_colorbar=dict(
    #         title="Criteria for Epidemiological Linkage",
    #         tickvals=[0, 10, 40, zmax],
    #         ticktext=["0-10 (identical)", "10 (low)", "40 (moderate)", f"{zmax} (high)"]
    #     )
    # )
    # Manually add dummy scatter traces for the legend
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(size=10, color='crimson'),
        name='0–10 Strong Linkage'
    ))

    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(size=10, color='yellow'),
        name='10–40 Mid Linkage'
    ))

    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(size=10, color='mistyrose'),
        name='>40 Weak Linkage'
    ))

    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(size=10, color='white'),
        name='Max distance'
    ))
    fig.update_layout(
    legend=dict(
        x=0,        # left (0.0), right (1.0)
        y=0.0,        # top (1.0), bottom (0.0)
        xanchor="left",
        yanchor="bottom",
        bgcolor="rgba(255,255,255,0.6)",  # semi-transparent background
        bordercolor="black",
        borderwidth=1
    ))
    fig.update_layout(title="SNP Count Heatmap")
    fig.update_yaxes(showgrid=True, tickangle=45)
    fig.update_xaxes(showgrid=True, tickangle=45)
    offline.plot(fig, filename=output_html, auto_open=False)
    print('complete')

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
    """Organize files by type based on the first word. Trees are managed via convert_to_phyloxml_trees."""
    with open(service_config) as file:
        data = json.load(file)
    work_dir = data["work_data_dir"]
    destination_dir = data["output_data_dir"]
    organize_files_by_type(work_dir, destination_dir)

@cli.command()
@click.argument("kchooser_report")
def parse_kchooser_report(kchooser_report):
    """ Parse the kChooser report for the optimum K value for the kSNP4 command."""
    parse_optimum_k(kchooser_report)

@cli.command()
@click.argument("ksnp_dist_report")
@click.argument("output_html")
def write_snp_distance_heatmap(ksnp_dist_report):
    """ Write a heatmap displaying snp level differences between genomes using the output from kSNPdist. kSNPDist only works with the SNPs all matrix fasta file. This command uses the kSNPdist.report."""
    snp_distance_heatmap(ksnp_dist_report)

if __name__ == "__main__":
    cli()