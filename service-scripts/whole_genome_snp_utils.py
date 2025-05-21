from Bio import Phylo

import click
import json
import numpy as np
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

def add_to_report_dict(source_name, item):
    if source_name not in report_data:
        report_data[source_name] = []  # Initialize the list if the source doesn't exist
    report_data[source_name].append(item)


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

def create_metadata_table(metadata_json, tsv_out, html_output_path):
    with open(metadata_json) as f:
        json_data = json.load(f)

    # # df = pd.read_json(metadata)
    df = pd.json_normalize(json_data)
    # Put columns into an order
    df.columns = ['genome_id', 'genome_name','genbank_accessions',   'genome_status',
       'strain', 'isolation_country', 'collection_year', 'geographic_group',
       'state_province', 'host_common_name', 'host_group',
       'host_scientific_name']
    df.to_csv(tsv_out, index=False, sep="\t")
    table_json_data= tsv_to_html(df, json_data, html_output_path)
    return table_json_data

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

# def infer_ree_type(filename):
#     if "core_SNPs" in filename:
#         return "Core_SNPs"
#     elif "majority" in filename:
#         return "Majority_SNPs"
#     elif "SNPs_all" in filename:
#         return "All_SNPs"

def organize_files_by_type(work_dir, destination_dir):
    print(destination_dir)
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
        if firstword == "All" or filename.startswith("SNPs_all") or filename == "all_snp_distance_heatmap.html":
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
        if firstword == "core" or firstword == "nonCore" or filename == "core_snp_distance_heatmap.html":
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
        if filename == "tree_AlleleCounts.parsimony.tre" or filename == "tree_AlleleCounts.parsimony.tre.phyloxml":
            # print("found".format(file_path))
            pass
        else:
            if group == "All_SNPs" or group == "Core_SNPs":
                tree_dir = os.path.join(destination_dir, group, "Trees")
                if filename == "tree_tipAlleleCounts.SNPs_all.parsimony.tre":
                    group_dir = os.path.join(destination_dir, "All_SNPs")
                    shutil.copy(file_path, group_dir)
                if filename == "tree_tipAlleleCounts.core_SNPs.parsimony.tre":
                    shutil.copy(file_path, "Core_SNPs")
                if filename == "tree_AlleleCounts.SNPs_all.parsimony.NodeLabel.tre":
                    group_dir = os.path.join(destination_dir, "All_SNPs")
                    shutil.copy(file_path, group_dir)
                if filename == "tree_AlleleCounts.core_SNPs.parsimony.NodeLabel.tre":
                    shutil.copy(file_path, "Core_SNPs")
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


def snp_distance_heatmap(config, ksnp_dist_report, output_html):
    # Set up data
    df = pd.read_csv(ksnp_dist_report, sep='\t', header=None)
    # df.columns = ['value', 'genome1', 'genome2']
    df.columns = ["value", "genome1", "genome2"]
    pivoted = df.pivot(index="genome1", columns="genome2", values="value")

   # Define zmin/zmax and map color stops using SNP value thresholds
    zmin, zmax = 0, df["value"].max()

    # Normalize real values to 0–1 range for color scale
    def normalize(v): return (v - zmin) / (zmax - zmin)

    color_scale = [
        [normalize(0), "crimson"],
        [normalize(10), "yellow"],
        [normalize(40), "mistyrose"],
        [1.0, "white"]
    ]

    # Plot heatmap
    fig = px.imshow(pivoted,
                    labels=dict(x="Genome", y="Genome", color="SNP Distance"),
                    x=pivoted.columns,
                    y=pivoted.index,
                    color_continuous_scale=color_scale,
                    zmin=zmin,
                    zmax=zmax,
                    aspect='auto')

    # Manually add scatter traces for the legend
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

def tsv_to_html(df, json_data, html_output_path):
    """
    Convert a DataFrame (from TSV) to an interactive HTML table.

    Args:
        df (pd.DataFrame): The DataFrame containing the TSV data.
        json_data (list of dict): List of metadata dictionaries for the table.
        html_output_path (str): Path to save the generated HTML file.
    """
    try:
        # Rename DataFrame columns for better readability
        df.rename(
            columns={
                "structure_link_html": "Viewer",
                "genome_id": "Genome ID",
                "genome_name": "Genome Name",
                "genbank_accessions": "Genbank Accessions",
                "genome_status": "Genome Status",
                "strain": "Strain",
                "isolation_country": "Isolation Country",
                "collection_year": "Collection Year",
                "geographic_group": "Geographic Group",
                "state_province": "State Province",
                "host_common_name": "Host Common Name",
                "host_group": "Host Group",
                "host_scientific_name": "Host Scientific Name"
            },
            inplace=True
        )

        # Get all unique headers from the JSON data
        all_headers = sorted({key for row in json_data for key in row.keys()})

        # Missing data filled with N/As
        for row in json_data:
            for header in all_headers:
                row.setdefault(header, "N/A")

        # Convert the JSON data to a string for embedding
        json_data_str = json.dumps(json_data)

        # HTML Template with DataTables
        html_template = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Explore All SNP Distances</title>
            <link rel="stylesheet" href="https://cdn.datatables.net/1.13.4/css/jquery.dataTables.min.css">
            <style>
                body {{
                    font-family: Arial, sans-serif;
                    margin: 20px;
                }}
                table {{
                    width: 100%;
                }}
                th input {{
                    width: 100%;
                    box-sizing: border-box;
                }}
            </style>
        </head>
        <body>
            <h1>Interactive Table Viewer</h1>
            <table id="dataTable" class="display" style="width:100%">
                <thead>
                    <tr>
                        {"".join(f"<th>{header}</th>" for header in all_headers)}
                    </tr>
                </thead>
                <tbody>
                </tbody>
            </table>

            <!-- Embedded JSON Data -->
            <script>
                const tableData = {json_data_str};

                document.addEventListener('DOMContentLoaded', function () {{
                    const tableBody = document.querySelector('#dataTable tbody');

                    // Populate the rows
                    tableData.forEach(row => {{
                        const tr = document.createElement('tr');
                        {''.join(f"tr.appendChild(document.createElement('td')).textContent = row['{header}'];" for header in all_headers)}
                        tableBody.appendChild(tr);
                    }});

                    // Initialize DataTables
                    $('#dataTable').DataTable({{
                        pageLength: 10,
                        lengthMenu: [10, 25, 50, 100]
                    }});
                }});
            </script>

            <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
            <script src="https://cdn.datatables.net/1.13.4/js/jquery.dataTables.min.js"></script>
        </body>
        </html>
        """

        # Write the HTML output to the specified file path
        with open(html_output_path, "w") as file:
            file.write(html_template)

        print(f"HTML file successfully created at {html_output_path}")
        return json_data
    except Exception as e:
        print(f"Error: {e}")

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
        if "tree" in filename:
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
@click.argument("config")
@click.argument("ksnp_dist_report")
@click.argument("output_html")
@click.argument("metadata_json")
def write_snp_distance_heatmap(config, ksnp_dist_report, output_html, metadata_json):
    """ Write a heatmap displaying snp level differences between genomes using the output from kSNPdist. kSNPDist only works with the SNPs all matrix fasta file. This command uses the kSNPdist.report."""
    snp_distance_heatmap(config, ksnp_dist_report, output_html)

@cli.command()
@click.argument("metadata_json")
@click.argument("tsv_out")
@click.argument("html_output_path")
def create_metadata_table_output(metadata_json, tsv_out, html_output_path):
    """ Write an interactive view and flat CSV file with all the metadata about this genome gorup. """
    table_json_data = create_metadata_table(metadata, tsv_out, html_output_path)
    trying_together(table_json_data)



if __name__ == "__main__":
    cli()