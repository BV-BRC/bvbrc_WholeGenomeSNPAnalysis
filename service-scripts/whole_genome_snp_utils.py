from Bio import Phylo
import click
import json
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import re
import shutil
import subprocess
import sys

from Bio import SeqIO
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

def add_to_report_dict(report_data, source_name, item):
    if source_name not in report_data:
        report_data[source_name] = []  # Initialize the list if the source doesn't exist
    report_data[source_name].append(item)
    return report_data


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

def cluster_heatmap_data(genome_ids, snp_matrix):
    # Convert matrix to distance format
    dist_array = squareform(snp_matrix)
    # linkage_result = linkage(dist_array, method="average")
    linkage_result = linkage(dist_array, method="single")
    idx = leaves_list(linkage_result)
    clustered_matrix = [[snp_matrix[i][j] for j in idx] for i in idx]
    clustered_labels = [genome_ids[i] for i in idx]
    return clustered_labels, clustered_matrix


def create_genome_length_bar_plot(clean_data_dir):
    genome_lengths = []
    for filename in os.listdir(clean_data_dir):
        if filename.endswith(".fasta") or filename.endswith(".fa") or filename.endswith(".fna"):
            file_path = os.path.join(clean_data_dir, filename)
            total_length = sum(len(record.seq) for record in SeqIO.parse(file_path, "fasta"))
            display_name = os.path.splitext(filename)[0].replace("_", ".")
            genome_lengths.append({"Genome": display_name, "Length": total_length})
    # Bar Plot
    fig = px.bar(genome_lengths, 
                 x="Genome", 
                 y="Length", 
                 hover_data=["Genome", "Length"],
                 title="Genome Lengths by Genome")
    fig.update_layout(
        xaxis_title="Genome",
        yaxis_title="Genome Length (bp)",
        xaxis_tickangle=-45
    )
    fig.write_html("genome_length_barplot.html", include_plotlyjs=False)
    barplot_html = read_plotly_html("genome_length_barplot.html")
    return barplot_html


def create_metadata_table(metadata_json, tsv_out):
    with open(metadata_json) as f:
        metadata = json.load(f)
    # Convert genome_id: replace '.' with '_'
    for record in metadata:
        # match the style of the genome ids in the heatmap
        record["genome_id"] = record["genome_id"].replace("_", ".")
    metadata_df = pd.json_normalize(metadata)

    # Get all unique headers from the JSON data
    all_headers = sorted({key for row in metadata for key in row.keys()})

    # Missing data filled with N/As
    for row in metadata:
        for header in all_headers:
            row.setdefault(header, "N/A")
    # export as TSV
    metadata_df.to_csv(tsv_out, index=False, sep="\t")
    return metadata, metadata_df


def define_html_template(input_genome_table, barplot_html, snp_distribution_html, homoplastic_snps_html, heatmap_html, majority_threshold, metadata_json_string):
    majority_percentage = majority_threshold * 100
    html_template = """
            <!DOCTYPE html>
            <html lang="en">
            <head>
                <meta charset="UTF-8">
                <meta name="viewport" content="width=device-width, initial-scale=1.0">
                <title>Whole Genome SNP Analysis Report</title>
                <style>
                    .plot-container {{
                        display: block;
                        text-align: center;
                    }}
                    .plot {{
                        width: 80%;
                        margin: 0 auto;
                    }}
                    body {{ 
                        font-family: Roboto, sans-serif; 
                        color: black; 
                    }}
                    header {{
                        display: flex;
                        justify-content: space-between;
                        align-items: center;
                        padding: 10px 20px;
                    }}
                    header > a img {{
                        max-width: 225px;
                        max-height: 225px;
                        width: auto;
                        height: auto;
                    }}
                    .title {{
                        font-size: 36px;  
                        font-family: 'Roboto', sans-serif;
                        font-weight: bold;
                        color: black;
                    }}
                    .warning {{ 
                        color: black; 
                    }}
                    .heatmap-controls,
                    .linkage-controls {{
                    display: flex;
                    flex-wrap: wrap;
                    gap: 10px;
                    align-items: center;
                    margin-bottom: 1em;
                    }}
                    .heatmap-controls h4,
                    .linkage-controls h4 {{
                    flex-basis: 100%;
                    margin: 0;
                    }}
                    table, th, td {{ 
                        border: 1px solid black; border-collapse: collapse; 
                    }}
                    th, td {{ 
                        padding: 5px; text-align: left; 
                    }}
                    img {{ 
                        width: 100%; max-width: 600px; height: auto; 
                    }}
                    .image-row {{
                        display: flex;
                        flex-wrap: wrap;
                        justify-content: flex-start;
                    }}
                    .image-container {{
                        width: 33%; /* Each image container takes up one-third of the row */
                        padding: 5px; /* Padding around the images */
                        box-sizing: border-box;
                    }}
                    .img {{
                        width: 100%;
                        max-width: 600px;
                    /* DataTables styles */
                    table#dataTable {{
                        width: 100%;
                    }}
                    }}
                    #svgContainer {{
                        width: 100%;
                        overflow-x: auto;
                        margin-top: 8px;
                    }}
                    th input {{
                        width: 100%;
                        box-sizing: border-box;
                    }}
                    table {{
                    width: 100%;
                    }}
                    th input {{
                    width: 100%;
                    box-sizing: border-box;
                    }}
                </style>
                <!-- Plotly.js v3.0.1  — last updated June 2025 --> 
                <script src="https://cdn.plot.ly/plotly-3.0.1.min.js"></script>
                <!-- DataTables CSS -->
                <link rel="stylesheet" href="https://cdn.datatables.net/1.13.4/css/jquery.dataTables.min.css">
                <header>
                    <div class="title">Whole Genome SNP Analysis Report</div>
                        </a>
                </header>
            <p> Explore the results of your Whole Genome SNP analysis in this interactive report. This service uses a reference-free tool for identifying Single Nucleotide Polymorphisms (SNPs) called kSNP4.  This reference free tool identifies SNPs through a k-mer based analysis.  Note, you can download the images of the plots in this report by clicking on the camera icon in the upper left-hand corner. Poor quality assemblies can impact the performance of this service.</p>
            <h3>About the Analysis Workflow</h3>
            <p> 
            kSNP4 identifies SNPs and does phylogenetic analysis without genome alignment or the use of reference genomes. They estimate phylogenetic trees based on 3 methods: parsimony, neighbor joining (NJ), and maximum likelihood (ML). The tool creates three groups, All SNPs, Core SNPs and SNPs within a majority threshold.
             
             Our pipeline begins by selecting the optimum k-mer size by reviewing your genomes with a series of add-length k-mers to seek the shortest length for which each k-mer occurs only once in the median-sized genome.
             
             While kSNP4 runs it identifies SNPs across the three genome groups (All, Core, Majority).  Then constructs phylogenetic trees based on the SNPs. Static images of the trees are available in this report. The trees are also viewable in the website's phylogenetic tree viewer with the ability to map metadata to the tree. Other files available are multiple alignment of all SNP positions and related files. A companion program, kSNPdist creates a SNP distance matrix for each genome group.</p>
            <h2>Getting to Know the Input Data</h2>
            <h3>Input Genomes</h3>
            <p>
            {input_genome_table}
            <p>
            <div class="plot-container">
                <div class="plot" id="Genome Length Box Plot">
                {barplot_html}
                </div>
            </div>
            <h3>Input Genome Metadata</h3>
            <!-- DataTables CSS -->
            <link rel="stylesheet" href="https://cdn.datatables.net/1.13.4/css/jquery.dataTables.min.css">
            <table id="dataTable" class="display" style="width:100%">
                <thead id="tableHead"></thead>
                <tbody id="tableBody"></tbody>
            </table>

            <!-- Embedded JSON Data -->
            <script>
                const tableData = {metadata_json_string};
            </script>

            <!-- jQuery -->
            <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
            <!-- DataTables -->
            <script src="https://cdn.datatables.net/1.13.4/js/jquery.dataTables.min.js"></script>
            <script>
            function populateTable(data) {{
                const tableHead = document.getElementById('tableHead');
                const tableBody = document.getElementById('tableBody');

                // Clear previous content
                tableHead.innerHTML = '';
                tableBody.innerHTML = '';

                // Extract headers from the first row
                const headers = Object.keys(data[0]);

                // Create header row and filter row
                const headRow = document.createElement('tr');
                const filterRow = document.createElement('tr');

                headers.forEach((header, index) => {{
                    // Header row
                    const th = document.createElement('th');
                    th.textContent = header;
                    headRow.appendChild(th);

                    // Filter row
                    const filterCell = document.createElement('th');
                    const filterInput = document.createElement('input');
                    filterInput.type = 'text'; // Accept both numeric and string filters
                    filterInput.placeholder = `Filter ${{header}}`;
                    filterInput.dataset.column = index;
                    filterInput.addEventListener('keyup', function () {{
                        const columnIndex = parseInt(this.dataset.column);
                        const filterValue = this.value.trim();
                        const isNumeric = data.every(row => !isNaN(parseFloat(row[header])));

                        if (filterValue) {{
                            if (isNumeric) {{
                            // Numeric threshold filtering
                                const parsedValue = parseFloat(filterValue);
                                if (!isNaN(parsedValue)) {{
                                    $.fn.dataTable.ext.search = [];
                                    $.fn.dataTable.ext.search.push((settings, row) => {{
                                    const cellValue = parseFloat(row[columnIndex]);
                                    return cellValue >= parsedValue; // Filter values >= threshold
                                    }});
                                }}
                                }} else {{
                                    // String filtering (case-insensitive substring match)
                                    $.fn.dataTable.ext.search = [];
                                    $.fn.dataTable.ext.search.push((settings, row) => {{
                                        const cellValue = row[columnIndex].toLowerCase();
                                        return cellValue.includes(filterValue.toLowerCase());
                                    }});
                                }}
                            }} else {{
                                // Clear filter if the input is empty
                                $.fn.dataTable.ext.search = [];
                            }}

                            // Redraw the table to apply the new filter
                            $('#dataTable').DataTable().draw();
                            }});
                        filterCell.appendChild(filterInput);
                        filterRow.appendChild(filterCell);
                    }});

                    tableHead.appendChild(headRow);
                    tableHead.appendChild(filterRow);

                    // Create data rows
                    data.forEach(row => {{
                        const tr = document.createElement('tr');
                        headers.forEach(header => {{
                        const td = document.createElement('td');
                        td.innerHTML = row[header]; // Render HTML links directly
                        tr.appendChild(td);
                        }});
                        tableBody.appendChild(tr);
                    }});

                    // Initialize DataTables
                    $('#dataTable').DataTable({{
                        pageLength: 10,  // Rows per page
                        lengthMenu: [10, 25, 50, 100],  // Pagination options
                        orderCellsTop: true,  // Keep filters at the top
                        initComplete: function () {{
                        // Automatically focus on the filter input of the first column for usability
                        $('#dataTable thead input').first();
                        }}
                    }});
                }}

                // Populate the table on page load
                document.addEventListener('DOMContentLoaded', function () {{
                populateTable(tableData);
                // Ensure all links open in a new tab
                const links = document.querySelectorAll('#dataTable a');
                links.forEach(link => {{
                    link.target = '_blank';
                }});
            }});
            </script>
            <br>
            <h2>Reviewing Identified SNPs</h2>
            <h3>SNPs Distribution</h3> 
            <div class="plot-container">
                <div class="plot" id="SNP Distribution Bar Plot">
                {snp_distribution_html}
                </div>
            </div>
            <p> This service defines three subsets of SNPs </p>
            <ul style="list-style-type: disc; padding-left: 25;">
            <li><b>Total SNPs</b> includes every SNP identified from all genomes given to the service regardless of how many genomes the SNP is present it.</li>
            <li><b>Core SNPs</b> are present in every genome analyzed.</li>
            <li><b>Majority SNPs</b> are present in at least {majority_percentage} percent of the genomes analyzed.</li>
            </ul>
            <h3>Homoplastic SNPs</h3>
            {homoplastic_snps_html}
            <br>
            <p>A SNP is considered Homoplastic which that SNP occurs in two or more unrelated positions on the same tree.</p>
            <ul style="list-style-type: disc; padding-left: 25;">
            <li><b>Parsimony</b> the parsimony tree method is a good fit for small datasets with close relatives and low divergency. It estimated by creating a consensus of up to 100 equally parsimonious trees. It seeks a tree topology that explains the observed sequence data with the smallest possible number of evolutionary changes.</li>
            <li><b>Maximum Likelihood</b> the maximum likelihood tree is a good fit for datasets with substantial divergency and complex substitution patterns. It is constructed by finding the tree topology that has the highest likelihood of producing the observed sequence data. </li>
            <li><b>Neighbor Joining</b> this tree method is a good fit for exploratory analysis, especially with very large datasets. This is a distance-based method that constructs a tree by iteratively finding pairs of taxa (neighbors) that minimize  the total branch length at each step.  It uses distance matrix (pairwise genetic distances between sequences).</li>
            </ul>
            <p>Please visit the kSNP4 documentation for more information about the many trees created by this service.<p>
            {heatmap_html}
            <h3>References</h3>

            <ol type="1">
            <li>Olson RD, Assaf R, Brettin T, Conrad N, Cucinell C, Davis JJ, Dempsey DM, Dickerman A, Dietrich EM, Kenyon RW, Kuscuoglu M, Lefkowitz EJ, Lu J, Machi D, Macken C, Mao C, Niewiadomska A, Nguyen M, Olsen GJ, Overbeek JC, Parrello B, Parrello V, Porter JS, Pusch GD, Shukla M, Singh I, Stewart L, Tan G, Thomas C, VanOeffelen M, Vonstein V, Wallace ZS, Warren AS, Wattam AR, Xia F, Yoo H, Zhang Y, Zmasek CM, Scheuermann RH, Stevens RL.</li>
            <li>Hall, B. G., & Nisbet, J. (2023). Building phylogenetic trees from genome sequences with KSNP4. Molecular Biology and Evolution, 40(11). https://doi.org/10.1093/molbev/msad235 </li>
            </ol>
        </body>
        </html>
        """.format(input_genome_table = input_genome_table, barplot_html=barplot_html, snp_distribution_html=snp_distribution_html,  homoplastic_snps_html=homoplastic_snps_html, \
                heatmap_html=heatmap_html, majority_percentage=majority_percentage, majority_threshold=majority_threshold, metadata_json_string=metadata_json_string)
    return html_template


def interactive_threshold_heatmap(service_config, metadata_json, majority_threshold):
    with open(service_config) as file:
        data = json.load(file)
    work_dir = data["work_data_dir"]

    # Paths for both data types for all three SNP subsets
    file_paths = {
        "all":      {"report_path": os.path.join(work_dir, "all_kSNPdist.report"),
                     "matrix_path": os.path.join(work_dir, "all_kSNPdist.matrix")},
        "core":     {"report_path": os.path.join(work_dir, "core_kSNPdist.report"),
                     "matrix_path": os.path.join(work_dir, "core_kSNPdist.matrix")},
        "majority": {"report_path": os.path.join(work_dir, "majority_kSNPdist.report"),
                     "matrix_path": os.path.join(work_dir, "majority_kSNPdist.matrix")},
    }

    any_file_found = any(
        os.path.exists(p)
        for subset in file_paths.values()
        for p in subset.values()
    )
    if not any_file_found:
        msg = "kSNP4 distance matrices are not available; cannot create heatmap. Please review the genome lengths for any outlier genomes, or re-run this job from the Jobs page.\n"
        sys.stderr.write(msg)
        heatmap_html = (
            '<div style="display: flex; justify-content: center; margin: 16px 0;">'
            '<div style="border: 2px solid #b00020; background-color: #fff0f0; color: #b00020; '
            'padding: 16px 20px; border-radius: 6px; max-width: 700px; text-align: center;">'
            '<strong>&#9888; kSNP4 distance matrices are not available.</strong> '
            'Please review the genome lengths for any outlier genomes, or re-run this job from the Jobs page.'
            '</div>'
            '</div>'
        )
        return heatmap_html, ""

    def load_subset(report_path, matrix_path):
        ids_rep, mat_rep, ids_mat, mat_mat = "null", "null", "null", "null"
        if os.path.exists(report_path):
            gids, mat = read_ksnp_distance_report(report_path)
            cl, cm = cluster_heatmap_data(gids, mat)
            ids_rep = json.dumps(cl)
            mat_rep = json.dumps(cm)
        if os.path.exists(matrix_path):
            gids, mat = read_ksnp_distance_matrix(matrix_path)
            cl, cm = cluster_heatmap_data(gids, mat)
            ids_mat = json.dumps(cl)
            mat_mat = json.dumps(cm)
        return ids_rep, mat_rep, ids_mat, mat_mat

    all_ids_rep,      all_mat_rep,      all_ids_mat,      all_mat_mat      = load_subset(**file_paths["all"])
    core_ids_rep,     core_mat_rep,     core_ids_mat,     core_mat_mat     = load_subset(**file_paths["core"])
    majority_ids_rep, majority_mat_rep, majority_ids_mat, majority_mat_mat = load_subset(**file_paths["majority"])
    # format the metadata into a string for the report
    metadata_json_string, metadata_df = create_metadata_table(metadata_json, "metadata.tsv")
    heatmap_template = """
    <!-- Plotly.js v3.0.1  — last updated June 2025 -->
    <script src="https://cdn.plot.ly/plotly-3.0.1.min.js"></script>
    <h3>SNP Distance Heatmap and Metadata</h3>
    <p>The kSNPdist tool measures genetic differences between genomes by comparing single nucleotide
    polymorphisms (SNPs), which are small variations in DNA sequence. It calculates how similar or different
    each pair of genomes is based only on SNP positions shared between both genomes, while ignoring missing
    data. The analysis produces two outputs: a matrix summarizing the genetic distances between all genomes
    in the dataset, and a report listing the number of SNP differences observed for every genome pair.
    Choose your view according to your goal:</p>
    <ul>
      <li><strong>Distance matrix (kSNPdist.matrix):</strong> Best for getting a broad overview of how
      closely related genomes are across an entire dataset. Because the results are organized in a matrix
      format, they are useful for clustering analyses, heatmaps, phylogenetic studies, and identifying
      groups of similar strains quickly.</li>
      <li><strong>Pairwise SNP report (kSNPdist.report):</strong> Best when you need the exact number of
      SNP differences between specific genome pairs. This format is easier for detailed comparisons,
      outbreak investigations, validation work, or reporting precise genomic distances between selected
      samples.</li>
    </ul>
    <p>Use the <strong>Data Source</strong> selector to toggle between views. The Heatmap, Close Genome
    Pairs, and Distance Matrix Table will all update to reflect the selected data source and SNP subset.</p>

    <!-- SNP subset + data source selectors — apply to all tabs -->
    <div style="display:flex; flex-wrap:wrap; gap:24px; align-items:center; margin-bottom:12px;">
      <label>Choose SNP Subset:
        <select id="matrixSelector" onchange="onMatrixChange()">
          <option value="1">All SNPs</option>
          <option value="2">Core SNPs</option>
          <option value="3">Majority SNPs</option>
        </select>
      </label>
      <label>Data Source:
        <select id="dataSourceSelector" onchange="onMatrixChange()">
          <option value="report">Pairwise SNP report (kSNPdist.report)</option>
          <option value="matrix">Distance matrix (kSNPdist.matrix)</option>
        </select>
      </label>
    </div>

    <!-- Single distance threshold — drives Close Genome Pairs and Distance Matrix Table -->
    <div style="display:flex; flex-wrap:wrap; gap:16px; align-items:center; margin-bottom:12px;">
      <label style="font-weight:bold;">Distance threshold (Close Pairs / Distance Table):
        <input type="number" id="linkageThreshold" min="0" style="width:70px; margin-left:6px;">
      </label>
    </div>

    <!-- Tab buttons -->
    <div style="display:flex; gap:0; margin-bottom:16px; border-bottom:2px solid #ccc;">
      <button id="tabHeatmap" onclick="switchView('heatmap')"
        style="padding:8px 20px; cursor:pointer; border:1px solid #ccc; border-bottom:none;
               background:#fff; font-size:14px; font-weight:bold; border-radius:4px 4px 0 0;
               margin-bottom:-2px; border-bottom:2px solid #fff;">
        Heatmap
      </button>
      <button id="tabClosePairs" onclick="switchView('closePairs')"
        style="padding:8px 20px; cursor:pointer; border:1px solid #ccc; border-bottom:none;
               background:#f5f5f5; font-size:14px; color:#555; border-radius:4px 4px 0 0;
               margin-bottom:-2px;">
        Close Genome Pairs
      </button>
      <button id="tabDistTable" onclick="switchView('distTable')"
        style="padding:8px 20px; cursor:pointer; border:1px solid #ccc; border-bottom:none;
               background:#f5f5f5; font-size:14px; color:#555; border-radius:4px 4px 0 0;
               margin-bottom:-2px;">
        Distance Matrix Table
      </button>
    </div>

    <!-- ====================== HEATMAP VIEW ====================== -->
    <div id="heatmapViewSection">
      <div class="heatmap-controls">
        <h4>Filter and Sort the Data:</h4>
        <label>Reorder Heatmap:
          <select id="metadataFieldSelect" onchange="recolorHeatmap()">
            <!-- options populated dynamically -->
          </select>
        </label>
      </div>
      <div class="linkage-controls" style="display:flex; flex-wrap:wrap; gap:10px; align-items:center;">
        <h4>Recolor Heatmap According to Linkage Thresholds:</h4>
        <label>Strong Linkage Thresholds:
          <input type="number" id="t0" value="0" disabled style="width:40px;">
          <input type="number" id="t1a" value="10" style="width:40px;">
        </label>
        <label>Mid Linkage Thresholds:
          <input type="number" id="t1b" value="10" style="width:40px;">
          <input type="number" id="t2a" value="40" style="width:40px;">
        </label>
        <label>Weak Linkage Thresholds:
          <input type="number" id="t2b" value="40" style="width:40px;">
          <input type="number" id="t3" placeholder="Max" disabled style="width:40px;">
        </label>
        <button style="padding:8px 8px; font-size:14px;" onclick="applyThresholdColoring()">Recolor</button>
      </div>
      <div id="heatmap" style="width:100%;"></div>
      <!-- Genome comparison panel (appears on cell click) -->
      <div id="comparisonPanel" style="display:none; margin-top:18px; border:1px solid #ccc;
           border-radius:6px; padding:16px; background:#fafafa;">
        <h3 id="comparisonTitle" style="margin-top:0;"></h3>
        <div style="display:flex; gap:24px; flex-wrap:wrap;">
          <div style="flex:1; min-width:220px;">
            <h4 id="genome1Label" style="margin-bottom:6px;"></h4>
            <table id="meta1Table" style="width:100%; border-collapse:collapse;"></table>
          </div>
          <div style="flex:1; min-width:220px;">
            <h4 id="genome2Label" style="margin-bottom:6px;"></h4>
            <table id="meta2Table" style="width:100%; border-collapse:collapse;"></table>
          </div>
        </div>
      </div>
    </div><!-- end heatmapViewSection -->

    <!-- ====================== CLOSE GENOME PAIRS ====================== -->
    <div id="closePairsViewSection" style="display:none;">
      <h3>Close Genome Pairs</h3>
      <p>This view lists each unique genome pair sorted by SNP distance, making it easy to identify
      the most closely related samples. Only the upper triangle of the distance matrix is shown (each
      pair appears once). Refine results by changing the distance threshold above. Click any column
      header to sort.</p>
      <div style="margin-bottom:12px;">
        <span id="cpCount" style="color:#555; font-size:13px;"></span>
      </div>
      <div style="overflow:auto; max-height:400px; border:1px solid #ccc; border-radius:4px;">
        <table id="cpTable" style="border-collapse:collapse; width:100%;">
          <thead>
            <tr style="position:sticky; top:0; background:#f0f0f0; z-index:2;">
              <th onclick="sortClosePairs('g1')"
                  style="padding:6px 12px; border:1px solid #ddd; cursor:pointer; user-select:none;">
                Genome A <span id="cpSort_g1"></span></th>
              <th onclick="sortClosePairs('g2')"
                  style="padding:6px 12px; border:1px solid #ddd; cursor:pointer; user-select:none;">
                Genome B <span id="cpSort_g2"></span></th>
              <th onclick="sortClosePairs('dist')"
                  style="padding:6px 12px; border:1px solid #ddd; cursor:pointer; user-select:none;">
                SNP Distance <span id="cpSort_dist">&#9650;</span></th>
            </tr>
          </thead>
          <tbody id="cpBody"></tbody>
        </table>
      </div>
    </div><!-- end closePairsViewSection -->

    <!-- ====================== DISTANCE MATRIX TABLE ====================== -->
    <div id="distanceTableViewSection" style="display:none;">
      <h3>Distance Matrix Table</h3>
      <p>The table below presents the full pairwise SNP distance matrix in a searchable, color-coded
      format. Each cell shows the SNP distance between two genomes. Use the <strong>search box</strong>
      to filter rows by genome ID, or set a distance threshold above to highlight cells at or below
      that value.</p>
      <div style="display:flex; flex-wrap:wrap; gap:12px; align-items:center; margin-bottom:10px;">
        <label style="display:flex; align-items:center; gap:6px;">
          Search genome ID:
          <input type="text" id="dmSearch" placeholder="Filter rows by ID..."
                 style="padding:4px 8px; width:200px; font-size:13px;"
                 oninput="buildDistTable()">
        </label>
        <span id="dmRowCount" style="color:#555; font-size:13px;"></span>
      </div>
      <div style="display:flex; gap:16px; align-items:center; margin-bottom:10px; font-size:12px; flex-wrap:wrap;">
        <strong>Color key:</strong>
        <span id="dmKeyIdentical" style="background:#2ca02c; color:white; padding:2px 8px; border-radius:3px;">Identical (0)</span>
        <span id="dmKeyMatching" style="background:#a8d5a2; color:#333; padding:2px 8px; border-radius:3px;">Matching threshold</span>
        <span id="dmKeyAbove"    style="background:#f4a460; color:#333; padding:2px 8px; border-radius:3px;">Above threshold</span>
      </div>
      <div id="dmTableWrapper"
           style="overflow:auto; max-height:600px; border:1px solid #ccc; border-radius:4px;">
        <table id="dmTable" style="border-collapse:collapse; white-space:nowrap;"></table>
      </div>
    </div><!-- end distanceTableViewSection -->

    <!-- ====================== PHYLOGENETIC TREES ====================== -->
    <h3 style="margin-top:32px;">Phylogenetic Trees</h3>
    <p>View a static image of the trees according to each SNP subset and each tree building method.
    An interactive view annotated with metadata is available in the job results in the directory
    labeled "Trees". Access the interactive tree viewer via files ending with ".tre" or ".phyloxml"
    in your job results directory.</p>
    <div style="display:flex; flex-wrap:wrap; gap:16px; align-items:center; margin-bottom:12px;">
      <label>Choose SNP Subset:
        <select id="treeSNPSelector" onchange="updateSVG()">
          <option value="1">All</option>
          <option value="2">Core</option>
          <option value="3">Majority</option>
        </select>
      </label>
      <label>Tree Building Method:
        <select id="methodSelector" onchange="updateSVG()">
          <option value="ML">Maximum Likelihood</option>
          <option value="NJ">Neighbor Joining</option>
          <option value="parsimony">Parsimony</option>
        </select>
      </label>
    </div>
    <div id="svgContainer"></div>

    <script>
        // ===== Embedded data =====
        // Pairwise SNP report (kSNPdist.report) — raw integer SNP counts
        const labelsAll_rep      = {all_ids_rep};
        const matrixAll_rep      = {all_mat_rep};
        const labelsCore_rep     = {core_ids_rep};
        const matrixCore_rep     = {core_mat_rep};
        const labelsMaj_rep      = {majority_ids_rep};
        const matrixMaj_rep      = {majority_mat_rep};
        // Distance matrix (kSNPdist.matrix) — proportional float distances
        const labelsAll_mat      = {all_ids_mat};
        const matrixAll_mat      = {all_mat_mat};
        const labelsCore_mat     = {core_ids_mat};
        const matrixCore_mat     = {core_mat_mat};
        const labelsMaj_mat      = {majority_ids_mat};
        const matrixMaj_mat      = {majority_mat_mat};

        const metadata      = {metadata_json_string};

        const idToMeta = {{}};
        metadata.forEach(obj => {{ idToMeta[obj.genome_id] = obj; }});

        let currentLabels = [];
        let currentMatrix = [];
        // 'viridis' on initial load; 'threshold' after Recolor is clicked
        let heatmapColorMode = 'viridis';

        // ===== Return the active matrix based on SNP subset + data source selectors =====
        function getActiveMatrix() {{
            const subset = document.getElementById('matrixSelector').value;
            const src    = document.getElementById('dataSourceSelector').value;
            const sets = {{
                "1": {{ report: {{ labels: labelsAll_rep,  matrix: matrixAll_rep  }},
                        matrix: {{ labels: labelsAll_mat,  matrix: matrixAll_mat  }} }},
                "2": {{ report: {{ labels: labelsCore_rep, matrix: matrixCore_rep }},
                        matrix: {{ labels: labelsCore_mat, matrix: matrixCore_mat }} }},
                "3": {{ report: {{ labels: labelsMaj_rep,  matrix: matrixMaj_rep  }},
                        matrix: {{ labels: labelsMaj_mat,  matrix: matrixMaj_mat  }} }},
            }};
            const chosen = sets[subset][src];
            if (!chosen || !chosen.labels) {{
                return {{ labels: [], matrix: [] }};
            }}
            return {{ labels: chosen.labels.slice(), matrix: chosen.matrix.map(r => r.slice()) }};
        }}

        // ===== SNP subset change — refresh all visible views =====
        function onMatrixChange() {{
            recolorHeatmap();
            updateSVG();
            const cpSec = document.getElementById('closePairsViewSection');
            if (cpSec && cpSec.style.display !== 'none') buildClosePairs();
            const dtSec = document.getElementById('distanceTableViewSection');
            if (dtSec && dtSec.style.display !== 'none') buildDistTable();
        }}

        // ===== Sync paired threshold inputs (validate only; no auto-recolor) =====
        function syncThresholdInputs() {{
            const t1a = document.getElementById('t1a');
            const t1b = document.getElementById('t1b');
            const t2a = document.getElementById('t2a');
            const t2b = document.getElementById('t2b');

            function validateThresholds() {{
                const t1 = parseFloat(t1a.value);
                const t2 = parseFloat(t2a.value);
                if (t1 > t2) {{
                    alert("Strong threshold must be ≤ Mid threshold. Resetting.");
                    t1a.value = t2;
                    t1b.value = t2;
                }}
            }}

            t1a.addEventListener('input', () => {{ t1b.value = t1a.value; validateThresholds(); }});
            t1b.addEventListener('input', () => {{ t1a.value = t1b.value; validateThresholds(); }});
            t2a.addEventListener('input', () => {{ t2b.value = t2a.value; validateThresholds(); }});
            t2b.addEventListener('input', () => {{ t2a.value = t2b.value; validateThresholds(); }});
        }}

        document.addEventListener('DOMContentLoaded', function () {{
            syncThresholdInputs();
            const ltEl = document.getElementById('linkageThreshold');
            if (ltEl) ltEl.addEventListener('input', function () {{
                buildClosePairs();
                buildDistTable();
            }});
        }});

        // ===== Populate metadata reorder dropdown =====
        (function populateMetadataFields() {{
            const select = document.getElementById('metadataFieldSelect');
            const allKeys = Object.keys(metadata[0]).filter(k => k !== "id");
            const defaultOption = document.createElement('option');
            defaultOption.value = "";
            defaultOption.textContent = "Hierarchical Clustering";
            select.appendChild(defaultOption);
            allKeys.forEach(field => {{
                const opt = document.createElement('option');
                opt.value = field;
                opt.textContent = field;
                select.appendChild(opt);
            }});
        }})();

        // ===== Binning helpers for threshold coloring =====
        function assignBin(val, t1, t2, maxVal) {{
            if (val === 0)     return 0;
            if (val < t1)     return 1;
            if (val < t2)     return 2;
            if (val < maxVal) return 3;
            return 4;
        }}

        function getColorScale() {{
            return [
                [0.0,  '#440154'],
                [0.25, '#3b528b'],
                [0.5,  '#21918c'],
                [0.75, '#5ec962'],
                [1.0,  '#fde725']
            ];
        }}

        // ===== Reorder matrix by metadata field =====
        function reorderByField(fieldName, labelsArr, matrixArr) {{
            const arr = metadata.map(obj => ({{ id: obj.genome_id, value: obj[fieldName] }}));
            arr.sort((a, b) => {{
                if (a.value < b.value) return -1;
                if (a.value > b.value) return 1;
                if (a.id < b.id) return -1;
                if (a.id > b.id) return 1;
                return 0;
            }});
            const newLabels = arr.map(o => o.id);
            const indexMap = {{}};
            labelsArr.forEach((lbl, idx) => {{ indexMap[lbl] = idx; }});
            const n = labelsArr.length;
            const newMatrix = [];
            for (let i = 0; i < n; i++) {{
                const origRow = indexMap[newLabels[i]];
                const row = [];
                for (let j = 0; j < n; j++) row.push(matrixArr[origRow][indexMap[newLabels[j]]]);
                newMatrix.push(row);
            }}
            return {{ newLabels, newMatrix }};
        }}

        // ===== Render metadata key/value table =====
        function renderMetaTable(tableEl, metaObj) {{
            tableEl.innerHTML = '';
            for (const [field, val] of Object.entries(metaObj)) {{
                const tr  = document.createElement('tr');
                const tdK = document.createElement('td');
                const tdV = document.createElement('td');
                tdK.style.cssText = 'padding:4px 8px 4px 0; font-weight:bold; white-space:nowrap; vertical-align:top;';
                tdV.style.cssText = 'padding:4px 0; vertical-align:top;';
                tdK.textContent = field;
                tdV.textContent = val;
                tr.appendChild(tdK);
                tr.appendChild(tdV);
                tableEl.appendChild(tr);
            }}
        }}

        // ===== Heatmap cell click → comparison panel =====
        function onHeatmapClick(eventData) {{
            if (!eventData || !eventData.points || eventData.points.length === 0) return;
            const pt   = eventData.points[0];
            const id1  = pt.y;
            const id2  = pt.x;
            const i1   = currentLabels.indexOf(id1);
            const i2   = currentLabels.indexOf(id2);
            const dist = (i1 >= 0 && i2 >= 0) ? currentMatrix[i1][i2] : '';
            const meta1 = idToMeta[id1] || {{ genome_id: id1 }};
            const meta2 = idToMeta[id2] || {{ genome_id: id2 }};
            document.getElementById('comparisonTitle').textContent = `SNP Distance: ${{dist}}`;
            document.getElementById('genome1Label').innerHTML = `<a href="https://www.bv-brc.org/view/Genome/${{id1}}" target="_blank">${{id1}}</a>`;
            document.getElementById('genome2Label').innerHTML = `<a href="https://www.bv-brc.org/view/Genome/${{id2}}" target="_blank">${{id2}}</a>`;
            renderMetaTable(document.getElementById('meta1Table'), meta1);
            renderMetaTable(document.getElementById('meta2Table'), meta2);
            const panel = document.getElementById('comparisonPanel');
            panel.style.display = 'block';
            panel.scrollIntoView({{ behavior: 'smooth', block: 'nearest' }});
        }}

        // ===== Draw / redraw heatmap =====
        function recolorHeatmap() {{
            const t1 = parseInt(document.getElementById('t1a').value);
            const t2 = parseInt(document.getElementById('t2a').value);
            const metaField = document.getElementById('metadataFieldSelect').value;

            let {{ labels: genomeLabels, matrix: snpMatrix }} = getActiveMatrix();

            if (genomeLabels.length === 0) {{
                document.getElementById('heatmap').innerHTML =
                    '<p style="color:#666; font-style:italic; padding:16px;">No distance data available for this SNP subset and data source combination.</p>';
                return;
            }}

            if (metaField) {{
                const {{ newLabels, newMatrix }} = reorderByField(metaField, genomeLabels, snpMatrix);
                genomeLabels = newLabels;
                snpMatrix = newMatrix;
            }}

            currentLabels = genomeLabels;
            currentMatrix = snpMatrix;

            const hoverText = snpMatrix.map((row, i) =>
                row.map((val, j) => `${{genomeLabels[i]}} vs ${{genomeLabels[j]}}<br>SNP Distance: ${{val}}<br><i>Click to show metadata</i>`)
            );

            let zData, colorscale, colorbarConfig, extraRange;
            if (heatmapColorMode === 'threshold') {{
                const maxVal = Math.max(...snpMatrix.flat());
                zData        = snpMatrix.map(row => row.map(val => assignBin(val, t1, t2, maxVal)));
                colorscale   = getColorScale();
                extraRange   = {{ zmin: 0, zmax: 4 }};
                colorbarConfig = {{
                    tickvals: [0, 1, 2, 3, 4],
                    ticktext: ['Zero', 'Strong', 'Mid', 'Weak', 'Max Value'],
                    title: 'Linkage Strength'
                }};
            }} else {{
                zData        = snpMatrix;
                colorscale   = 'Viridis';
                extraRange   = {{}};
                colorbarConfig = {{ title: 'SNP Distance' }};
            }}

            const titleStr = heatmapColorMode === 'threshold'
                ? `SNP Heatmap (Thresholds: ${{t1}}, ${{t2}})` + (metaField ? ` – Reordered by "${{metaField}}"` : '')
                : 'SNP Distance Heatmap'                        + (metaField ? ` – Reordered by "${{metaField}}"` : '');

            const traceData = [Object.assign({{
                z: zData, x: genomeLabels, y: genomeLabels,
                type: 'heatmap', colorscale: colorscale,
                text: hoverText, hoverinfo: 'text', colorbar: colorbarConfig
            }}, extraRange)];

            const layout = {{
                title: titleStr,
                xaxis: {{ type: 'category', tickangle: 45 }},
                yaxis: {{ type: 'category', tickangle: 45 }}
            }};

            Plotly.newPlot('heatmap', traceData, layout);
            document.getElementById('heatmap').on('plotly_click', onHeatmapClick);
        }}

        // ===== Switch to threshold coloring (Recolor button) =====
        function applyThresholdColoring() {{
            heatmapColorMode = 'threshold';
            recolorHeatmap();
        }}

        // ===== Tab switching =====
        function switchView(view) {{
            const sections = {{
                heatmap:    document.getElementById('heatmapViewSection'),
                closePairs: document.getElementById('closePairsViewSection'),
                distTable:  document.getElementById('distanceTableViewSection'),
            }};
            const tabs = {{
                heatmap:    document.getElementById('tabHeatmap'),
                closePairs: document.getElementById('tabClosePairs'),
                distTable:  document.getElementById('tabDistTable'),
            }};
            Object.keys(sections).forEach(k => {{ if (sections[k]) sections[k].style.display = 'none'; }});
            Object.keys(tabs).forEach(k => {{
                if (tabs[k]) {{
                    tabs[k].style.background   = '#f5f5f5';
                    tabs[k].style.color        = '#555';
                    tabs[k].style.borderBottom = '';
                    tabs[k].style.fontWeight   = '';
                }}
            }});
            if (sections[view]) sections[view].style.display = '';
            if (tabs[view]) {{
                tabs[view].style.background   = '#fff';
                tabs[view].style.color        = '';
                tabs[view].style.borderBottom = '2px solid #fff';
                tabs[view].style.fontWeight   = 'bold';
            }}
            if (view === 'closePairs') buildClosePairs();
            if (view === 'distTable')  buildDistTable();
        }}

        // ===== Shared color helper for Close Pairs and Distance Matrix =====
        function getDmColor(val, t) {{
            if (val === 0) return {{ bg: '#2ca02c', fg: 'white' }};
            if (t !== null && val <= t) return {{ bg: '#a8d5a2', fg: '#333' }};
            return {{ bg: t !== null ? '#f4a460' : '#ffffff', fg: '#333' }};
        }}

        function getDmThreshold() {{
            const el  = document.getElementById('linkageThreshold');
            const raw = el ? el.value.trim() : '';
            return raw !== '' ? parseInt(raw) : null;
        }}

        // ===== Close Genome Pairs =====
        let _cpPairs     = [];
        let _cpSortField = 'dist';
        let _cpSortAsc   = true;

        function buildClosePairs() {{
            const t = getDmThreshold();
            const {{ labels, matrix }} = getActiveMatrix();
            const n = labels.length;
            _cpPairs = [];
            for (let i = 0; i < n; i++) {{
                for (let j = i + 1; j < n; j++) {{
                    const d = matrix[i][j];
                    if (t === null || d <= t) _cpPairs.push({{ g1: labels[i], g2: labels[j], dist: d }});
                }}
            }}
            _renderClosePairs(t);
        }}

        function sortClosePairs(field) {{
            if (_cpSortField === field) {{ _cpSortAsc = !_cpSortAsc; }}
            else {{ _cpSortField = field; _cpSortAsc = (field === 'dist'); }}
            _renderClosePairs(getDmThreshold());
        }}

        function _renderClosePairs(t) {{
            ['g1', 'g2', 'dist'].forEach(f => {{
                const el = document.getElementById('cpSort_' + f);
                if (el) el.textContent = '';
            }});
            const sortEl = document.getElementById('cpSort_' + _cpSortField);
            if (sortEl) sortEl.textContent = _cpSortAsc ? ' ▲' : ' ▼';

            const sorted = _cpPairs.slice().sort((a, b) => {{
                const av = a[_cpSortField], bv = b[_cpSortField];
                const cmp = av < bv ? -1 : av > bv ? 1 : 0;
                return _cpSortAsc ? cmp : -cmp;
            }});

            const tbody = document.getElementById('cpBody');
            tbody.innerHTML = '';
            if (sorted.length === 0) {{
                const tr = document.createElement('tr');
                const td = document.createElement('td');
                td.colSpan = 3;
                td.style.cssText = 'padding:14px; text-align:center; color:#666; font-style:italic;';
                td.textContent = 'No pairs found within the specified threshold. Try increasing the threshold.';
                tr.appendChild(td); tbody.appendChild(tr);
            }} else {{
                sorted.forEach(pair => {{
                    const {{ bg, fg }} = getDmColor(pair.dist, t);
                    const tr   = document.createElement('tr');
                    const tdg1 = document.createElement('td');
                    const tdg2 = document.createElement('td');
                    const tdd  = document.createElement('td');
                    tdg1.style.cssText = 'padding:4px 12px; border:1px solid #ddd; font-size:12px;';
                    tdg2.style.cssText = 'padding:4px 12px; border:1px solid #ddd; font-size:12px;';
                    tdd.style.cssText  = 'padding:4px 12px; border:1px solid #ddd; font-size:12px; text-align:center; font-weight:bold; background:' + bg + '; color:' + fg + ';';
                    tdg1.textContent = pair.g1;
                    tdg2.textContent = pair.g2;
                    tdd.textContent  = pair.dist;
                    tr.appendChild(tdg1); tr.appendChild(tdg2); tr.appendChild(tdd);
                    tbody.appendChild(tr);
                }});
            }}
            document.getElementById('cpCount').textContent =
                sorted.length === 1 ? '1 pair found' : sorted.length + ' pairs found';
        }}

        // ===== Full Distance Matrix Table =====
        function buildDistTable() {{
            const searchVal = document.getElementById('dmSearch') ? document.getElementById('dmSearch').value.trim().toLowerCase() : '';
            const t = getDmThreshold();
            const {{ labels, matrix }} = getActiveMatrix();
            const n = labels.length;

            const kmEl = document.getElementById('dmKeyMatching');
            const kaEl = document.getElementById('dmKeyAbove');
            if (kmEl) kmEl.textContent = t !== null ? `Matching threshold (≤${{t}})` : 'Matching threshold';
            if (kaEl) kaEl.textContent = t !== null ? `Above threshold (>${{t}})`          : 'Above threshold';

            const visibleRows = [];
            labels.forEach((lbl, i) => {{
                if (!searchVal || lbl.toLowerCase().includes(searchVal)) visibleRows.push(i);
            }});

            const table    = document.getElementById('dmTable');
            const fragment = document.createDocumentFragment();

            const thead  = document.createElement('thead');
            const hRow   = document.createElement('tr');
            const corner = document.createElement('th');
            corner.style.cssText = 'position:sticky; top:0; left:0; z-index:4; background:#fff; padding:4px 8px; border:1px solid #ddd; min-width:120px;';
            corner.textContent = '';
            hRow.appendChild(corner);
            labels.forEach(lbl => {{
                const th = document.createElement('th');
                th.style.cssText = 'position:sticky; top:0; z-index:2; background:#f8f8f8; border:1px solid #ddd; padding:2px; font-size:10px; font-weight:normal; writing-mode:vertical-rl; transform:rotate(180deg); height:110px; vertical-align:bottom; text-align:left;';
                th.textContent = lbl;
                hRow.appendChild(th);
            }});
            thead.appendChild(hRow);
            fragment.appendChild(thead);

            const tbody = document.createElement('tbody');
            visibleRows.forEach(i => {{
                const tr    = document.createElement('tr');
                const rowTh = document.createElement('th');
                rowTh.style.cssText = 'position:sticky; left:0; z-index:1; background:#f8f8f8; border:1px solid #ddd; padding:3px 8px; font-size:11px; font-weight:normal; white-space:nowrap; text-align:left;';
                rowTh.textContent = labels[i];
                tr.appendChild(rowTh);
                labels.forEach((lbl, j) => {{
                    const val        = matrix[i][j];
                    const {{ bg, fg }} = getDmColor(val, t);
                    const td         = document.createElement('td');
                    td.style.cssText = 'background:' + bg + '; color:' + fg + '; padding:3px 5px; border:1px solid #ddd; text-align:center; font-size:11px; min-width:28px;';
                    td.textContent   = val;
                    td.title         = labels[i] + ' vs ' + lbl + ': ' + val + ' SNP differences';
                    tr.appendChild(td);
                }});
                tbody.appendChild(tr);
            }});
            if (visibleRows.length === 0) {{
                const emptyTr = document.createElement('tr');
                const emptyTd = document.createElement('td');
                emptyTd.colSpan = n + 1;
                emptyTd.style.cssText = 'padding:16px; text-align:center; color:#666; font-style:italic;';
                emptyTd.textContent   = 'No genomes match the current filters.';
                emptyTr.appendChild(emptyTd);
                tbody.appendChild(emptyTr);
            }}
            fragment.appendChild(tbody);
            table.innerHTML = '';
            table.appendChild(fragment);

            document.getElementById('dmRowCount').textContent =
                visibleRows.length === n
                    ? 'Showing all ' + n + ' genomes'
                    : 'Showing ' + visibleRows.length + ' of ' + n + ' genomes';
        }}

        // ===== Tree viewer =====
        function updateSVG() {{
            const data_input = document.getElementById('treeSNPSelector').value;
            const method = document.getElementById('methodSelector').value;
            const fp_map = {{
                "1": "SNPs_all",
                "2": "core_SNPs",
                "3": "SNPs_in_majority{majority_threshold}"
            }};
            const tree_type = fp_map[data_input] || "SNPs_all";
            const svgPath = `report_supporting_documents/tree.${{tree_type}}.${{method}}.tre.svg`;
            fetch(svgPath)
                .then(response => {{
                    if (!response.ok) throw new Error('SVG not found');
                    return response.text();
                }})
                .then(svgContent => {{
                    document.getElementById('svgContainer').innerHTML = svgContent;
                }})
                .catch(() => {{
                    document.getElementById('svgContainer').innerHTML =
                        '<p style="color:#555;">Tree not available for this combination. ' +
                        'Ensure the <code>report_supporting_documents</code> directory is in ' +
                        'the same folder as this report file.</p>';
                }});
        }}

        // Initial render — Viridis by default
        recolorHeatmap();
        updateSVG();
    </script>
    """.format(
    all_ids_rep=all_ids_rep,
    all_mat_rep=all_mat_rep,
    core_ids_rep=core_ids_rep,
    core_mat_rep=core_mat_rep,
    majority_ids_rep=majority_ids_rep,
    majority_mat_rep=majority_mat_rep,
    all_ids_mat=all_ids_mat,
    all_mat_mat=all_mat_mat,
    core_ids_mat=core_ids_mat,
    core_mat_mat=core_mat_mat,
    majority_ids_mat=majority_ids_mat,
    majority_mat_mat=majority_mat_mat,
    metadata_json_string=metadata_json_string,
    majority_threshold=majority_threshold,
    )
    return heatmap_template, metadata_json_string


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


def generate_table_html_2(kchooser_df, table_width='75%'):
    # Generate table headers
    headers = ''.join(f'<th>{header}</th>' for header in kchooser_df.columns)
    rows = ''

    # Generate table rows
    for _, row in kchooser_df.iterrows():
        row_html = ''
        for column in kchooser_df.columns:
            cell_value = row[column]
            # replacing is numeric function with pd.api.types.is_numeric_dtype
            if pd.api.types.is_numeric_dtype(type(cell_value)):
                # Apply number formatting
                if isinstance(cell_value, (int, np.integer)):
                    formatted_value = f"{cell_value:,}"  # Comma formatting for integers
                elif isinstance(cell_value, (float, np.floating)):
                    formatted_value = f"{cell_value:,.2f}"  # Two decimals for floats
                else:
                    formatted_value = str(cell_value)
                row_html += f'<td style="text-align: center;">{formatted_value}</td>'
            else:
                row_html += f'<td>{cell_value}</td>'
        rows += f'<tr>{row_html}</tr>'

    # Construct the complete HTML table with specified width
    table_html = f'''
    <div style="text-align: center;">
    <table style="margin: auto; width: {table_width}; border-collapse: collapse; border: 1px solid black;">
        <thead>
            <tr>{headers}</tr>
        </thead>
        <tbody>
            {rows}
        </tbody>
    </table>
    </div>
    '''
    return table_html


def infer_output_subtype(filename):
    if "core_SNPs" in filename:
        return "Core_SNPs"
    elif "majority" in filename:
        return "Majority_SNPs"
    elif "SNPs_all" in filename:
        return "All_SNPs"


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


def make_genome_bar_chart(data, report_data, majority_threshold):
    if "COUNT_coreSNPs" not in report_data or "COUNT_SNPs" not in report_data:
        msg = "SNP count files not found; skipping SNP distribution chart.\n"
        sys.stderr.write(msg)
        return "<p>SNP distribution chart not available: SNP count data is missing.</p>"
    try:
        majority_snps_value = report_data["COUNT_coreSNPs"][0]["Number SNPs in at least a fraction {} of genomes".format(majority_threshold)]
        core_snps_value = report_data["COUNT_coreSNPs"][0]["Number core SNPs"]
        total_snps_value = report_data["COUNT_SNPs"][0]["Number_SNPs"]
    except (KeyError, IndexError) as e:
        msg = "SNP count data is incomplete ({}); skipping SNP distribution chart.\n".format(e)
        sys.stderr.write(msg)
        return "<p>SNP distribution chart not available: SNP count data is incomplete.</p>"
    
    categories = ["Total SNPs", "Majority SNPs", "Core_SNPs"] 
    values = [total_snps_value, majority_snps_value, core_snps_value]
    
    # Plot in greatest to least
    # Zip and sort by value
    sorted_pairs = sorted(zip(values, categories), reverse=False)
    sorted_values, sorted_categories = zip(*sorted_pairs)
    
    #  Create figure with both stacked bar rows
    fig = go.Figure(go.Bar(
        x=sorted_values,
        y=sorted_categories,
        orientation="h"
    ))

    fig.update_layout(
    title='SNP Distribution by SNP subset',
    xaxis_title='Count',
    yaxis_title='SNP subset',
    template='plotly_white'
    )
    fig.write_html("snp_distribution_bar_chart_offline.html", include_plotlyjs=False)
    snp_distribution = read_plotly_html("snp_distribution_bar_chart_offline.html")
    return snp_distribution


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
        if firstword == "All" or filename.startswith("SNPs_all") or filename == "all_snp_distance_heatmap.html":
            All_SNPs_dir = os.path.join(destination_dir, "All_SNPs")
            os.makedirs(All_SNPs_dir, exist_ok=True)
            # Add extensions to SNPs_all and SNPs_all_matrix for the website
            if filename == "SNPs_all":
                new_path = os.path.join(All_SNPs_dir, (filename + ".tsv"))
                shutil.copy(file_path, new_path)
            elif filename == "SNPs_all_matrix":
                new_path = os.path.join(All_SNPs_dir, (filename + ".txt"))
                shutil.copy(file_path, new_path)
            # else copy all files to the all SNPs dir
            else:
                shutil.copy(file_path, All_SNPs_dir)
        if firstword == "annotate" and os.path.getsize(file_path) > 0:
            shutil.copy(file_path, intermediate_dir)
        if firstword == "ClusterInfo.SNPs" or firstword == "ClusterInfo.core":
            group = infer_output_subtype(filename)
            cluster_dir = os.path.join(destination_dir, group, "Cluster_Information")
            os.makedirs(cluster_dir, exist_ok=True)
            shutil.copy(file_path, cluster_dir)
        if firstword == "core" or firstword == "nonCore" or filename == "core_snp_distance_heatmap.html":
            core_snp_dir = os.path.join(destination_dir, "Core_SNPs")
            os.makedirs(core_snp_dir, exist_ok=True)
            # Add extensions to core_SNPs and core_SNPs_matrix for the website
            if filename == "core_SNPs":
                new_path = os.path.join(core_snp_dir, (filename + ".tsv"))
                shutil.copy(file_path, new_path)
            elif filename == "core_SNPs_matrix":
                new_path = os.path.join(core_snp_dir, (filename + ".txt"))
                shutil.copy(file_path, new_path)
            elif filename.startswith("core_kSNPdist"):
                # shutil.copy(file_path, work_dir)
                pass
            # else copy all files to the core SNPs dir
            else:
                shutil.copy(file_path, core_snp_dir)
        if firstword == "COUNT" or firstword == "tip" or firstword == "Node" or firstword == "NJ.dist.matrix":
            intermediate_dir = os.path.join(destination_dir, "Intermediate_Files")
            os.makedirs(intermediate_dir, exist_ok=True)
            shutil.copy(file_path, intermediate_dir)
        if firstword == "Homoplasy":
            group = infer_output_subtype(filename)
            homoplasy_dir = os.path.join(destination_dir, group, "Homoplasy")
            os.makedirs(homoplasy_dir, exist_ok=True)
            shutil.copy(file_path, homoplasy_dir)
        if "SNPs_in_majority" in filename and "matrix" in filename:
            group = infer_output_subtype(filename)
            majority_dir = os.path.join(destination_dir, group)
            os.makedirs(majority_dir, exist_ok=True)
            if len(filename) == 26:
                new_path = os.path.join(majority_dir, (filename + ".txt"))
                shutil.copy(file_path, new_path)
            else:
                shutil.copy(file_path, majority_dir)
        # Capture specifically SNPs_in_majority0.5 where 5 could be any integer 0-9
        if filename.startswith("SNPs_in_majority0.") or filename == "majority_snp_distance_heatmap.html":
            group = infer_output_subtype(filename)
            majority_dir = os.path.join(destination_dir, group)
            os.makedirs(majority_dir, exist_ok=True)
            # shutil.copy(file_path, majority_dir)
            if filename.startswith("SNPs_in_majority") and len(filename) == 19:
                new_path = os.path.join(majority_dir, (filename + ".tsv"))
                shutil.copy(file_path, new_path)
            # elif "SNPs_in_majority" in filename and len(filename) == 26:
            #     new_path = os.path.join(majority_dir, (filename + ".txt"))
            #     shutil.copy(file_path, new_path)
        if firstword.startswith("VCF"):
            VCFs_dir = os.path.join(destination_dir, "VCFs")
            os.makedirs(VCFs_dir, exist_ok=True)
            shutil.copy(file_path, VCFs_dir)
    # Trees get their own loop
    clean_tree_dir = os.path.join(work_dir,"clean_trees")
    for filename in os.listdir(clean_tree_dir):
        file_path = os.path.join(clean_tree_dir, filename)
        group = infer_output_subtype(filename)
        if filename == "tree_AlleleCounts.parsimony.tre" or filename == "tree_AlleleCounts.parsimony.tre.phyloxml":
            all_snps_dir = os.path.join(destination_dir, "All_SNPs", "Trees")
            shutil.copy(file_path, all_snps_dir)
        else:
            if group == "All_SNPs" or group == "Core_SNPs" or group == "Majority_SNPs":
                tree_dir = os.path.join(destination_dir, group, "Trees")
                os.makedirs(tree_dir, exist_ok=True)
                if os.path.splitext(file_path)[-1].lower() == ".tre":
                    newick_dir = os.path.join(tree_dir, "Newick_Files")
                    os.makedirs(newick_dir, exist_ok=True)
                    shutil.copy(file_path, newick_dir)
                else:
                    shutil.copy(file_path, tree_dir)
            else:
                tree_dir = os.path.join(destination_dir, "All_SNPs", "Trees")
                if os.path.splitext(file_path)[-1].lower() == ".tre":
                    newick_dir = os.path.join(tree_dir, "Newick_Files")
                    os.makedirs(newick_dir, exist_ok=True)
                    shutil.copy(file_path, newick_dir)
                else:
                    shutil.copy(file_path, tree_dir)


def parse_core_snps(file_path, filename):
    if os.path.getsize(file_path) > 0:
        with open(file_path, "r") as file:
            content = file.read()
            # Search via regex patterns
            core_snp_match = re.search(r'Number core SNPs:\s*(\d+)', content)
            non_core_snp_match = re.search(r'Number non-core SNPs:\s*(\d+)', content)
            fraction_snp_match = re.search(r'Number SNPs in at least a fraction ([\d.]+) of genomes:\s*(\d+)', content)

            count_data = {
                "core_SNPs": int(core_snp_match.group(1)) if core_snp_match else None,
                "non_core_SNPs": int(non_core_snp_match.group(1)) if non_core_snp_match else None,
                "fraction": float(fraction_snp_match.group(1)) if fraction_snp_match else None,
                "genome_count": int(fraction_snp_match.group(2)) if fraction_snp_match else None
            }
            add_to_report_dict(filename, count_data)

def parse_intermediate_files(report_data, work_dir):
    for filename in os.listdir(work_dir):
        file_path = os.path.join(work_dir, filename)
        firstword = filename.split("_")[0]
        if firstword == "COUNT":
            cs_data = {}
            with open(file_path) as file:
                for line in file:
                    if line.strip():  # Ensure the line is not empty
                        key, value = line.split(": ", 1)
                        cs_data[key.strip()] = int(value.strip())
                        report_data = add_to_report_dict(report_data, filename, cs_data)
    return report_data


def parse_kchooser_report(report_data, kchooser_report):
    kchooser_data = {}
    with open(kchooser_report, "r") as file:
        text = file.read()
    # Extract number of genomes
    match = re.search(r'There were (\d+) genomes', text)
    if match:
        kchooser_data["Total Genomes"] = int(match.group(1))

    # Extract median genome and its length
    match = re.search(r'The median length genome was (\S+)', text)
    if match:
        kchooser_data["Median Genome"] = match.group(1).replace("_", ".")

    match = re.search(r'Its length is (\d+)', text)
    if match:
        kchooser_data["Median Genome Length"] = int(match.group(1))
    # Extract shortest genome and its length
    match = re.search(r'The shortest genomes is (\S+) its length is (\d+)', text)
    if match:
        kchooser_data["Shortest Genome"] = match.group(1).replace("_", ".")
        kchooser_data["Shortest Genome Length"] = int(match.group(2))
    # return kchooser_data
    report_data = add_to_report_dict(report_data, "kchooser_report", kchooser_data)
    return report_data


def parse_node_file(report_data, node_file, filename):
    with open(node_file, "r") as file:
        content = file.read()   
        node_pattern = re.search(r'node:\s+(\S+)', content)
        targets_pattern = re.search(r'NumberTargets:\s+(\d+)', content)
        snps_pattern = re.search(r'NumberSNPs:\s+(\d+)', content)

        node_info = {
            "node_pattern": int(node_pattern.group(1)) if node_pattern else None,
            "targets_pattern": int(targets_pattern.group(1)) if targets_pattern else None,
            "snps_pattern": int(snps_pattern.group(1)) if snps_pattern else None,
            }
        report_data = add_to_report_dict(filename, node_info)
    return report_data


def parse_optimum_k(kchooser_report):
    with open(kchooser_report, 'r') as file:
        for line in file:
            match = re.search(r'The optimum value of k is (\d+)', line)
            if match:
                print(match.group(1))
                return
    print("Optimum value of k not found")
  

def read_ksnp_distance_matrix(ksnp_dist_matrix):
    df = pd.read_csv(ksnp_dist_matrix, sep='\t', header=0, index_col=None)
    genome_ids_raw = list(df.columns)
    df.index = genome_ids_raw
    # Round to 4 decimal places for readability
    matrix = [[round(v, 4) for v in row] for row in df.values.tolist()]
    genome_ids = [gid.replace("_", ".") for gid in genome_ids_raw]
    return genome_ids, matrix


def read_ksnp_distance_report(ksnp_dist_report):
    df = pd.read_csv(ksnp_dist_report, sep='\t', header=None)
    df.columns = ["value", "genome1", "genome2"]
    pivot_df = df.pivot(index="genome1", columns="genome2", values="value")
    # Ensure it's symmetric by filling missing cells
    pivot_df = pivot_df.combine_first(pivot_df.T)

    # Ensure genomes are in same order for rows/columns
    genomes = sorted(set(pivot_df.index) | set(pivot_df.columns))
    pivot_df = pivot_df.reindex(index=genomes, columns=genomes)

    # Fill diagonals or missing cells as needed
    pivot_df = pivot_df.fillna(0)    
    pivot_df.index = pivot_df.index.astype(str)
    pivot_df.columns = pivot_df.columns.astype(str)

    # give the data as lists of lists, replacing _ with . to match metadata IDs
    genome_ids = [gid.replace("_", ".") for gid in pivot_df.columns.tolist()]
    snpMatrix = pivot_df.values.tolist()
    return genome_ids, snpMatrix


def run_newick_to_phyloxml(clean_nwk):
        # Run phyloxml command
        result = subprocess.run(["p3x-newick-to-phyloxml", "--verbose", "-l", "genome_id", "-g", "collection_year,host_common_name,isolation_country,strain,genome_name,genome_id,accession,subtype,lineage,host_group,collection_date,geographic_group,geographic_location", clean_nwk])
        msg = "{}".format(result)
        sys.stderr.write(msg)


def read_plotly_html(plot_path):
    # Read the content from 'Variant_Plot_Interactive.html'
    with open(plot_path, 'r') as file:
        plotly_html_content = file.read()
    # Extract everything within the <body> tags
    extracted_content = re.findall(r'<body>(.*?)</body>', plotly_html_content, re.DOTALL)

    # Assuming extracted_content contains our needed Plotly graph initialization scripts
    plotly_graph_content = extracted_content[0] if extracted_content else ''
    return plotly_graph_content


def run_p3x_tree_to_svg(file_path, tree_svg_dir):
        # Run tree to svg command 
        subprocess.run(["p3x-tree-to-svg", file_path])
        tree_file_path = file_path + ".svg"
        shutil.copy(tree_file_path, tree_svg_dir)
  

def write_homoplastic_snp_table(report_data):
    # Initialize containers for the data
    all_snps_data = []
    core_snps_data = []
    majority_snps_data = []
    # Iterate over the dictionary
    for key, value in report_data.items():
        # get all SNPs
        if "COUNT_Homoplastic_SNPs.SNPs_all." in key:
            method = key.split('.')[-1]
            homoplastic_count = value[0]['Number_Homoplastic_SNPs']
            if method == "parsimony":
                method = "Parsimony"
            if method == "ML":
                method = "Maximum Likelihood"
            if method == "NJ":
                method = "Neighbor Joining"
            all_snps_data.append({"Method": method, "Number_Homoplastic_SNPs": homoplastic_count})
        # get core snps
        elif "COUNT_Homoplastic_SNPs.core_SNPs." in key:
            method = key.split('.')[-1]
            homoplastic_count = value[0]['Number_Homoplastic_SNPs']
            if method == "parsimony":
                method = "Parsimony"
            if method == "ML":
                method = "Maximum Likelihood"
            if method == "NJ":
                method = "Neighbor Joining"
            homoplastic_count = value[0]['Number_Homoplastic_SNPs']
            core_snps_data.append({"Method": method, "Number_Homoplastic_SNPs": homoplastic_count})
        elif "COUNT_Homoplastic_SNPs.SNPs_in_majority" in key:
            method = key.split('.')[-1]
            homoplastic_count = value[0]['Number_Homoplastic_SNPs']
            if method == "parsimony":
                method = "Parsimony"
            if method == "ML":
                method = "Maximum Likelihood"
            if method == "NJ":
                method = "Neighbor Joining"
            homoplastic_count = value[0]['Number_Homoplastic_SNPs']
            majority_snps_data.append({"Method": method, "Number_Homoplastic_SNPs": homoplastic_count})
    # Convert to DataFrames — skip any SNP set that produced no data
    merged_df = None
    if all_snps_data:
        merged_df = pd.DataFrame(all_snps_data).rename(columns={"Number_Homoplastic_SNPs": "All SNPs"})
    else:
        sys.stderr.write("No All SNP homoplastic data found; skipping All SNPs column.\n")
    if core_snps_data:
        df_c = pd.DataFrame(core_snps_data).rename(columns={"Number_Homoplastic_SNPs": "Core SNPs"})
        merged_df = pd.merge(merged_df, df_c, on="Method", how="outer") if merged_df is not None else df_c
    else:
        sys.stderr.write("No Core SNP homoplastic data found; skipping Core SNPs column.\n")
    if majority_snps_data:
        df_m = pd.DataFrame(majority_snps_data).rename(columns={"Number_Homoplastic_SNPs": "Majority SNPs"})
        merged_df = pd.merge(merged_df, df_m, on="Method", how="outer") if merged_df is not None else df_m
    if merged_df is None or merged_df.empty:
        return "<p>No homoplastic SNP data available.</p>"
    homoplastic_snps_html = generate_table_html_2(merged_df, table_width='75%')
    return homoplastic_snps_html


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
        # newick to phyloxml driver cocd
        clean_tree_dir = os.path.join(work_dir, "clean_trees")
        os.makedirs(clean_tree_dir, exist_ok=True)
        if "tree" in filename:
            clean_nwk_path = os.path.join(clean_tree_dir, filename)
            edit_newick_genome_id(file_path, clean_nwk_path)
            run_newick_to_phyloxml(clean_nwk_path)

@cli.command()
@click.argument("service_config")
def run_tree_to_svg(service_config):
    """Convert static svg images of nine basic trees for the report"""
    with open(service_config) as file:
        data = json.load(file)
        majority_threshold = data["params"]["majority-threshold"]
    ### start organize output files ###
    work_dir = data["work_data_dir"]
    output_dir = data["output_data_dir"]
    # make the output file directory
    tree_svg_dir = os.path.join (output_dir, "report_supporting_documents") 
    os.makedirs(tree_svg_dir, exist_ok=True)   
    if not os.path.exists(work_dir):
        sys.stderr.write("Work directory, {}, does not exist".format(work_dir))
        return
    tree_filenames = [
                "tree.SNPs_all.ML.tre",
                "tree.SNPs_all.NJ.tre",
                "tree.SNPs_all.parsimony.tre",
                "tree.core_SNPs.ML.tre",
                "tree.core_SNPs.NJ.tre",
                "tree.core_SNPs.parsimony.tre",
                "tree.SNPs_in_majority{}.ML.tre".format(majority_threshold),
                "tree.SNPs_in_majority{}.NJ.tre".format(majority_threshold),
                "tree.SNPs_in_majority{}.parsimony.tre".format(majority_threshold)
                ]
    for tree_filename in tree_filenames:
        file_path = os.path.join(work_dir, tree_filename)
        if os.path.exists(file_path) == True and os.path.getsize(file_path) > 0:
            run_p3x_tree_to_svg(file_path, tree_svg_dir)
        else:
            sys.stderr.write("service did not generate {}".format(tree_filename))


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
def find_optimum_k(kchooser_report):
    """ Parse the kChooser report for the optimum K value for the kSNP4 command."""
    parse_optimum_k(kchooser_report)

@cli.command()
@click.argument("service_config")
@click.argument("html_report_path")
def write_html_report(service_config, html_report_path):
    """Write an interactive report summarizing all outputs"""
    # run the functions here 
    report_data = {}
    with open(service_config) as file:
        data = json.load(file)
    clean_data_dir = data["clean_data_dir"]
    work_dir = data["work_data_dir"]
    majority_threshold = data["params"]["majority-threshold"]
    metadata_json = os.path.join(os.getcwd(), "genome_metadata.json")
    
    kchooser_report = os.path.join(clean_data_dir, "Kchooser4_ksnp4_input_file.report")
    report_data = parse_kchooser_report(report_data, kchooser_report)
    kchooser_df = pd.DataFrame.from_dict(report_data["kchooser_report"])
    
    report_data = parse_intermediate_files(report_data, work_dir)
    homoplastic_snps_html = write_homoplastic_snp_table(report_data)
    barplot_html = create_genome_length_bar_plot(clean_data_dir)
    snp_distribution_html = make_genome_bar_chart(data, report_data, majority_threshold)
    input_genome_table = generate_table_html_2(kchooser_df, table_width='75%')
    # SNP Counts 
    heatmap_html, metadata_json_string = interactive_threshold_heatmap(service_config, metadata_json, majority_threshold)
    html_template = define_html_template(input_genome_table, barplot_html, snp_distribution_html, \
                    homoplastic_snps_html, heatmap_html, \
                    majority_threshold, metadata_json_string)
    with open(html_report_path, 'w') as file:
        file.write(html_template)
    sys.stderr.write("Generated HTML report at {}.".format(html_report_path))
    print("let's go")


if __name__ == "__main__":
    cli()
