import base64
import click
import glob
import json
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import re
import sys

from pathlib import Path

def add_to_report_dict(report_data, source_name, item):
    if source_name not in report_data:
        report_data[source_name] = []  # Initialize the list if the source doesn't exist
    report_data[source_name].append(item)
    return report_data

def define_html_template(input_genome_table, identified_snps, pie_charts, homoplastic_all_snps_html, homoplastic_core_snps_html, homoplastic_html, heatmap_html):
    core_snps_html = pie_charts[0]
    majority_snps_html = pie_charts[1]
    html_template = """
            <!DOCTYPE html>
            <html lang="en">
            <head>
                <meta charset="UTF-8">
                <meta name="viewport" content="width=device-width, initial-scale=1.0">
                <style>
                    .plot-container {{
                        display: flex;
                        justify-content: space-around;
                        }}
                    .plot {{
                            width: 50%; /* Each plot takes up half of the container width */
                        }}
                        .caption {{
                            text-align: center;
                            font-size: 14px;
                            margin-top: 20px;
                            color: black;
                            font-family: Roboto;
                        }}
                    body {{ font-family: Roboto, sans-serif; color: black; }}
                        header {{
                            display: flex;
                            justify-content: space-between;
                            align-items: center;
                            padding: 10px 20px;
                        }}
                        header > a img {{
                            max-width: 225px;  /* Maximum width */
                            max-height: 225px;  /* Maximum height */
                            width: auto;
                            height: auto;
                        }}
                        .title {{
                            font-size: 36px;  /* Adjust the size of the title text */
                            font-family: 'Roboto', sans-serif;
                            font-weight: bold;
                            color: black;
                        }}
                        .warning {{ color: black; }}
                        table, th, td {{ border: 1px solid black; border-collapse: collapse; }}
                        th, td {{ padding: 5px; text-align: left; }}
                        img {{ width: 100%; max-width: 600px; height: auto; }}
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
                            width: 100%; /* Make the image expand to fill the container */
                            max-width: 600px; /* Maximum width of the image */
                </style>
                <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
            </head>
            <body>
                <header>
                    <div class="title">Whole Genome SNP Analysis Report</div>
                        </a>
                </header>
            <p> 
            Explore SNP results more text here...NICOLE
            <h2>Getting to Know the Input Data</h2>
            <h3>Input Genomes</h3>
            <p>
            {input_genome_table}
            <p>
            <h3>Identified SNPs</h3> 
            <div class="plot-container">
                <div class="plot" id="plot1">
                    {core_snps_html} <!-- Direct embedding of the plot content -->
                </div>
                <div class="plot" id="plot2">
                    {majority_snps_html} <!-- Direct embedding of the plot content -->
                </div>
            </div>
            <h3>Homoplastic SNPs</h3> 
            {homoplastic_html}
            {homoplastic_all_snps_html}
            </div>
            {homoplastic_core_snps_html}
            </div>
            {heatmap_html}
        </body>
        </html>
        """.format(input_genome_table = input_genome_table, identified_snps=identified_snps, core_snps_html=core_snps_html, \
                majority_snps_html=majority_snps_html,  homoplastic_all_snps_html=homoplastic_all_snps_html, \
                homoplastic_core_snps_html=homoplastic_core_snps_html, homoplastic_html=homoplastic_html, heatmap_html=heatmap_html)
    return html_template

def generate_table_html_2(kchooser_df, table_width='95%'):
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
    <table style="width: {table_width}; border-collapse: collapse; border: 1px solid black;">
        <thead>
            <tr>{headers}</tr>
        </thead>
        <tbody>
            {rows}
        </tbody>
    </table>
    '''
    return table_html

def generate_plotly_homoplastic_table(df1, df2, table_titles=("Homoplastic SNPs in All SNPS", "Homoplastic SNPS in Core SNPs")):
    # Create subplots with 1 row and 2 columns
    fig = make_subplots(rows=1, cols=2, 
                        subplot_titles=table_titles,
                        specs=[[{"type": "table"}, {"type": "table"}]])

    # Generate Table 1
    fig.add_trace(go.Table(
        header=dict(values=list(df1.columns),
                    fill_color='lightgrey',
                    align='center',
                    font=dict(size=16, color='black')
                    ),
        cells=dict(values=[df1[col] for col in df1.columns],
                   fill_color='lightpink',
                   align='center',
                   font=dict(size=13, color='black')
            )
        ), row=1, col=1)

    # Generate Table 2
    fig.add_trace(go.Table(
        header=dict(values=list(df2.columns),
                    fill_color='lightgrey',
                    align='center',
                    font=dict(size=16, color='black')),
        cells=dict(values=[df2[col] for col in df2.columns],
                   fill_color='lightpink',
                   align='center',
                   font=dict(size=13, color='black')
            )
        ), row=1, col=2)


    # Update layout for better display
    fig.update_layout(height=400, width=1800, showlegend=False, 
                      title=dict(
            text="Side by Side Comparison of Homoplasstic SNPs",   # Main title text
            font=dict(size=20, color='black', family='Arial'),  # Font settings
            x=0.5,   # Center alignment
            xanchor='center'
        ),
        font=dict(size=14),
        annotations=[  # Styling subtitles
            dict(font=dict(size=16, color='black', family='Arial'))
        ]
    )
    fig.write_html("homoplastic_snp_table.html", include_plotlyjs=False)
    homoplastic_html = read_plotly_html("homoplastic_snp_table.html")
    return homoplastic_html

def make_count_piecharts(identified_snps_df):
    # CoreSNPS pie chart asks you to pass data as a list and categories as a list
    core_snps_df = identified_snps_df[["Number core SNPs", "Number non-core SNPs"]]
    core_snps_df.rename(
        columns={
            "Number core SNPs" : "Core SNPs Present in All Genomes",
            "Number non-core SNPs" : "Non-Core SNPs Not Present in All Genomes"
        },
        inplace=True
    )
    categories = core_snps_df.columns
    counts = core_snps_df.values
    counts = list(counts[0])

    fig = go.Figure()
    fig.add_trace(go.Pie(
    labels=categories,
    values=counts,
    hoverinfo='label+percent',
    textinfo='value',
    marker=dict(colors=["palevioletred","mistyrose"])
    ))
    fig.update_layout(height=600, width=800, title_text="Core SNPs")
    fig.write_html("core_snps_pie_chart.html", include_plotlyjs=False)

    # majority SNPS
    majority_snps_df = identified_snps_df.iloc[:, -2:]
    categories = majority_snps_df.columns
    counts = majority_snps_df.values
    counts = list(counts[0])
    non_core_snps = counts[1] - counts[0]
    counts = [counts[0], non_core_snps]
    fig = go.Figure()
    fig.add_trace(go.Pie(
    labels=categories,
    values=counts,
    hoverinfo='label+percent',
    textinfo='value',
    marker=dict(colors=["palevioletred","mistyrose"])
    ))
    fig.update_layout(height=600, width=800, title_text="SNPs in the Majority of Genomes")
    fig.write_html("majority_snp_pie_chart.html", include_plotlyjs=False) # change for prod look at wastewater to see how you made this work 
    
    core_snps_html = read_plotly_html("core_snps_pie_chart.html")
    majority_snps_html = read_plotly_html("majority_snp_pie_chart.html")
    return core_snps_html, majority_snps_html

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
        kchooser_data["Median Genome"] = match.group(1)

    match = re.search(r'Its length is (\d+)', text)
    if match:
        kchooser_data["Median Genome Length"] = int(match.group(1))
    # Extract shortest genome and its length
    match = re.search(r'The shortest genomes is (\S+) its length is (\d+)', text)
    if match:
        kchooser_data["Shortest Genome"] = match.group(1)
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

def parse_intermediate_files(report_data, work_dir):
    count_files = {}
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

def read_plotly_html(plot_path):
    # Read the content from 'Variant_Plot_Interactive.html'
    with open(plot_path, 'r') as file:
        plotly_html_content = file.read()
    # Extract everything within the <body> tags
    extracted_content = re.findall(r'<body>(.*?)</body>', plotly_html_content, re.DOTALL)

    # Assuming extracted_content contains our needed Plotly graph initialization scripts
    plotly_graph_content = extracted_content[0] if extracted_content else ''
    return plotly_graph_content

def write_homoplastic_snp_table(report_data):
    # Initialize containers for the data
    all_snps_data = []
    core_snps_data = []
    # Iterate over the dictionary
    for key, value in report_data.items():
        # get all SNPS
        if "COUNT_Homoplastic_SNPs.SNPs_all." in key:
            method = key.split('.')[-1]
            homoplastic_count = value[0]['Number_Homoplastic_SNPs']
            all_snps_data.append({"Method": method, "Number_Homoplastic_SNPs": homoplastic_count})
        # get core snps
        elif "COUNT_Homoplastic_SNPs.core_SNPs." in key:
            method = key.split('.')[-1]
            homoplastic_count = value[0]['Number_Homoplastic_SNPs']
            core_snps_data.append({"Method": method, "Number_Homoplastic_SNPs": homoplastic_count})

    # Create DataFrames
    df_all_snps = pd.DataFrame(all_snps_data)
    df_core_snps = pd.DataFrame(core_snps_data)
    print(df_all_snps)
    print(df_all_snps.columns)
    homoplastic_all_snps = generate_table_html_2(df_all_snps, table_width='95%')
    homoplastic_core_snps = generate_table_html_2(df_core_snps, table_width='95%')
    homoplastic_html = generate_plotly_homoplastic_table(df_all_snps, df_core_snps, table_titles=("Homoplastic SNPs in All SNPS", "Homoplastic SNPS in Core SNPs"))
    return homoplastic_all_snps, homoplastic_core_snps, homoplastic_html

@click.group()
def cli():
    """ This script supports the Whole Genome SNP service with multiple commands."""
    pass


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
    # destination_dir = data["output_data_dir"]
    kchooser_report = os.path.join(clean_data_dir, "Kchooser4_ksnp4_input_file.report")

    report_data = parse_kchooser_report(report_data, kchooser_report)

    report_data = parse_intermediate_files(report_data, work_dir)
    # print(report_data)

    # nicole here
    # which table is better?
    # homoplastic_html = write_homoplastic_snp_table(report_data)
    homoplastic_all_snps_html, homoplastic_core_snps_html, homoplastic_html = write_homoplastic_snp_table(report_data)

    # nicole new dev 
    # metadata table to df?
    # Load the JSON data
    with open('genome_metadata.json') as f:
        data = json.load(f)

    # Convert to DataFrame
    meta_data_df = pd.DataFrame(data)

    # Display the DataFrame
    print(meta_data_df)

    # heatmap_for_report.html
    heatmap_html = read_plotly_html("heatmap_for_report.html")

    # below this is good to go
    identified_snps_df = pd.DataFrame.from_dict(report_data["COUNT_coreSNPs"])
    identified_snps_df["Total Snps"] = report_data["COUNT_SNPs"][0]["Number_SNPs"]
    pie_charts = make_count_piecharts(identified_snps_df)
    kchooser_df = pd.DataFrame.from_dict(report_data["kchooser_report"])
    input_genome_table = generate_table_html_2(kchooser_df, table_width='95%')
    # snp counts as a table or pie
    identified_snps_table = generate_table_html_2(identified_snps_df, table_width='95%')
    html_template = define_html_template(input_genome_table, identified_snps_table, pie_charts, \
                    homoplastic_all_snps_html, homoplastic_core_snps_html, homoplastic_html, heatmap_html)
    with open(html_report_path, 'w') as file:
        file.write(html_template)
    sys.stderr.write("Generated HTML report at {}.".format(html_report_path))
    print("let's go")

### DEV TO DO###
"""
path to working directory 
/home/nbowers/bvbrc-dev/dev_container/holder_modules/tmp_dev_ksnp4/writing_report

path to this file /home/nbowers/bvbrc-dev/dev_container/modules/bvbrc_WholeGenomeSNPAnalysis/service-scripts/tmp_report.py

python3 /home/nbowers/bvbrc-dev/dev_container/modules/bvbrc_WholeGenomeSNPAnalysis/service-scripts/tmp_report.py write-html-report config.json test.html

1. check the import statements at the top
2. fix lable in core snp pie chart - try plotly subplots
3. fix decimal points in the genome IDs
4. if we go with plotly tables, fix the numbesr the same way as the other tables
5. which method is used to make the results in the heatmap from core snps from ksnpdist? NJ, ML, or parsimony?
6. !! add date column to the metadata table or find out why it ismissing???

LAST return to making the metadata options in the heatmap
"""
if __name__ == "__main__":
    cli()