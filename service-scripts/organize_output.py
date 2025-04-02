import json
import os
import re
import shutil
import sys


def organize_files_by_type(source_dir, destination_dir):
    """
    Organize files by type based on the first word.
    """
    if not os.path.exists(source_dir):
        sys.stderr.write("Work directory, {}, does not exist".format(source_dir))
        return
    not_copied = []
    for filename in os.listdir(source_dir):
        file_path = os.path.join(source_dir, filename)
        intermediate_dir = os.path.join(destination_dir, "Intermediate_Files")
        os.makedirs(intermediate_dir, exist_ok=True)
        firstword = filename.split("_")[0]

        # Skip directories
        if not os.path.isfile(file_path):
            continue
        
        # Sort files into directories according to the first word
        if firstword == "All_SNPs" or firstword == "nonCore":
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
        if firstword == "tree":
            tree_dir = os.path.join(destination_dir, "Trees")
            os.makedirs(tree_dir, exist_ok=True)
            shutil.copy(file_path, tree_dir)
        if firstword.startswith("VCF"):
            VCFs_dir = os.path.join(destination_dir, "VCFs")
            os.makedirs(VCFs_dir, exist_ok=True)
            shutil.copy(file_path, VCFs_dir)


# driver code
# Load the JSON data
with open("config.json") as file:
    data = json.load(file)
### start organize output files ###
source_dir = data["work_data_dir"]
destination_dir = data["output_data_dir"]
organize_files_by_type(source_dir, destination_dir)
### end organize output files ###