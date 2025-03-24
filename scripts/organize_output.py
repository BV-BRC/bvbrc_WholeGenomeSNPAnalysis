import os
import shutil

# Load the JSON data
with open('work/config.json') as f:
    data = json.load(f)
source_folder = data
destination_folder = data["output_data_dir"]


# organize by
#trees
    # ML maximum liklihood tree 
    # and NJ Neighbor jonining tree
# SNPS
#   SNPS/ML
#   SNPS/NJ
#   SNPS/node
# counts
    # 
# # Ensure destination folders exist
# os.makedirs(f"{destination_folder}/ML", exist_ok=True)
# os.makedirs(f"{destination_folder}/NJ", exist_ok=True)

# # Loop through files in the source folder
# for file in os.listdir(source_folder):
#     file_path = os.path.join(source_folder, file)

#     if "ML" in file:
#         shutil.move(file_path, f"{destination_folder}/ML/{file}")
#     elif "NJ" in file:
#         shutil.move(file_path, f"{destination_folder}/NJ/{file}")

# if not empty move annotated list 



## elsewhere parse kchooser report for the report 