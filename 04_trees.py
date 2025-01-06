import shutil
import os
import subprocess
from datetime import datetime

# Get the current date in YYYYMMDD format
tdate = datetime.now().strftime("%Y%m%d")

trimmed_folder = f"/Users/capelastegui.f/git/bovine_tree_figure/Data/03_trimmed_cleaned/{tdate}"
tree_folder = "/Users/capelastegui.f/git/bovine_tree_figure/Data/04_trees"

output_dir = f"/Users/capelastegui.f/git/bovine_tree_figure/Data/04_trees/{tdate}"

os.makedirs(output_dir, exist_ok=True)

for file_name in os.listdir(trimmed_folder):
    if file_name.endswith("clean.fasta") or file_name.endswith("clean.fa"):  # Only process fasta files
        # Construct the full file path
        file_path = os.path.join(trimmed_folder, file_name)
        shutil.copyfile(f"/Users/capelastegui.f/git/bovine_tree_figure/Data/03_trimmed_cleaned/{tdate}/{file_name}",
                        f"/Users/capelastegui.f/git/bovine_tree_figure/Data/04_trees/{tdate}/{file_name}")


# =============================================================================
# Run iqtree from python
# =============================================================================

def run_iqtree(input_fasta):
    """
    Runs IQTree on a FASTA file and saves the result to a new file.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to the output aligned FASTA file.
    """
    # Construct the MAFFT command. --auto lets MAFFT decide the best algorithm
    output_path = os.path.splitext(input_fasta)[0] + "_tree"  # Output filename

    iqtree_command = ["/Users/capelastegui.f/git/bovine_tree_figure/Data/04_trees/iqtree-2.2.2.6-MacOSX/bin/iqtree2",
                      "-s",
                      input_fasta,
                      "-redo"
                      ]
    print("IQTree made for")

    try:
        subprocess.run(iqtree_command, check=True)
        print(f"Successfully ran IQTree on {input_fasta}. Output saved to {output_path}.")
    except subprocess.CalledProcessError as e:
        print(f"Error running IQTree on {input_fasta}: {e}")


trimmed_alignment_folder = f"/Users/capelastegui.f/git/bovine_tree_figure/Data/04_trees/{tdate}/"

# Loop through the files in the alignment folder
for file_name in os.listdir(trimmed_alignment_folder):
    if file_name.endswith("clean.fasta") or file_name.endswith("clean.fa"):  # Only process fasta files
        # Construct the full file path
        base_name, ext = os.path.splitext(file_name)
        first_part = base_name.split('_')[0]
        input_fasta = os.path.join(trimmed_alignment_folder, file_name)

        run_iqtree(input_fasta)

