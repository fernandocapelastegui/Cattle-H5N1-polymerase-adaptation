#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
big Created on Thu Sep 26 21:56:40 2024

@author: capelastegui.f
"""
import os
from Bio import AlignIO
from pytrimal import AutomaticTrimmer, Alignment, ManualTrimmer
from sqlalchemy.dialects.mssql.information_schema import sequences
from Bio import SeqIO
from datetime import datetime

# Get the current date in YYYYMMDD format
tdate = datetime.now().strftime("%Y%m%d")

alignment_folder = f"/Users/capelastegui.f/git/bovine_tree_figure/Data/02_aligned_data/{tdate}"

all_alignments = {}

# Loop through the files in the alignment folder
for file_name in os.listdir(alignment_folder):
    if file_name.endswith(".fasta") or file_name.endswith(".fa"):  # Only process fasta files
        # Construct the full file path
        file_path = os.path.join(alignment_folder, file_name)

        # Read the alignment and store it in the dictionary with the filename as the key
        alignment = Alignment.load(file_path)  # Reads alignment as an Alignment object
        all_alignments[file_name] = alignment

# set the function for the trimming:

trimmer = ManualTrimmer(gap_threshold=0.85, conservation_percentage=70)

#Create an empty list for trimmed alignments in preparation
all_alignments_trimmed = {}

#This is where we want the trimmed alignment to go. Create the dir if it doesnt exist
output_dir = f"/Users/capelastegui.f/git/bovine_tree_figure/Data/03_trimmed_cleaned/{tdate}"

# Create the directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Loop through the dictionary of alignments by filename
for file_name, alignment in all_alignments.items():
    #trim the alignment
    trimmed = trimmer.trim(alignment)
    all_alignments_trimmed[file_name] = trimmed  # Store the trimmed alignment

    # Generate a new filename by appending '_trimmed' before the file extension
    base_name, ext = os.path.splitext(file_name)
    trimmed_filename = f"{base_name}_trimmed{ext}"

    # Construct the full output path
    output_path = os.path.join(f"/Users/capelastegui.f/git/bovine_tree_figure/Data/03_trimmed_cleaned/{tdate}",
                               trimmed_filename)

    # Export the trimmed alignment to a FASTA file
    trimmed.dump(output_path, format="fasta")


########
#cut out very missing sequences and further cleaning
########

trimmed_folder = f"/Users/capelastegui.f/git/bovine_tree_figure/Data/03_trimmed_cleaned/{tdate}"
cleaned_folder = f"/Users/capelastegui.f/git/bovine_tree_figure/Data/03_trimmed_cleaned/{tdate}"
# Create output directory if it doesn't exist


# Function to remove sequences with gaps or 'N's above a certain threshold
def remove_sequences_with_gaps(trimmed_alignment):
    filtered_sequences = []

    # Iterate over the records in the trimmed alignment directly
    for record in trimmed_alignment:  # Directly iterate over the alignment object
        # Count the number of missing positions (gaps or ambiguous positions like 'N')
        gap_count = record.seq.count("-")
        n_count = record.seq.count("X")
        missing_count = gap_count + n_count  # Sum both counts

        # If the total count of missing positions is below or equal to the threshold, keep the sequence
        if missing_count <=30:
            filtered_sequences.append(record)
    return filtered_sequences  # Return filtered sequences


# Loop through the files in the trimmed folder
for file_name in os.listdir(trimmed_folder):
    if file_name.endswith(".fasta") or file_name.endswith(".fa"):  # Only process fasta files
        # Construct the full file path
        file_path = os.path.join(trimmed_folder, file_name)
        print(f"Processing file: {file_path}")

        # Read the alignment
        alignment = AlignIO.read(file_path, "fasta")
        # Remove sequences based on gaps and N's
        cleaned_seqs = remove_sequences_with_gaps(alignment)
        # Remove duplicates by sequence content
        unique_seqs = []
        seen_seqs = set()
        for record in cleaned_seqs:
            seq_str = str(record.seq)  # Convert the sequence to a string for comparison
            if seq_str not in seen_seqs:
                seen_seqs.add(seq_str)
                unique_seqs.append(record)

        # Create a new filename by appending "_clean" before the file extension
        base_name, ext = os.path.splitext(file_name)  # Split file name and extension
        clean_file_name = f"{base_name}_clean{ext}"  # Append '_clean' to the base name
        clean_file_path = os.path.join(cleaned_folder, clean_file_name)

        # Write the cleaned sequences to a new file
        if cleaned_seqs:  # Check if there are sequences left to write
            SeqIO.write(unique_seqs, clean_file_path, "fasta")
            print(f"Saved cleaned sequences to: {clean_file_path}")
        else:
            print(f"No sequences left after filtering in {file_name}.")
