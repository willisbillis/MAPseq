#!/bin/bash

# Define the base directory where MAPseq directories are located
base_dir="/home/Projects/Scharer_sc/LRA.MAPseq"

# Find all directories matching the pattern {MAPseq_project}* on the first level
find "$base_dir" -maxdepth 1 -type d -name "LRA0[0-9][0-9]" | while read -r dir; do
  # Extract the MAPseq directory name
  mapseq_dir=$(basename "$dir")

  echo "Processing directory: $mapseq_dir"

  # Concatenate fastq files for RNA
  # Check if a matching {MAPseq_project}*sh directory exists
  shallow_dir="$base_dir/${mapseq_dir}sh/data/${mapseq_dir}sh_RNA/outs"
  if [[ -d "$shallow_dir" ]]; then
    # Define the full directory path
    full_dir="$dir/data/${mapseq_dir}_RNA/outs"

    echo "  Concatenating RNA fastq files..."

    # Check if orig.fqs directory exists
    if [[ -d "$full_dir/orig.fqs" ]]; then
      echo "    orig.fqs directory exists. Skipping concatenation for RNA."
      continue
    fi

    # Get all files in the full directory
    full_sample_names=($(ls $full_dir/*gz | awk -F '_S' '{print $1}' | sort -u))
    # Get basenames of full_sample_names
    full_sample_basenames=($(for name in "${full_sample_names[@]}"; do echo $(basename $name); done))

    # Get all files in the shallow directory
    shallow_sample_names=($(ls $shallow_dir/*gz | awk -F '_S' '{print $1}' | sort -u))
    # Get basenames of shallow_sample_names
    shallow_sample_basenames=($(for name in "${shallow_sample_names[@]}"; do echo $(basename $name); done))

    # Get the unique union of samples from both directories
    unset combined_samples
    combined_samples+=("${full_sample_basenames[@]}" "${shallow_sample_basenames[@]}")
    IFS=$' '
    union_samples=($(printf "%s\n" "${combined_samples[@]}" | sort -u | tr '\n' ' '))
    unset IFS

    for sample in "${union_samples[@]}"; do
      # Check if both full and shallow files exist
      full_file_r1="$full_dir/${sample}*_R1_001.fastq.gz"
      full_file_r2="$full_dir/${sample}*_R2_001.fastq.gz"
      shallow_file_r1="$shallow_dir/${sample}*_R1_001.fastq.gz"
      shallow_file_r2="$shallow_dir/${sample}*_R2_001.fastq.gz"

      # Resolve symlinks to get the original file paths
      full_file_r1=$(readlink -f $full_file_r1)
      full_file_r2=$(readlink -f $full_file_r2)
      shallow_file_r1=$(readlink -f $shallow_file_r1)
      shallow_file_r2=$(readlink -f $shallow_file_r2)

      # Concatenate multiple shallow files if necessary
      if [[ $(echo $shallow_file_r1 | wc -w) -gt 1 ]]; then
        # Concatenate all files in shallow_file_r1 into a single file
        for file in $shallow_file_r1; do
          cat "$file" >> "$shallow_dir/${sample}_S1_R1_001.fastq.gz"
        done
        shallow_file_r1="$shallow_dir/${sample}_S1_R1_001.fastq.gz"
      fi
      if [[ $(echo $shallow_file_r2 | wc -w) -gt 1 ]]; then
        # Concatenate all files in shallow_file_r2 into a single file
        for file in $shallow_file_r2; do
          cat "$file" >> "$shallow_dir/${sample}_S1_R2_001.fastq.gz"
        done
        shallow_file_r2="$shallow_dir/${sample}_S1_R2_001.fastq.gz"
      fi

      if [[ -f "$full_file_r1" && -f "$shallow_file_r1" && -f "$full_file_r2" && -f "$shallow_file_r2" ]]; then
        # Concatenate full and shallow files for each sample
        echo "    Concatenating sample: $sample"
        cat $full_file_r1 $shallow_file_r1 > "$full_dir/${sample}.mgd_S1_R1_001.fastq.gz"
        cat $full_file_r2 $shallow_file_r2 > "$full_dir/${sample}.mgd_S1_R2_001.fastq.gz"
      else
        echo "WARNING: Missing files for sample $sample in directory $mapseq_dir (RNA). Creating files from available data."
        # Check if shallow files exist and copy them if they do
        if [[ -f "$shallow_file_r1" && -f "$shallow_file_r2" ]]; then
          cp "$shallow_file_r1" "$full_dir/${sample}.mgd_S1_R1_001.fastq.gz"
          cp "$shallow_file_r2" "$full_dir/${sample}.mgd_S1_R2_001.fastq.gz"
        # Otherwise, copy the full files if they exist
        elif [[ -f "$full_file_r1" && -f "$full_file_r2" ]]; then
          cp "$full_file_r1" "$full_dir/${sample}.mgd_S1_R1_001.fastq.gz"
          cp "$full_file_r2" "$full_dir/${sample}.mgd_S1_R2_001.fastq.gz"
        fi
      fi

      # Create orig.fqs directory and move original fastq files
      mkdir -p "$full_dir/orig.fqs"
      mv "$full_dir"/${sample}*_R1_001.fastq.gz "$full_dir/orig.fqs"
      mv "$full_dir"/${sample}*_R2_001.fastq.gz "$full_dir/orig.fqs"

      # Rename mgd files
      mv $full_dir/orig.fqs/*.mgd* $full_dir
      rename 's/.mgd_S1/_S1/' $full_dir/*.mgd*.gz
    done
  fi

  # Concatenate fastq files for ATAC
  # Check if a matching {MAPseq_project}*sh directory exists
  shallow_dir="$base_dir/${mapseq_dir}sh/data/${mapseq_dir}sh_ATAC/outs"
  if [[ -d "$shallow_dir" ]]; then
    # Define the full directory path
    full_dir="$dir/data/${mapseq_dir}_ATAC/outs"

    echo "  Concatenating ATAC fastq files..."

    # Check if orig.fqs directory exists
    if [[ -d "$full_dir/orig.fqs" ]]; then
      echo "    orig.fqs directory exists. Skipping concatenation for ATAC."
      continue
    fi

    ## Get all files in the full directory
    full_sample_names=($(ls $full_dir/*gz | awk -F '_S' '{print $1}' | sort -u))
    # Get basenames of full_sample_names
    full_sample_basenames=($(for name in "${full_sample_names[@]}"; do echo $(basename $name); done))

    # Get all files in the shallow directory
    shallow_sample_names=($(ls $shallow_dir/*gz | awk -F '_S' '{print $1}' | sort -u))
    # Get basenames of shallow_sample_names
    shallow_sample_basenames=($(for name in "${shallow_sample_names[@]}"; do echo $(basename $name); done))

    # Get the unique union of samples from both directories
    combined_samples+=("${full_sample_basenames[@]}" "${shallow_sample_basenames[@]}")
    IFS=$' '
    union_samples=($(printf "%s\n" "${combined_samples[@]}" | sort -u | tr '\n' ' '))
    unset IFS

    for sample_path in "${union_samples[@]}"; do
      # Remove existing concatenated files if they exist
      if [[ -f "$shallow_dir/${sample}_S1_R1_001.fastq.gz" ]]; then
        rm "$shallow_dir/${sample}_S1_R1_001.fastq.gz"
      fi
      if [[ -f "$shallow_dir/${sample}_S1_R2_001.fastq.gz" ]]; then
        rm "$shallow_dir/${sample}_S1_R2_001.fastq.gz"
      fi
      if [[ -f "$shallow_dir/${sample}_S1_I1_001.fastq.gz" ]]; then
        rm "$shallow_dir/${sample}_S1_I1_001.fastq.gz"
      fi
      if [[ -f "$shallow_dir/${sample}_S1_I2_001.fastq.gz" ]]; then
        rm "$shallow_dir/${sample}_S1_I2_001.fastq.gz"
      fi

      # Check if both full and shallow files exist
      full_file_r1="$full_dir/${sample}*_R1_001.fastq.gz"
      full_file_r2="$full_dir/${sample}*_R2_001.fastq.gz"
      full_file_i1="$full_dir/${sample}*_I1_001.fastq.gz"
      full_file_i2="$full_dir/${sample}*_I2_001.fastq.gz"
      shallow_file_r1="$shallow_dir/${sample}*_R1_001.fastq.gz"
      shallow_file_r2="$shallow_dir/${sample}*_R2_001.fastq.gz"
      shallow_file_i1="$shallow_dir/${sample}*_I1_001.fastq.gz"
      shallow_file_i2="$shallow_dir/${sample}*_I2_001.fastq.gz"

      # Resolve symlinks to get the original file paths
      full_file_r1=$(readlink -f $full_file_r1)
      full_file_r2=$(readlink -f $full_file_r2)
      full_file_i1=$(readlink -f $full_file_i1)
      full_file_i2=$(readlink -f $full_file_i2)
      shallow_file_r1=$(readlink -f $shallow_file_r1)
      shallow_file_r2=$(readlink -f $shallow_file_r2)
      shallow_file_i1=$(readlink -f $shallow_file_i1)
      shallow_file_i2=$(readlink -f $shallow_file_i2)

      # Concatenate multiple shallow files if necessary
      if [[ $(echo $shallow_file_r1 | wc -w) -gt 1 ]]; then
        # Concatenate all files in shallow_file_r1 into a single file
        for file in $shallow_file_r1; do
          cat "$file" >> "$shallow_dir/${sample}_S1_R1_001.fastq.gz"
        done
        shallow_file_r1="$shallow_dir/${sample}_S1_R1_001.fastq.gz"
      fi
      if [[ $(echo $shallow_file_r2 | wc -w) -gt 1 ]]; then
        # Concatenate all files in shallow_file_r2 into a single file
        for file in $shallow_file_r2; do
          cat "$file" >> "$shallow_dir/${sample}_S1_R2_001.fastq.gz"
        done
        shallow_file_r2="$shallow_dir/${sample}_S1_R2_001.fastq.gz"
      fi
      if [[ $(echo $shallow_file_i1 | wc -w) -gt 1 ]]; then
        cat $shallow_file_i1 > "$shallow_dir/${sample}_S1_I1_001.fastq.gz"
        shallow_file_i1="$shallow_dir/${sample}_S1_I1_001.fastq.gz"
      fi
      if [[ $(echo $shallow_file_i2 | wc -w) -gt 1 ]]; then
        cat $shallow_file_i2 > "$shallow_dir/${sample}_S1_I2_001.fastq.gz"
        shallow_file_i2="$shallow_dir/${sample}_S1_I2_001.fastq.gz"
      fi

      if [[ -f "$full_file_r1" && -f "$shallow_file_r1" && -f "$full_file_r2" && -f "$shallow_file_r2" && -f "$full_file_i1" && -f "$shallow_file_i1" && -f "$full_file_i2" && -f "$shallow_file_i2" ]]; then
        # Concatenate full and shallow files for each sample
        echo "    Concatenating sample: $sample"
        cat $full_file_r1 $shallow_file_r1 > "$full_dir/${sample}.mgd_S1_R1_001.fastq.gz"
        cat $full_file_r2 $shallow_file_r2 > "$full_dir/${sample}.mgd_S1_R2_001.fastq.gz"
        cat $full_file_i1 $shallow_file_i1 > "$full_dir/${sample}.mgd_S1_I1_001.fastq.gz"
        cat $full_file_i2 $shallow_file_i2 > "$full_dir/${sample}.mgd_S1_I2_001.fastq.gz"
      else
        echo "WARNING: Missing files for sample $sample in directory $mapseq_dir (ATAC). Creating files from available data."
        # Check if shallow files exist and copy them if they do
        if [[ -f "$shallow_file_r1" && -f "$shallow_file_r2" && -f "$shallow_file_i1" && -f "$shallow_file_i2" ]]; then
          cp "$shallow_file_r1" "$full_dir/${sample}.mgd_S1_R1_001.fastq.gz"
          cp "$shallow_file_r2" "$full_dir/${sample}.mgd_S1_R2_001.fastq.gz"
          cp "$shallow_file_i1" "$full_dir/${sample}.mgd_S1_I1_001.fastq.gz"
          cp "$shallow_file_i2" "$full_dir/${sample}.mgd_S1_I2_001.fastq.gz"
        # Otherwise, copy the full files if they exist
        elif [[ -f "$full_file_r1" && -f "$full_file_r2" && -f "$full_file_i1" && -f "$full_file_i2" ]]; then
          cp "$full_file_r1" "$full_dir/${sample}.mgd_S1_R1_001.fastq.gz"
          cp "$full_file_r2" "$full_dir/${sample}.mgd_S1_R2_001.fastq.gz"
          cp "$full_file_i1" "$full_dir/${sample}.mgd_S1_I1_001.fastq.gz"
          cp "$full_file_i2" "$full_dir/${sample}.mgd_S1_I2_001.fastq.gz"
        fi
      fi
    
      # Create orig.fqs directory and move original fastq files
      mkdir -p "$full_dir/orig.fqs"
      mv "$full_dir"/${sample}*_001.fastq.gz "$full_dir/orig.fqs"
      mv "$full_dir"/${sample}*_001.fastq.gz "$full_dir/orig.fqs"

      # Rename mgd files
      mv $full_dir/orig.fqs/*.mgd* $full_dir
      rename 's/.mgd_S1/_S1/' $full_dir/*.mgd*.gz
    done
  fi
done
