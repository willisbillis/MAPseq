#!/bin/bash

function find_install() {
    if [ -z "$(which $1)" ]; then
        echo "Error: $1 not found or installed!"
        return 1
    fi
}

echo "Checking to see if python is installed..."
find_install python

echo "Checking to see if Cell Ranger is installed and working..."
find_install cellranger
cellranger testrun --id=tiny

echo "Checking to see if Cell Ranger ATAC is installed and working..."
find_install cellranger-atac
cellranger-atac testrun --id=tiny

echo "Checking to see if kallisto is installed..."
find_install kallisto

echo "Checking to see if bustools is installed..."
find_install bustools

echo "Copying src/create_run.sh function to ~/bin to add to PATH..."
chmod +x src/create_run.sh
mkdir -p ~/bin
cp src/create_run.sh ~/bin