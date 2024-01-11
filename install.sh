#!/bin/bash

mkdir -p logs
INSTALL_OUTPUT_FILE=logs/install_output.log

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
cellranger testrun --id=tiny >> $INSTALL_OUTPUT_FILE
rm -r tiny

echo "Checking to see if Cell Ranger ATAC is installed and working..."
find_install cellranger-atac
cellranger-atac testrun --id=tiny >> $INSTALL_OUTPUT_FILE
rm -r tiny

echo "Checking to see if kallisto is installed..."
find_install kallisto

echo "Checking to see if bustools is installed..."
find_install bustools

if [ -z "$(echo $PATH | grep -o ~/.local/bin)" ]; then
    echo "Error: please add user local bin to PATH with the command below"
    echo "export PATH=\$PATH:~/.local/bin >> ~/.bashrc && source ~/.bashrc"
    return 1
fi
echo "Copying src/create_run.sh functions to ~/.local/bin to add to PATH..."
chmod +x src/create_run.sh
mkdir -p ~/.local/bin
cp src/create_run.sh ~/.local/bin