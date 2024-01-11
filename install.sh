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
echo "Appending 'source $PWD/src/create_run.sh' to ~/.bashrc to add bash functions in this repository..."
chmod +x src/create_run.sh
echo "MAPSEQ_REPO_PATH=$PWD" >> ~/.bashrc
echo "source $MAPSEQ_REPO_PATH/src/create_run.sh" >> ~/.bashrc
source ~/.bashrc
echo "Creating dummy project config file in LRAXXX_template for example..."
cat $MAPSEQ_REPO_PATH/src/project_config_header.txt >> $1/project_config.txt
echo "PROJECT_NAME=LRAXXX" >> $1/project_config.txt
echo "PROJECT_PATH=$PWD/LRAXXX" >> $1/project_config.txt
cat $MAPSEQ_REPO_PATH/src/project_config_body.txt >> $1/project_config.txt

echo -e "Install complete! \n\n"
echo "To create a new directory for a single run, use the command below with an appropriately named run name."
echo -e "\ncreate_run <RUN_NAME>\n"
echo -e "This directory (except for configuration and reference files) will mirror the local repository version of this pipeline.\n"
echo "To create a new directory for aggregating all completed runs, use the command below with an appropriately named run name."
echo -e "\ncreate_aggr_run <AGGR_RUN_NAME>\n"
echo "This directory (except for configuration and reference files) will mirror the local repository version of this pipeline, too."