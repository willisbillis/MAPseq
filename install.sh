#!/bin/bash

# Check if user running install script can write to their ~/.bashrc
bashrc_permissions=$(stat -c "%A" $(realpath ~/.bashrc))
# Do you own your .bashrc?
if ! [[ $(stat -c "%U" $(realpath ~/.bashrc)) == $USER ]]; then
  # If not, are you part of the group that owns your .bashrc?
  if ! [[ $(stat -c "%G" $(realpath ~/.bashrc)) == $(id -gn $USER) ]]; then
    # If not, are others able to write to your .bashrc? (unlikely)
    if ! [[ "${bashrc_permissions:8:1}" == "w" ]]; then
        echo "Unable to write to ~/.bashrc. Please modify permissions and try again."
        exit 1
    fi
  else
    # If you are part of the group that owns your .bashrc, are there group write priviledges? (also unlikely)
    if ! [[ "${bashrc_permissions:5:1}" == "w" ]]; then
        echo "Unable to write to ~/.bashrc. Please modify permissions and try again."
        exit 1
    fi
  fi
fi

if [[ -n ${MAPSEQ_REPO_PATH+x} ]]; then
    echo "Found existing installation (located at $MAPSEQ_REPO_PATH). Exiting..."
    exit 1
fi

mkdir -p logs
INSTALL_OUTPUT_FILE=logs/install_output.log

function find_install() {
    if [ -z "$(which $1)" ]; then
        echo "Error: $1 not found or installed!"
        exit 1
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

echo "Appending 'source $PWD/src/create_run.sh' to ~/.bashrc to add bash functions in this repository..."
echo -e "\n# MAPSEQ INSTALLATION EDITS - DO NOT MANUALLY CHANGE (Installed $(date))" >> ~/.bashrc
echo "export MAPSEQ_REPO_PATH=$PWD" >> ~/.bashrc
echo "source $PWD/src/create_run.sh" >> ~/.bashrc

echo -e "Installation complete! Restart shell (by closing and opening new command line session or running $SHELL) for changes to take effect. \n\n"
echo "To create a new directory for a single run, use the command below with an appropriately titled run name."
echo -e "\ncreate_ms_run <RUN_NAME>\n"
echo -e "This directory (except for configuration and reference files) will mirror the local repository version of this pipeline.\n"
echo "To create a new directory for aggregating all completed runs, use the command below with an appropriately titled run name."
echo -e "\ncreate_ms_aggr_run <AGGR_RUN_NAME>\n"
echo "This directory (except for configuration and reference files) will mirror the local repository version of this pipeline, too."