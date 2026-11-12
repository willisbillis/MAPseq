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

echo "Checking to see if kallisto | bustools (kb) is installed..."
find_install kb

echo "Adding mapseq_functions.sh to /etc/profile.d/ to enable bash functions from this repository..."
if ! profile_tmp=$(mktemp); then
    echo "Error: Unable to prepare MAPseq shell functions." >&2
    exit 1
fi
trap 'rm -f -- "$profile_tmp"' EXIT
if ! { cat "$PWD/src/mapseq_functions.sh" && printf '\nexport MAPSEQ_REPO_PATH=%q\n' "$PWD"; } > "$profile_tmp" ||
   ! sudo install -m 644 "$profile_tmp" /etc/profile.d/mapseq_functions.sh; then
    echo "Error: Unable to install MAPseq shell functions." >&2
    exit 1
fi

echo "To create a new directory for a single run, use the command below with an appropriately titled run name."
echo -e "\ncreate_ms_run <RUN_NAME>\n"
echo -e "This directory (except for configuration and reference files) will mirror the local repository version of this pipeline.\n"
echo "To create a new directory for aggregating all completed runs, use the command below with an appropriately titled run name."
echo -e "\ncreate_ms_aggr_run <AGGR_RUN_NAME>\n"
echo "This directory (except for configuration and reference files) will mirror the local repository version of this pipeline, too."
echo -e "\n\nInstallation complete! Open a new login shell or run 'source /etc/profile.d/mapseq_functions.sh' for changes to take effect."