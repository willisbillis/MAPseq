#!/bin/bash

# Check if user running uninstall script can write to their ~/.bashrc
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
    echo "Found existing installation (located at $MAPSEQ_REPO_PATH)."
    grep -v "^# MAPSEQ INSTALLATION\|^export MAPSEQ_REPO_PATH\|^source $MAPSEQ_REPO_PATH" ~/.bashrc > ~/.bashrc.tmp
    mv ~/.bashrc.tmp ~/.bashrc
    source ~/.bashrc
    echo "Uninstallation complete."
else
    echo "No installation found (environment variable MAPSEQ_REPO_PATH not set). Exiting..."
    exit 1
fi