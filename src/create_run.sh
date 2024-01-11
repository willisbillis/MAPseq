#!/bin/bash

function create_run() {
    if [ $# -ne 1 ]; then
        echo "Error: Exactly 1 arguments required!"
        return 1
    fi
    if [ -z "$1" ]; then
        echo "Error: Argument is empty!"
        return 1
    fi

    ln -s LRAXXX_template $1
}

function create_aggr_run() {
    if [ $# -ne 1 ]; then
        echo "Error: Exactly 1 arguments required!"
        return 1
    fi
    if [ -z "$1" ]; then
        echo "Error: Argument is empty!"
        return 1
    fi

    ln -s LRA_all_template $1
}