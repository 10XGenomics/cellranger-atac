#!/bin/bash

#
# Copyright (c) 2020 10x Genomics, Inc. All rights reserved.
#
# Environment setup for package cellranger-atac-1.2.0.
# Source this file before running.
#

# Determine path to this script; resolve symlinks
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
    DIR="$( cd -P "$( dirname "$SOURCE" )" > /dev/null && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
DIR="$( cd -P "$( dirname "$SOURCE" )" > /dev/null && pwd )"

#
# Source user's own environment first.
#
# Only source .bashrc if we're being sourced from shell10x script.
# Otherwise, we could end up in an infinite loop if user is
# sourcing this file from their .bashrc.
# Note: .bash_profile is for login shells only.
if [ ! -z $_RUN10X ] && [ -e ~/.bashrc ]; then
    source ~/.bashrc
fi


#
# Add module binary paths to PATH
#
export PATH="$DIR/bin:$PATH"
export PATH="$DIR/lib/bin:$PATH"
export PATH="$DIR/tenkit/bin:$PATH"
export PATH="$DIR/tenkit/lib/bin:$PATH"

#
# Module-specific env vars
#
export MROPATH="$DIR/tenkit/mro:$MROPATH"
export MROPATH="$DIR/mro/atac:$MROPATH"
export PYTHONPATH="$DIR/lib/cellranger/lib/python:$PYTHONPATH"
export PYTHONPATH="$DIR/lib/python:$PYTHONPATH"
export PYTHONPATH="$DIR/tenkit/lib/python:$PYTHONPATH"
