#!/bin/bash

# Check if nvcc is available
nvcc_path=$(command -v nvcc 2>/dev/null)

# If nvcc is not found in the PATH, try locating in common directories
if [ -z "$nvcc_path" ]; then
    if [ -x "/usr/local/cuda/bin/nvcc" ]; then
        nvcc_path="/usr/local/cuda/bin/nvcc"
    fi
fi

# If nvcc is found
if [ -n "$nvcc_path" ]; then
    echo "CUDA is installed !"   
else
    echo "CUDA is not installed or nvcc is not accessible. Please use the installation script..."
fi
