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

# Vérification des paquets libthrust-dev et libcub-dev (Debian/Ubuntu)
if command -v dpkg &>/dev/null; then
    for pkg in libthrust-dev libcub-dev; do
        if dpkg -l | grep -q "^ii  $pkg "; then
            echo "$pkg is installed. ok !"
        else
            echo "$pkg is not installed. Please install it with sudo apt install $pkg"
        fi
    done
else
    echo "Unable to check packages, dpkg not available."
fi
