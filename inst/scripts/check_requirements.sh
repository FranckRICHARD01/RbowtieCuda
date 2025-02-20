#!/bin/bash

# Checking libthrust-dev, libcub-dev packages and cmake (Debian/Ubuntu)
if command -v dpkg &>/dev/null; then
    for pkg in libthrust-dev libcub-dev cmake; do
        if dpkg -l | grep -q "^ii  $pkg "; then
            echo "$pkg is installed. ok !"
        else
            echo "$pkg is not installed. Please install it with sudo apt install $pkg"
        fi
    done
else
    echo "Unable to check packages, dpkg not available."
fi

# Check if nvcc is available
nvcc_path=$(command -v nvcc 2>/dev/null)

# If nvcc is not found in the PATH, try locating in common directories
if [ -z "$nvcc_path" ]; then
    if [ -x "/usr/local/cuda/bin/nvcc" ]; then
        nvcc_path="/usr/local/cuda/bin/nvcc"
    fi
fi

# Check if nvidia-smi is available
nvidia_smi_path=$(command -v nvidia-smi 2>/dev/null)

# If nvidia-smi is not found in the PATH, try locating in common directories
if [ -z "$nvidia_smi_path" ]; then
    if [ -x "/usr/bin/nvidia-smi" ]; then
        nvidia_smi_path="/usr/bin/nvidia-smi"
    fi
fi

# If nvidia-smi and nvcc is found
if [ -n "$nvcc_path" ] && [ -n "$nvidia_smi_path" ]; then
    echo "CUDA is installed !"   
else
    echo "CUDA is not installed, nvcc or nvidia-smi are not accessibles. Please use the installation script..."
fi

# nvcc and nvidia-smi paths
echo "nvcc path: $nvcc_path"
echo "nvidia-smi path: $nvidia_smi_path"
echo ""
echo "nvidia-smi output:"
nvidia-smi


