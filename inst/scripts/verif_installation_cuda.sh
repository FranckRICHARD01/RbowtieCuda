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
    # Get CUDA version
    version=$($nvcc_path --version | grep "release" | awk '{print $5}' | sed 's/,//')
    
    # Compare version with 12.4
    if [ "$(printf '%s\n' "12.4" "$version" | sort -V | head -n1)" == "12.4" ] && [ "$version" != "12.4" ]
    then
        echo "The CUDA version is higher than 12.4. This version is not compatible with RbowtieCuda. Please install a CUDA version lower or equal to 12.4. : $version"
    else
        echo "CUDA version less than or equal to 12.4. Congratulations! Your installation is ok! : $version"
    fi
else
    echo "CUDA is not installed or nvcc is not accessible. Please use the installation script..."
fi
