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

# Get the CUDA version (only major and minor)
cuda_version=$($nvcc_path --version | grep "release" | sed -E ' s/.*release ([0-9.]+),.*/\1/')
cuda_version_short=$(echo "$cuda_version" | cut -d'.' -f1,2)

# Get the driver-supported CUDA version
driver_cuda_version=$($nvidia_smi_path | grep -oP 'CUDA Version: \K[0-9]+\.[0-9]+')

echo ""
echo "Detected CUDA version: $cuda_version_short"
echo "Driver-supported CUDA version: $driver_cuda_version"

# Check for incompatibility between CUDA and the driver
if [[ $(echo "$cuda_version_short > $driver_cuda_version" | bc -l) -eq 1 ]]; then
    echo "⚠️ Incompatibility detected: CUDA $cuda_version_short was compiled, but the driver only supports up to CUDA $driver_cuda_version."
    echo "You should either downgrade CUDA to $driver_cuda_version or upgrade the NVIDIA driver."
fi

# Get the GCC version (only major and minor)
gcc_version=$(gcc --version | head -n1 | awk '{print $3}')
gcc_version_short=$(echo "$gcc_version" | cut -d'.' -f1,2)

# Define exact compatibility mappings (only using major.minor for CUDA and GCC)
declare -A cuda_gcc_compat=(
    ["13.0"]="15.9"
    ["12.8"]="14.9"
    ["12.9"]="14.9"
    ["12.4"]="13.2"
    ["12.5"]="13.2"
    ["12.6"]="13.2"
    ["12.7"]="13.2"
    ["12.1"]="12.2"
    ["12.2"]="12.2"
    ["12.3"]="12.2"
    ["12.0"]="12.1"
    ["11.4"]="11.0"
    ["11.5"]="11.0"
    ["11.6"]="11.0"
    ["11.7"]="11.0"
    ["11.8"]="11.0"
    ["11.0"]="9.0"
    ["11.1"]="9.0"
    ["11.2"]="9.0"
    ["11.3"]="9.0"
)

# Find the maximum supported GCC version

expected_gcc=${cuda_gcc_compat[$cuda_version_short]}

if [[ -n "$expected_gcc" ]]; then
    if [[ "$(echo "$gcc_version_short <= $expected_gcc" | bc -l)" -eq 1 ]]; then
        echo "Compatible: CUDA $cuda_version_short with GCC $gcc_version_short ok !"
    else
        echo "Incompatible: CUDA $cuda_version_short requires GCC $expected_gcc max, but found $gcc_version_short."
        echo "You should install GCC $expected_gcc or lower."

        # Suggest installation commands based on the OS
        if [[ -f /etc/debian_version ]]; then
            echo "On Debian/Ubuntu, you can install it with:"
            echo "  sudo apt install gcc-$expected_gcc g++-$expected_gcc"
            echo "Then, set it as the default version with:"
            echo "  sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-$expected_gcc 100"
            echo "  sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-$expected_gcc 100"
        elif [[ -f /etc/redhat-release ]]; then
            echo "On RHEL/CentOS/Fedora, you can install it with:"
            echo "  sudo dnf install gcc-$expected_gcc g++-$expected_gcc"
        elif [[ -f /etc/arch-release ]]; then
            echo "On Arch Linux, you can install it with:"
            echo "  sudo pacman -S gcc$expected_gcc"
        else
            echo "Check your package manager for installing GCC $expected_gcc."
        fi
    fi
else
    echo "No compatibility information found for CUDA $cuda_version_short. Check NVIDIA specifications."
fi

# nvcc and nvidia-smi paths
echo ""
echo "nvcc path: $nvcc_path"
echo "nvidia-smi path: $nvidia_smi_path"
echo ""
echo "nvidia-smi output:"
nvidia-smi


