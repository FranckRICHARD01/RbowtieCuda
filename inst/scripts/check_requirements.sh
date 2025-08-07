#!/usr/bin/env bash

if ! command -v dpkg &>/dev/null; then
    echo "Error: 'dpkg' not found. Version comparisons may be unreliable." >&2
fi

if command -v dpkg &>/dev/null; then
    for pkg in libthrust-dev libcub-dev cmake; do
        if dpkg -l | grep -q "^ii  $pkg "; then
            echo "$pkg is installed"
        else
            echo "$pkg is missing; install with: sudo apt install $pkg"
        fi
    done
else
    echo "Cannot verify packages: dpkg unavailable"
fi

nvcc_path=$(command -v nvcc 2>/dev/null)
[ -z "$nvcc_path" ] && [ -x "/usr/local/cuda/bin/nvcc" ] && nvcc_path=/usr/local/cuda/bin/nvcc

nvidia_smi_path=$(command -v nvidia-smi 2>/dev/null)
[ -z "$nvidia_smi_path" ] && [ -x "/usr/bin/nvidia-smi" ] && nvidia_smi_path=/usr/bin/nvidia-smi

if [ -n "$nvcc_path" ] && [ -n "$nvidia_smi_path" ]; then
    echo "CUDA detected"
else
    echo "CUDA not found or tools missing; please run installer"
fi

full_cuda_ver=$("$nvcc_path" --version | grep release | awk '{print $NF}' | tr -d 'V,')
cuda_short=$(echo "$full_cuda_ver" | cut -d. -f1,2)

driver_cuda_ver=$("$nvidia_smi_path" | grep -oP 'CUDA Version: \K[0-9]+\.[0-9]+')

echo
echo "Detected CUDA compiler version: '$cuda_short'"
echo "Driver supports up to CUDA:        '$driver_cuda_ver'"

ver_cmp() {
    if command -v dpkg &>/dev/null; then
        dpkg --compare-versions "$1" "$2" "$3"; return
    fi
    awk -v a="$1" -v b="$3" 'function cmp(x,y,   i,X,Y,n,m){n=split(x,X,".");m=split(y,Y,".");for(i=1;i<=n||i<=m;i++){if((X[i]+0)<(Y[i]+0))return -1; if((X[i]+0)>(Y[i]+0))return 1;}return 0}
        BEGIN{exit !(cmp(a,b) '"$2"' 0)}'
}

if ver_cmp "$cuda_short" gt "$driver_cuda_ver"; then
    echo "*** Incompatibility: compiler $cuda_short vs driver $driver_cuda_ver ***"
    echo "*** Please downgrade CUDA or update NVIDIA driver ***"
else
    echo "Compiler '$cuda_short' is compatible with driver"
fi

full_gcc_ver=$(gcc --version | head -n1 | awk '{print $3}')
gcc_short=$(echo "$full_gcc_ver" | cut -d. -f1,2)

declare -A cuda_gcc_map=(
    ["13.0"]="15.9" ["12.8"]="14.9" ["12.9"]="14.9"
    ["12.4"]="13.2" ["12.5"]="13.2" ["12.6"]="13.2" ["12.7"]="13.2"
    ["12.1"]="12.2" ["12.2"]="12.2" ["12.3"]="12.2" ["12.0"]="12.1"
    ["11.4"]="11.0" ["11.5"]="11.0" ["11.6"]="11.0" ["11.7"]="11.0" ["11.8"]="11.0"
    ["11.0"]="9.0"  ["11.1"]="9.0"  ["11.2"]="9.0"  ["11.3"]="9.0"
)

expected_gcc=${cuda_gcc_map[$cuda_short]}

echo
echo "Detected GCC version: '$gcc_short'"
if [ -n "$expected_gcc" ]; then
    echo "Max recommended GCC for CUDA '$cuda_short' is '$expected_gcc'"
    if ver_cmp "$gcc_short" le "$expected_gcc"; then
        echo "GCC '$gcc_short' is compatible"
    else
        echo "*** Incompatibility: CUDA '$cuda_short' needs GCC ≤ '$expected_gcc', found '$gcc_short' ***"
        maj=${expected_gcc%%.*}; ((maj--))
        echo "Install GCC $maj or lower"
        if [ -f /etc/debian_version ]; then
            echo "On Debian/Ubuntu: sudo apt install gcc-$maj g++-$maj"
            echo "Use update-alternatives to switch versions"
        elif [ -f /etc/redhat-release ]; then
            echo "On RHEL/CentOS/Fedora: sudo dnf install gcc-$maj g++-$maj"
        elif [ -f /etc/arch-release ]; then
            echo "On Arch: sudo pacman -S gcc$maj"
        else
            echo "Refer to your distro’s package manager for GCC $maj"
        fi
    fi
else
    echo "No GCC compatibility entry for CUDA $cuda_short"
fi

echo
echo "nvcc path:       $nvcc_path"
echo "nvidia-smi path: $nvidia_smi_path"
echo 
"$nvidia_smi_path"
