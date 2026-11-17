#!/usr/bin/env bash
#
# Build-time environment diagnostic, called by src/nvbio/CMakeLists.txt.
#
# CONTRACT with CMake (see issue #4): this script must exit non-zero ONLY for
# things that genuinely block COMPILATION:
#   - no CUDA compiler (nvcc)
#   - missing build packages (cmake, libthrust-dev, libcub-dev)
# A missing, failing or GPU-less driver is a RUN-TIME requirement only.
# Build-only hosts (CI runners, r-universe, GitHub Actions, cluster login
# nodes, container builds) legitimately have no GPU: CMake then falls back to
# the oldest supported SASS + PTX JIT. Never exit non-zero for GPU/driver
# problems, and always end with an explicit `exit 0`.

# The GPU/driver section must not abort the script: report and continue.
set +e

hard_error=0

if command -v dpkg &>/dev/null; then
    for pkg in libthrust-dev libcub-dev cmake; do
        if dpkg -l "$pkg" 2>/dev/null | awk '$1 == "ii" {found=1} END {exit !found}'; then
            echo "$pkg is installed"
        else
            echo "$pkg is missing; install with: sudo apt install $pkg"
            hard_error=1
        fi
    done
else
    echo "Cannot verify packages: dpkg unavailable (non-Debian system?). Skipping those checks."
fi

nvcc_path=$(command -v nvcc 2>/dev/null)
[ -z "$nvcc_path" ] && [ -x "/usr/local/cuda/bin/nvcc" ] && nvcc_path=/usr/local/cuda/bin/nvcc

if [ -z "$nvcc_path" ]; then
    echo "CUDA not found; please run installer."
    exit 1
fi

full_cuda_ver=$("$nvcc_path" --version | grep release | awk '{print $NF}' | tr -d 'V,')
cuda_short=$(echo "$full_cuda_ver" | cut -d. -f1,2)

# --- GPU / driver section: purely informational -----------------------------
nvidia_smi_path=$(command -v nvidia-smi 2>/dev/null)
[ -z "$nvidia_smi_path" ] && [ -x "/usr/bin/nvidia-smi" ] && nvidia_smi_path=/usr/bin/nvidia-smi

driver_cuda_ver=""
gpu_state="none"

if [ -z "$nvidia_smi_path" ]; then
    echo "Warning: no NVIDIA driver/nvidia-smi found. Building anyway; a GPU is required at run time."
else
    "$nvidia_smi_path" -L >/dev/null 2>&1
    gpu_status=$?
    if [ "$gpu_status" -ne 0 ]; then
        echo "nvidia-smi failed (exit $gpu_status): no usable GPU on this host (build-only is fine)."
        echo "A GPU is required at run time only; CMake will use the oldest supported SASS + PTX JIT fallback."
    else
        gpu_state="usable"
        # The 'CUDA Version:' banner lives in the full table, not in -L.
        driver_cuda_ver=$("$nvidia_smi_path" 2>/dev/null | grep -oP 'CUDA Version: \K[0-9]+\.[0-9]+' | head -n1)
    fi
fi

ver_cmp() {
    if command -v dpkg &>/dev/null; then
        dpkg --compare-versions "$1" "$2" "$3"; return
    fi
    awk -v a="$1" -v b="$3" 'function cmp(x,y,   i,X,Y,n,m){n=split(x,X,".");m=split(y,Y,".");for(i=1;i<=n||i<=m;i++){if((X[i]+0)<(Y[i]+0))return -1; if((X[i]+0)>(Y[i]+0))return 1;}return 0}
        BEGIN{exit !(cmp(a,b) '"$2"' 0)}'
}

echo
echo "Detected CUDA compiler version: '$cuda_short'"
echo "Driver supports up to CUDA:        '$driver_cuda_ver'"

if [ -z "$driver_cuda_ver" ]; then
    if [ -z "$nvidia_smi_path" ]; then
        echo "No NVIDIA driver detected: the driver/toolkit pairing cannot be checked here."
        echo "If you install one, pick a driver at least as new as CUDA '$cuda_short' (sudo apt install nvidia-driver-XXX)."
    else
        echo "Could not read the driver CUDA version (no usable GPU on this host): skipping that comparison."
    fi
elif ver_cmp "$cuda_short" gt "$driver_cuda_ver"; then
    echo "*** Incompatibility: compiler $cuda_short vs driver $driver_cuda_ver ***"
    echo "*** Please downgrade CUDA or update NVIDIA driver (this blocks RUN time, not the build) ***"
else
    echo "Compiler '$cuda_short' is compatible with driver"
fi

full_gcc_ver=$(gcc --version | head -n1 | awk '{print $3}')
gcc_short=$(echo "$full_gcc_ver" | cut -d. -f1,2)

# Max host GCC officially supported per CUDA toolkit. Unlisted CUDA versions
# (newer toolkits) are settled by the compile probe in src/nvbio/CMakeLists.txt
# instead of being guessed here.
declare -A cuda_gcc_map=(
    ["13.2"]="15.9" ["13.1"]="15.9" ["13.0"]="15.9" ["12.8"]="14.9" ["12.9"]="14.9"
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
        echo "nvcc still builds if it is given a supported host compiler (g++-NN);"
        echo "src/nvbio/CMakeLists.txt probes this and picks one automatically."
        maj=${expected_gcc%%.*}; ((maj--))
        echo "To do it by hand: install GCC $maj or lower and configure with"
        echo "  cmake -DCMAKE_CUDA_HOST_COMPILER=\$(command -v g++-$maj) .."
        if [ -f /etc/debian_version ]; then
            echo "On Debian/Ubuntu: sudo apt install gcc-$maj g++-$maj"
        elif [ -f /etc/redhat-release ]; then
            echo "On RHEL/CentOS/Fedora: sudo dnf install gcc-$maj g++-$maj"
        elif [ -f /etc/arch-release ]; then
            echo "On Arch: sudo pacman -S gcc$maj"
        else
            echo "Refer to your distro’s package manager for GCC $maj"
        fi
    fi
else
    echo "No compatibility entry for CUDA '$cuda_short': trusting the compile probe at configure time."
fi

echo
echo "nvcc path:       $nvcc_path"
echo "nvidia-smi path: ${nvidia_smi_path:-<not found>}"
echo "GPU usable now:  $gpu_state (only needed at run time)"

if [ "$hard_error" -ne 0 ]; then
    echo
    echo "*** Missing build dependency (see above): compilation cannot proceed. ***"
    exit 1
fi

echo "Build environment OK."
exit 0
