#!/bin/bash

### steps ####
# verify the system has a cuda-capable gpu
# download and install the nvidia cuda toolkit and cudnn
# setup environmental variables
# verify the installation
###

### to verify your gpu is cuda enable check
if ! lspci | grep -i nvidia; then
    echo "No NVIDIA GPU detected. Exiting."
    exit 1
fi

### continue ?
echo ""
echo "Do you want to delete all nvidia libraries and drivers before you start?
This allows you to start from scratch.
To be used only if previous executions of this script have failed to install cuda properly."

while true; do
    read -p "Do you wish to continue (y/n)? " yn
    case $yn in
        [Yy]* ) 
            echo "Removing old NVIDIA drivers and libraries..."
            sudo apt-get purge -y nvidia* libnvidia*
            sudo apt-get remove -y nvidia-*
            sudo rm -f /etc/apt/sources.list.d/*cuda*
            sudo apt-get autoremove -y && sudo apt-get autoclean -y
            sudo rm -rf /usr/local/cuda*
            echo "Old NVIDIA drivers and libraries have been removed."
            break
            ;;
        [Nn]* ) 
            break
            ;;
        * ) 
            echo "Please answer yes (Y) or no (N)."
            ;;
    esac
done

echo "Cuda installation..."

# system update
sudo apt-get update
sudo apt-get upgrade

# install headers
sudo apt install linux-headers-$(uname -r) 

# install other import packages
sudo apt install g++ freeglut3-dev build-essential libx11-dev libxmu-dev libxi-dev libglu1-mesa libglu1-mesa-dev 

echo "Choose an option:"
echo "1) install the default version of cuda on your distribution "
echo "2) install the version 12.8 of cuda (do not use if you're not sure)"

read -p "Your choice (1 or 2): " choice

case $choice in
    1)
        # installing default version of CUDA
        sudo apt install nvidia-cuda-dev nvidia-cuda-toolkit libthrust-dev libcub-dev
        ;;
    2)
        # NVidia repository
        wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/cuda-ubuntu2404.pin
        sudo mv cuda-ubuntu2404.pin /etc/apt/preferences.d/cuda-repository-pin-600
        sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/3bf863cc.pub
        sudo add-apt-repository "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/ /"
        sudo apt-get update

        # installing CUDA
        sudo apt install cuda-12-8 nvidia-driver-570 libthrust-dev libcub-dev
        ;;
    *)
        echo "Invalid choice. Please enter 1 or 2."
        ;;
esac


# setup your paths
echo ''
echo 'Please modify your .bashrc file located at the root of your account by adding the lines :'
echo '  export PATH=/usr/local/cuda-xx.x/bin:$PATH'
echo '  export LD_LIBRARY_PATH=/usr/local/cuda-xx.x/lib64:$LD_LIBRARY_PATH'
echo 'Replacing xx.x with the version of cuda installed (e.g. 12.6)'
echo ' '
echo 'then do :'
echo '  source ~/.bashrc'
echo '  sudo ldconfig'
echo ' '
echo 'Please reboot the machine and check that the installation has been completed 
successfully using the check_requirements.sh script.'

