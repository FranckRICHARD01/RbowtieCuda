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

# system update
sudo apt-get update
sudo apt-get upgrade

# install headers
sudo apt install linux-headers-$(uname -r) 

# install other import packages
sudo apt-get install g++ freeglut3-dev build-essential libx11-dev libxmu-dev libxi-dev libglu1-mesa libglu1-mesa-dev libthrust-dev libcub-dev

# NVidia repository
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/cuda-ubuntu2404.pin
sudo mv cuda-ubuntu2404.pin /etc/apt/preferences.d/cuda-repository-pin-600
sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/3bf863cc.pub
sudo add-apt-repository "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/ /"
sudo apt-get update

# installing CUDA
sudo apt install nvidia-cuda-dev nvidia-cuda-toolkit

# setup your paths
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
successfully using the verif_installation_cuda.sh script.'

