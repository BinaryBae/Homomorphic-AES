#!/bin/bash

# Exit on any error
set -e

# Step 0: Delete existing build
rm -rf build bin

# Step 1: Configure the project with CMake
cmake -S . -B ./build

# Step 2: Build the project
cd build
make

# Step 3: Run the executable
../bin/Homo-AES
