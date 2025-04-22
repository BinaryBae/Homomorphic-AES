# AES-TFHE Project

This project demonstrates the integration of AES encryption with the [TFHE (Torus Fully Homomorphic Encryption)](https://github.com/tfhe/tfhe) library using CMake for build configuration.
+
## Group Details
- **Patel Ayush Pravinkumar** (21114071)  
- **Rudransh** (21114085)  
- **Sanidhya Singh** (21114091)  
- **Shrey Gupta** (21112103)  

---

## Requirements

Make sure the following dependencies are installed before building the project:

- **CMake** (version 3.10 or above)
- **TFHE library**

### Installing Requirements

#### Install CMake

On Ubuntu/Debian:

```sh
sudo apt update
sudo apt install cmake
```
####  Install TFHE

Clone and build the TFHE library:

```sh
git clone https://github.com/tfhe/tfhe.git
cd tfhe
mkdir build && cd build
cmake ..
make
sudo make install
```
## Building and Running

Once all dependencies are installed, use the provided bash script to build and run the project.

### Run the Project

```sh
./build_and_run.sh
```
This script performs the following:

1. Configures the project using CMake.
2. Builds the project in a `build/` directory and generate the executable file in a `bin/` directory.
3. Runs the compiled executable located at `bin/Home-AES`.

---


