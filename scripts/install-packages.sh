#!/usr/bin/env bash
# Install system packages required to bootstrap vcpkg and build all dependencies.
set -euo pipefail

# -- vcpkg bootstrap requirements --
VCPKG_PACKAGES=(
    curl
    zip
    unzip
    tar
)

# -- Build-time requirements for vcpkg ports (gmsh, opencascade, etc.) --
PORT_PACKAGES=(
    autoconf
    autoconf-archive
    automake
    libtool
)

# -- X11/display libraries required by OpenCASCADE (gmsh[occ]) --
OCC_PACKAGES=(
    libx11-dev
    libxext-dev
    libxi-dev
    libxmu-dev
    libxt-dev
    libfreetype-dev
    libfontconfig1-dev
)

ALL_PACKAGES=("${VCPKG_PACKAGES[@]}" "${PORT_PACKAGES[@]}" "${OCC_PACKAGES[@]}")

echo "Updating package lists..."
sudo apt-get update -qq

echo "Installing packages: ${ALL_PACKAGES[*]}"
sudo apt-get install -y "${ALL_PACKAGES[@]}"

echo "All packages installed successfully."
