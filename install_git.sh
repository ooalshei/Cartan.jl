#!/bin/bash

if command -v git &> /dev/null; then
    exit 0
fi

detect_package_manager() {
    if command -v apt-get &> /dev/null; then
        echo "apt-get"
    elif command -v dnf &> /dev/null; then
        echo "dnf"
    elif command -v yum &> /dev/null; then
        echo "yum"
    elif command -v pacman &> /dev/null; then
        echo "pacman"
    elif command -v brew &> /dev/null; then
        echo "brew"
    else
        echo "unknown"
    fi
}

echo "Git not found. Detecting package manager..."
PACKAGE_MANAGER=$(detect_package_manager)
if [ "$PACKAGE_MANAGER" == "unknown" ]; then
    echo "No supported package manager found (apt-get, dnf, yum, pacman, brew)."
    echo "Please install manually."
    exit 1
fi

echo "Installing git..."
case "$PACKAGE_MANAGER" in
    "apt-get")
        sudo apt-get install -y git
        ;;
    "dnf")
        sudo dnf install -y git
        ;;
    "yum")
        sudo yum install -y git
        ;;
    "pacman")
        sudo pacman -Sy --noconfirm git
        ;;
    "brew")
        brew install git
        ;;
esac

if command -v git &> /dev/null; then
    echo "Installation complete."
else
    echo "Error: Git installation failed. Please install manually."
    exit 1
fi