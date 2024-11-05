#!/bin/bash

# Determine the platform
OS_NAME="$(uname -s)"

echo "Running build on Unix-based system: $OS_NAME"
if [[ "$OS_NAME" == "Linux" || "$OS_NAME" == "Darwin" ]]; then
  # Unix-based systems (Linux or macOS)
  bash build_unix.sh
elif [[ "$OS_NAME" == "MINGW"* || "$OS_NAME" == "CYGWIN"* || "$OS_NAME" == "MSYS"* ]]; then
  # Windows-like environment detected, use batch file
  cmd.exe /c build_windows.bat
else
  echo "Unsupported platform: $OS_NAME"
  exit 1
fi