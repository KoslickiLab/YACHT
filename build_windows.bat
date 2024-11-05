@echo off

REM Set up paths for directories
set SRC_DIR=src\cpp
set BIN_DIR=src\yacht

REM Create bin directory if it doesn't exist
if not exist %BIN_DIR% (
    mkdir %BIN_DIR%
)

REM Compile the main.cpp file using g++ from MinGW or another suitable compiler
g++ -std=c++17 -Wsign-compare -Wall -O3 -o %BIN_DIR%\run_yacht_train_core.exe %SRC_DIR%\main.cpp

REM Check if compilation succeeded
if %errorlevel% neq 0 (
    echo Compilation failed!
    exit /b %errorlevel%
)

echo Compilation successful. Executable created at %BIN_DIR%\run_yacht_train_core.exe