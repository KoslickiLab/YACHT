@echo off
setlocal enabledelayedexpansion

rem Set up paths for directories
set SRC_DIR=src\cpp
set BIN_DIR=src\yacht

rem Create bin directory if it doesn't exist
if not exist %BIN_DIR% (
    mkdir %BIN_DIR%
)

rem Compile source files into object files
set OBJ_FILES=
for %%f in (%SRC_DIR%\*.cpp) do (
    g++ -std=c++17 -Wall -O3 -Wsign-compare -c %%f -o %%~nf.o
    if %errorlevel% neq 0 (
        echo Compilation failed for %%f!
        exit /b %errorlevel%
    )
    set OBJ_FILES=!OBJ_FILES! %%~nf.o
)

rem Create yacht_train_core.exe
g++ %OBJ_FILES% -o %BIN_DIR%\yacht_train_core.exe
if %errorlevel% neq 0 (
    echo Linking failed for yacht_train_core.exe!
    exit /b %errorlevel%
)

rem Create yacht_run_compute_similarity.exe
g++ %OBJ_FILES% -o %BIN_DIR%\yacht_run_compute_similarity.exe
if %errorlevel% neq 0 (
    echo Linking failed for yacht_run_compute_similarity.exe!
    exit /b %errorlevel%
)

rem Clean up object files
del *.o

echo Compilation successful. Executables created in %BIN_DIR%.