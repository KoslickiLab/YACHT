#! /bin/bash
# Setup environment for running YACHT

# Check if the YACHT environment exists
ENV_NAME="yacht_env"
check=$(conda env list | cut -d" " -f 1 | grep -w $ENV_NAME | wc -l)

if [ $check -eq 1 ]; then
    echo "The environment '$ENV_NAME' already exists. Please activate it by running 'conda activate $ENV_NAME'."
else
    conda env create -f env/yacht_env.yml
    
    if [ $check -eq 1 ]; then
        echo "The environment '$ENV_NAME' has been successfully created and please activate it by running 'conda activate $ENV_NAME'."
    else
        echo "There was a problem creating the environment '$ENV_NAME'. Please check the error messages above."
    fi
fi