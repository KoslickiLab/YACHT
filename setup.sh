#! /bin/bash
# Setup environment for running YACHT


# Check if 'mamba' command is available
if command -v mamba &> /dev/null
then
    echo "Mamba is already installed."
else
    echo "Mamba is not installed. Installing Mamba..."
    
    # Install Mamba from the conda-forge channel
    conda install mamba -c conda-forge
    
    if command -v mamba &> /dev/null
    then
        echo "Mamba has been successfully installed."
    else
        echo "There was a problem installing Mamba. Please check the error messages above."
        exit 1
    fi
fi

# Check if the YACHT environment exists
ENV_NAME="yacht_env"
conda info --envs | grep -w $ENV_NAME

if [ $? -eq 0 ]; then
    echo "The environment '$ENV_NAME' already exists.Please activate it by running 'conda activate $ENV_NAME'.""
else
    mamba env create -f env/yacht_env.yml
    
    if [ $? -eq 0 ]; then
        echo "The environment '$ENV_NAME' has been successfully created and please activate it by running 'conda activate $ENV_NAME'."
    else
        echo "There was a problem creating the environment '$ENV_NAME'. Please check the error messages above."
    fi
fi