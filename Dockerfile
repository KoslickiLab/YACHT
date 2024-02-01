
FROM --platform=linux/x86_64 continuumio/miniconda3:latest
# Set the working directory in the container
WORKDIR /usr/src/app

# Copy the current directory contents into the container at /usr/src/app
COPY . /usr/src/app

RUN conda env create -f env/yacht_env.yml 

# Make RUN commands use the new environment:
# This will activate the conda environment before running any subsequent commands
SHELL ["conda", "run", "-n", "yacht_env", "/bin/bash", "-c"]

# Install the Python package in the environment
RUN pip install .

# The Docker container starts with this command, customize this according to your application's start command
# Note: The ENTRYPOINT is not needed as the SHELL directive already ensures the environment is activated
