# Settings of the database build scripts (override them in the environment)
: "${YACHT_DB_WORK:=$PWD/dbbuild}"   # downloads, sketches, trained and packaged databases (several TB)
: "${YACHT_ENV:=${CONDA_PREFIX:?activate the YACHT conda environment (env/yacht_env.yml + YACHT) or set YACHT_ENV}}"
SCRIPTS=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export YACHT_DB_WORK YACHT_ENV SCRIPTS
