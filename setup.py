from setuptools import setup, find_packages

# Import the version number
from yacht import __version__

setup(
    name='YACHT', 
    version=__version__,
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'pandas',
        'loguru',
        'tqdm',
        'biom-format',
        'numpy',
        'scipy'
    ],
    entry_points={
        'console_scripts': [
            'yacht_sketch_genomes = yacht.sketch_genomes:main',
            'yacht_sketch_sample = yacht.sketch_sample:main',
            'yacht_training = yacht.make_training_data_from_sketches:main',
            'yacht_run = yacht.run_YACHT:main',
            'yacht_standardize = yacht.standardize_yacht_output:main',
        ],
    },
    # Add other package metadata here
    author='KoslickiLab',
    author_email='your.email@example.com',
    description='A mathematically characterized hypothesis test for organism presence/absence in a metagenome',
    license='MIT',  
    keywords='bioinformatics microbial genomics',
    url='https://github.com/KoslickiLab/YACHT', 
)

