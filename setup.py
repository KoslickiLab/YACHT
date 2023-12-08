from setuptools import setup, find_packages

# Import the version number
from yacht import __version__

setup(
    name='yacht',
    version=__version__,
    include_package_data=True,
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'yacht = yacht:main',
        ],
    },
    python_requires='>=3.6',
    # Add other package metadata here
    author='Koslicki, D., White, S., Ma, C., & Novikov, A.',
    description='YACHT is a mathematically rigorous hypothesis test for the presence or absence of organisms in a metagenomic sample, based on average nucleotide identity (ANI).',
    license='MIT',
    url='https://github.com/KoslickiLab/YACHT'
)
