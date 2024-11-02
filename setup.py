from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
import os
import sys
import subprocess

# Import the version number
from src.yacht import __version__

class CustomBuild(build_ext):
    def run(self):
        # Custom build process for compiling the C++ core
        if sys.platform.startswith('win'):
            # Use Windows batch file to compile C++ code
            print("Running build for Windows...")
            subprocess.check_call(['cmd.exe', '/c', 'build_windows.bat'])
        else:
            # Use Unix-based shell script to compile C++ code
            print("Running build for Unix-based system...")
            subprocess.check_call(['bash', 'build_unix.sh'])

        # Run the usual build_ext logic (necessary to continue with setuptools)
        super().run()

setup(
    name='yacht',
    version=__version__,
    include_package_data=True,
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    cmdclass={
        'build_ext': CustomBuild,
    },
    entry_points={
        'console_scripts': [
            'yacht = yacht:main',
        ],
    },
    install_requires=[
        'pandas',
        'scipy',
        'sourmash',
        'loguru',
        'tqdm',
        'biom-format',
        'numpy',
        'setuptools',
        'requests',
        'scikit-learn',
        'pytest',
        'pytest-cov',
        'ruff',
    ],
    python_requires='>=3.6',
    # Add other package metadata here
    author='Koslicki, D., White, S., Ma, C., & Novikov, A.',
    description='YACHT is a mathematically rigorous hypothesis test for the presence or absence of organisms in a metagenomic sample, based on average nucleotide identity (ANI).',
    license='MIT',
    url='https://github.com/KoslickiLab/YACHT'
)
