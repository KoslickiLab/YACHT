from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
import os
import sys
import subprocess
import shutil

# Import the version number
#from src.yacht import __version__
sys.path.insert(0, os.path.abspath("src"))
from yacht import __version__


# Custom build class to run the C++ compilation step
class CustomBuildExt(build_ext):
    def run(self):
        # Run the custom build process for C++ code
        if sys.platform.startswith('win'):
            # Use the Windows batch file to compile C++ code
            print("Running Windows build script...")
            try:
                subprocess.check_call(['cmd.exe', '/c', 'build_windows.bat'])
            except subprocess.CalledProcessError as e:
                print(f"Error during Windows compilation: {e.output}")
                raise e
        else:
            # Use the Unix-based shell script to compile C++ code
            print("Running Unix-based build script...")
            try:
                subprocess.check_call(['bash', 'build_unix.sh'], stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                print(f"Error during Unix compilation: {e}")
                raise e

        # Move the compiled binary to the correct location for packaging
        compiled_binary = os.path.join('src', 'yacht', 'run_yacht_train_core')
        if os.path.exists(compiled_binary):
            destination = os.path.join(self.build_lib, 'yacht')
            os.makedirs(destination, exist_ok=True)
            shutil.move(compiled_binary, destination)
        else:
            print("Compiled binary not found after build step.")
            raise FileNotFoundError("The executable 'run_yacht_train_core' was not generated successfully.")

        # Run the usual build_ext logic (necessary to continue with setuptools)
        super().run()

class CustomInstall(install):
    def run(self):
        self.run_command('build_ext')
        super().run()

setup(
    name='yacht',
    version=__version__,
    include_package_data=True,
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    cmdclass={
        'build_ext': CustomBuildExt,
        'install': CustomInstall
    },
    entry_points={
        'console_scripts': [
            'yacht = yacht:main',
        ],
    },
    python_requires='>3.6,<3.12',
    # Add other package metadata here
    author='Koslicki, D., White, S., Ma, C., & Novikov, A.',
    description='YACHT is a mathematically rigorous hypothesis test for the presence or absence of organisms in a metagenomic sample, based on average nucleotide identity (ANI).',
    license='MIT',
    url='https://github.com/KoslickiLab/YACHT'
)