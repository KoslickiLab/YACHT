from setuptools import setup, find_packages

setup(
    name='yacht',
    version='1.0',
    author='Koslicki, D., White, S., Ma, C., & Novikov, A.',
    description='YACHT is a mathematically rigorous hypothesis test for the presence or absence of organisms in a metagenomic sample, based on average nucleotide identity (ANI).',
    packages=find_packages(),
    install_requires=[
        'sourmash>=4.8.3,<5',
        'scipy',
        'numpy',
        'pandas',
        'scikit-learn',
        'loguru',
        'tqdm',
        'biom-format',
        'pytaxonkit',
        'openpyxl'
    ],
    entry_points={
        'console_scripts': [
            'yacht = yacht.cli:main',
        ],
    },
)
