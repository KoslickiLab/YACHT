from setuptools import setup, find_packages

setup(
    name='yacht',
    version='1.2.3',
    include_package_data=True,
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'yacht = yacht:main',
        ],
    },
    python_requires='>=3.6',
    # Include your pip-installable dependencies here
    install_requires=[
        'sourmash>=4.8.3,<5',
        'scipy',
        'numpy',
        'pandas',
        'scikit-learn',
        'codecov',
        'pytest',
        'pytest-cov',
        'loguru',
        'requests',
        'maturin>=1,<2',
        'tqdm',
        'biom-format',
        'pytaxonkit',
        'openpyxl',
        'ruff',
        'sourmash_plugin_branchwater'
    ],
    # Add other package metadata here
    author='Koslicki, D., White, S., Ma, C., & Novikov, A.',
    description='YACHT is a mathematically rigorous hypothesis test for the presence or absence of organisms in a metagenomic sample, based on average nucleotide identity (ANI).',
    license='MIT',
    url='https://github.com/KoslickiLab/YACHT'
)
