from setuptools import setup, find_packages

setup(
    name='blast_python',
    version='0.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'blast_python = blast_python.main:main',
        ],
    },
)