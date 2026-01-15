from setuptools import setup, find_packages

setup(
    name='TExTra',
    version='1.1.0',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'TExTra=TExTra.cli:main',
        ],
    },
)

