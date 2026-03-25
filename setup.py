from setuptools import find_packages, setup

setup(
    name="OrganelleCleaner",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "networkx",
    ],
    entry_points={
        "console_scripts": [
            "organelle-cleaner=organelle_cleaner.cli:main",
        ],
    },
)
