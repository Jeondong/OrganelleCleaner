from setuptools import setup

setup(
    name="organelle-cleaner",
    version="0.1.0",
    py_modules=[
        "organelle_cleaner",
        "organelle_scoring",
        "blast_features",
        "parse_gfa",
        "graph_analysis",
        "sequence_features",
    ],
    install_requires=[
        "networkx",
    ],
    entry_points={
        "console_scripts": [
            "organelle-cleaner=organelle_cleaner:main",
        ],
    },
)
