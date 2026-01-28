from setuptools import setup, find_packages

setup(
    name="PlasmidNL typing",
    version="2.0",
    author="Gijs Teunis",
    author_email="gijs.teunis@rivm.nl",
    description="Pipeline to reconstruct plasmids from paired Illumina short-read data",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "type_plasmids = scripts.main:main",
        ],
    },
)
