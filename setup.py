"""Setup module for pyllelic.

This setup file assumes the user is a data scientist using
conda to install dependencies.
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# import subprocess
import pathlib
import pyllelic.__version__ as version

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

# Setup via pip
setup(
    name="pyllelic",
    version=version.__version__,
    author=version.__author__,
    author_email=version.__author_email__,
    description=version.__description__,
    license=version.__license__,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=version.__url__,
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    keywords="genomics, methylation, DNA sequencing",
    python_requires=">=3.6",
    install_requires=[
        "pandas",
        "numpy",
        "scipy",
        "plotly",
        "notebook",
        "xlsxwriter",
        "xlrd",
        "openpyxl",
        "tqdm",
        "pysam",
        "scikit-bio",
        "biopython",
        "ipywidgets",
        "jupyter_contrib_nbextensions",
    ],
)
