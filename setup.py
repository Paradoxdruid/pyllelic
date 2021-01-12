"""Setup module for pyllelic.

This setup file assumes the user is a data scientist using
conda to install dependencies.
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages  # , Command

# import subprocess
import pathlib

# from distutils.command.build import build as _build

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

# Get version info from __version__.py
version = {}
with open(here / "pyllelic" / "__version__.py") as f:
    exec(f.read(), version)

# Install non-python dependencies
# CUSTOM_COMMANDS = [
#     ["conda", "install", "-c", "bioconda", "emboss"],
#     ["conda", "install", "-c", "bioconda", "perl", "perl-cpanminus"],
#     ["cpan", "install", "Statistics::Lite"],
# ]


# class CustomCommands(Command):
#     """A setuptools Command class able to run arbitrary commands."""

#     def initialize_options(self):
#         pass

#     def finalize_options(self):
#         pass

#     def RunCustomCommand(self, command_list):
#         print("Running command: %s" % command_list)
#         p = subprocess.Popen(
#             command_list,
#             stdin=subprocess.PIPE,
#             stdout=subprocess.PIPE,
#             stderr=subprocess.STDOUT,
#         )
#         # Can use communicate(input='y\n'.encode()) if the command run requires
#         # some confirmation.
#         stdout_data, _ = p.communicate()
#         print("Command output: %s" % stdout_data)
#         if p.returncode != 0:
#             raise RuntimeError(
#                 "Command %s failed: exit code: %s" % (command_list, p.returncode)
#             )

#     def run(self):
#         for command in CUSTOM_COMMANDS:
#             self.RunCustomCommand(command)


# # This class handles the pip install mechanism.
# class build(_build):
#     """A build command class that will be invoked during package install.
#     The package built using the current setup.py will be staged and later
#     installed in the worker using `pip install package'. This class will be
#     instantiated during install for this specific scenario and will trigger
#     running the custom commands specified.
#     """

#     sub_commands = _build.sub_commands + [("CustomCommands", None)]


# Setup via pip
setup(
    name="pyllelic",
    version=version["__version__"],
    author=version["__author__"],
    author_email=version["__author_email__"],
    description=version["__description__"],
    license=version["__license__"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=version["__url__"],
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
        "xlsxwriter",
        "xlrd",
        "samtools",
        "pysam",
        "scikit-bio",
    ],
    # cmdclass={
    #     # Command class instantiated and run during pip install scenarios.
    #     "build": build,
    #     "CustomCommands": CustomCommands,
    # },
)
