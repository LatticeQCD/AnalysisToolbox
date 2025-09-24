from setuptools import setup, find_packages


def read_requirements():
    """Read requirements.txt and return a list of dependencies."""
    with open("/home/dclarke/GitHub/AnalysisToolbox/requirements.txt", "r") as fh:
        return fh.read().splitlines()


with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    name="latqcdtools",
    version="1.3.1",
    author="D. A. Clarke",
    author_email="clarke.davida@gmail.com",
    description="A set of Python tools for statistically analyzing correlated data. This includes aspects of lattice QCD applications related to QCD phenomenology.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/LatticeQCD/AnalysisToolbox",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    install_requires=read_requirements(),
)
