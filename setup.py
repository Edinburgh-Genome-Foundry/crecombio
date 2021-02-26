from setuptools import setup, find_packages

version = {}
with open("crecombio/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="crecombio",
    version=version["__version__"],
    author="Peter Vegh",
    description="Site-specific DNA recombination simulator",
    long_description=open("pypi-readme.rst").read(),
    license="MIT",
    url="https://github.com/Edinburgh-Genome-Foundry/crecombio",
    keywords="biology",
    packages=find_packages(exclude="docs"),
    install_requires=["biopython"],
)
