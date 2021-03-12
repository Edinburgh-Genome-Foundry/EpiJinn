from setuptools import setup, find_packages

version = {}
with open("epijinn/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="epijinn",
    version=version["__version__"],
    author="Peter Vegh",
    description="DNA methylation and restriction",
    long_description=open("pypi-readme.rst").read(),
    license="MIT",
    keywords="biology",
    packages=find_packages(exclude="docs"),
    install_requires=["biopython"],
    include_package_data=True,
)
