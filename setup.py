#!/usr/bin/env python
from setuptools import setup, find_packages

INSTALL_REQUIRES = open("requirements.txt").readlines()

setup(
    name="lagtraj",
    version="0.0.0",
    description="Python trajectory code for Lagrangian simulations",
    url="https://github.com/EUREC4A-UK/lagtraj",
    maintainer="Leif Denby",
    maintainer_email="l.c.denby@leeds.ac.uk",
    py_modules=["lagtraj"],
    packages=find_packages(include=['lagtraj']),
    package_data={"": ["*.csv", "*.yml", "*.html"]},
    include_package_data=True,
    install_requires=INSTALL_REQUIRES,
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    zip_safe=False,
)
