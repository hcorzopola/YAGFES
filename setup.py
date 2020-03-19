#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: hcorzopola
"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="yagfes",
    version="0.9.1",
    author="Héctor A. Corzo Pola",
    author_email="hcorzopola@gmail.com",
    description="Yet another Python 3 package for 2D groundwater flow modeling.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hcorzopola/YAGFES",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
