#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(
    name="wflib",
    version="00.00.01",
    description="waveform analyzer/simulator",
    author="Leonid Burmistrov",
    packages=find_packages(exclude=['tmp']),
    long_description=description,
)
