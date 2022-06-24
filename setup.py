from setuptools import setup, find_packages
from codecs import open
import os
import re

with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()

def get_version():
    version_file = os.path.join(os.path.abspath(os.path.dirname(__file__)),'wfutil','__version__.py')

    with open(version_file, "r") as f:
        lines = f.read()
        version = re.search(r"^_*version_* = ['\"]([^'\"]*)['\"]", lines, re.M).group(1)
        f.close()

        return version

wfutil_version = get_version()

setup(
    name='wfutil',
    packages=find_packages(),
    version=wfutil_version,
    author='Kristof Rozgonyi',
    author_email='rstofi@gmail.com',
    description='Collection of utility code to work with time and frequency arrays and time-frequency data matrices for creating wtarefall plots',
    #install_requires=requirements
    )

