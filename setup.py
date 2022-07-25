import setuptools
from setuptools import setup

with open("README.md", "r") as long_info:
    long_description = long_info.read()

setup(
    name='monda',
    version='0.2.0',
    author='',
    author_email='',
    description='A package for retrieval, quality control and analysis of Data from MONOCLE systems',
    url='https://github.com/monocle-h2020/MONDA',
    py_modules=["monda"],
    package_dir={'':'src'},
    packages=setuptools.find_packages(where="src"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=["Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",],
    python_requires=">=3.8",
    install_requires = ['markdown','requests','cartopy','matplotlib','numpy', 'pandas', 'datetime',
                        'argparse'],
)
