from setuptools import setup

with open("README.md", "r") as long_info:
    long_description = long_info.read()

setup(
    name='monda',
    version='0.0.1',
    description='A package for the retrieval, QC and analysis of Data from MONOCLE systems',
    py_modules=["monda"],
    package_dir={'':'src'},
    packages=['monda', 'monda.WISP', 'monda.WISP.data_access', 'monda.WISP.data_analysis', 'monda.So_Rad'],
    long_description=long_description,
    long_description_content_type="text/markdown"
)