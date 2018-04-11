import sys

if sys.version_info[0] > 2:
    raise "Invalid Python interpreter, must of 2.X"

from setuptools import setup, find_packages

setup(
    name="chop",
    version="0.2",
    packages=find_packages(),
    scripts=['chop.py', 'graph_functions.py', 'eval_mem.py', 'colorer.py'],

    install_requires=['docutils>=0.3', 'networkx==1.9', 'bitarray'],
    package_data={
        '': ['*.txt', '*.rst'],
    },
    entry_points={
        'console_scripts': [
            'chop = chop:main'
        ]
    },

    author="Tom Mokveld",
    author_email="T.O.Mokveld@tudelft.nl",
    description="(Un)phased Path extraction from graphs",
    license="PSF",
    keywords="population graph overlap context path extraction haplotype",
    url="",
)
