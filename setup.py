import sys

if sys.version_info[0] > 2:
    raise "Invalid Python interpreter, must be 2.X"

from setuptools import setup, find_packages

setup(
    name="chop",
    version="0.2",
    packages=find_packages(),
    scripts=['src/chop.py', 'src/graph_functions.py',
             'src/eval_mem.py', 'src/colorer.py'],

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
    description="(Un)phased path extraction from population graphs",
    license="MIT",
    keywords="population graph overlap context path extraction haplotype",
    url="https://github.com/tomokveld/CHOP",
)
