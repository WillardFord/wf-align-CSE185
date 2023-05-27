import os
from setuptools import setup, find_packages

# version-keeping code based on pybedtools
curdir = os.path.abspath(os.path.dirname(__file__))
MAJ = 0
MIN = 0
REV = 0
VERSION = '%d.%d.%d' % (MAJ, MIN, REV)
with open(os.path.join(curdir, 'wfalign/version.py'), 'w') as fout:
        fout.write(
            "\n".join(["",
                       "# THIS FILE IS GENERATED FROM SETUP.PY",
                       "version = '{version}'",
                       "__version__ = version"]).format(version=VERSION)
        )

setup(
    name='wf-align',
    version=VERSION,
    description='CSE185 Demo Project', # TODO add longer description
    author='Willard Ford',
    author_email='wwford@ucsd.edu',
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "wf-align=wfalign.wfalign:main"
        ],
    },
)