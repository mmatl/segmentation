"""
Setup of segmentation python codebase
Author: Matthew Matl
"""
from setuptools import setup

requirements = [
    'numpy',
    'matplotlib',
    'networkx',
    'meshpy',
    'visualization'
]


setup(name='segmentation',
      version='0.1.dev0',
      description='Mesh Segmentation Project Code',
      author='Matthew Matl',
      author_email='mmatl@berkeley.edu',
      package_dir = {'': '.'},
      packages=['segmentation'],
      install_requires=requirements,
      test_suite='test'
     )
