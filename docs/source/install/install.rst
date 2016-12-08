Installation Instructions
=========================

Dependencies
~~~~~~~~~~~~
The `segmentation` module's dependencies can automatically be installed with
`pip`. Note that `meshpy` and `segmentation`, which are also AutoLab repos,
should already be installed.


Cloning the Repository
~~~~~~~~~~~~~~~~~~~~~~
You can clone or download our source code from `Github`_. ::

    $ git clone git@github.com:mmatl/segmentation.git

.. _Github: https://github.com/mmatl/segmentation

Installation
~~~~~~~~~~~~
To install `segmentation` in your current Python environment, simply
change directories into the `segmentation` repository and run ::

    $ pip install -e .

or ::

    $ pip install -r requirements.txt

Alternatively, you can run ::

    $ pip install /path/to/segmentation

to install `segmentation` from anywhere.

Building Documentation
~~~~~~~~~~~~~~~~~~~~~~
Building `segmentation`'s documentation requires a few extra dependencies --
specifically, `sphinx`_ and a few plugins.

.. _sphinx: http://www.sphinx-doc.org/en/1.4.8/

To install the dependencies required, simply run ::

    $ pip install -r docs_requirements.txt

Then, go to the `docs` directory and run `make` with the appropriate target.
For example, ::

    $ cd docs/
    $ make html

will generate a set of web pages. Any documentation files
generated in this manner can be found in `docs/build`.

