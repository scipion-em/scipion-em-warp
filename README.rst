===========
Warp plugin
===========

This plugin provides a wrapper for some of the programs of the `Warp <https://github.com/dtegunov/warp>`_ software.

.. image:: https://img.shields.io/pypi/v/scipion-em-warp.svg
        :target: https://pypi.python.org/pypi/scipion-em-warp
        :alt: PyPI release

.. image:: https://img.shields.io/pypi/l/scipion-em-warp.svg
        :target: https://pypi.python.org/pypi/scipion-em-warp
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em-warp.svg
        :target: https://pypi.python.org/pypi/scipion-em-warp
        :alt: Supported Python versions

.. image:: https://img.shields.io/sonar/quality_gate/scipion-em_scipion-em-warp?server=https%3A%2F%2Fsonarcloud.io
        :target: https://sonarcloud.io/dashboard?id=scipion-em_scipion-em-warp
        :alt: SonarCloud quality gate

.. image:: https://img.shields.io/pypi/dm/scipion-em-warp
        :target: https://pypi.python.org/pypi/scipion-em-warp
        :alt: Downloads

Installation
-------------

You will need to use 3.0+ version of Scipion to be able to run these protocols. To install the plugin, you have two options:

a) Stable version

.. code-block::

   scipion installp -p scipion-em-warp

b) Developer's version

   * download repository

    .. code-block::

        git clone -b devel https://github.com/scipion-em/scipion-em-warp.git

   * install

    .. code-block::

       scipion installp -p /path/to/scipion-em-warp --devel


Verifying
---------

To check the installation, simply run one of the following Scipion tests:

* scipion tests warp.tests.test_protocols.TestDeconvolve2D
* scipion tests warp.tests.test_protocols_tomo.TestDeconvolve3D

Protocols
----------

* deconvolution 2D
* deconvolution 3D

References
-----------

1. Dimitry Tegunov and Patrick Cramer. Real-time cryo-electron microscopy data preprocessing with Warp. Nature Methods 16(11), 2019, p. 1146-1152.
