hoomd-DEP component
===================

.. description.rst

Based out of the `bevan lab`_, `hoomd-DEP`_ provides external potentials and forces to simulate dielectrophoresis (DEP) using either the HPMC or MD integrators within `hoomd-blue`_. It includes template C++ and Python modules, an example unit test, CMake scripts to build the component, and GitHub Actions workflows.

.. _bevan lab: https://bevan.jh.edu/
.. _hoomd-blue: https://hoomd-blue.readthedocs.io/en/latest/
.. _hoomd-DEP: https://hoomd-dep.readthedocs.io/en/latest/

.. install.rst

Prerequisites
*************

| hoomd-dep is easily set up using anaconda. Create a conda environment with the following packages installed.
| \- `cmake`_, `ninja`_, `pybind11`_: to compile the C++ code
| \- `openmpi`_, `eigen`_ `cereal`_: important C++ libraries
| \- `numpy`_: math in python

.. code-block:: bash

   $ conda create -n ENV_NAME
   $ conda activate ENV_NAME
   $ conda install -c conda-forge numpy cmake ninja pybind11 openmpi eigen ceeal

.. _numpy: https://numpy.org/doc/stable/
.. _cmake: https://cmake.org/cmake/help/latest/
.. _ninja: https://ninja-build.org/
.. _pybind11: https://pybind11.readthedocs.io/en/stable/
.. _openmpi: https://www.open-mpi.org/doc/v4.0/
.. _eigen: https://libeigen.gitlab.io/
.. _cereal: https://uscilab.github.io/cereal/

Building the Component
**********************

To build hoomd-DEP:

1. `Build and install hoomd-blue from source <https://hoomd-blue.readthedocs.io/en/latest/building.html>`_.

2. Obtain the source code source.

.. code-block:: bash

   [~]$ git clone https://github.com/johnemmanuelbond/hoomd-DEP.git

3. Configure.

.. code-block:: bash
   
   [~/hoomd-DEP]$ cmake -B build -S . -GNinja -D ENABLE_MPI=ON


4. Build the component.

.. code-block:: bash
   
   [~/hoomd-DEP/build]$ ninja

5. Install the component to your current hoomd installation.

.. code-block:: bash
   
      [~/hoomd-DEP/build]$ ninja install
    

Alternatively, you can link the component to your current hoomd build

.. code-block:: bash
   
   [~]$ ln -s ~/hoomd-DEP/build/src ~/hoomd-v7/build/hoomd/dep

Once installed, the template is available for import via:

.. code-block:: python
   
   import hoomd.dep

.. docs.rst

Documentation
=============

The full documentation is available at https://hoomd-dep.readthedocs.io/en/latest/

.. <!-- ## Creating a new component

.. To create a new component:

.. 1. Fork [hoomd-component-template](https://github.com/glotzerlab/hoomd-component-template/).
.. 2. Address all **TODO** comments (including those in `.github/`)
.. 3. Add C++ and Python files to `src/`.
.. 4. Add unit tests in `src/pytest`.
.. 5. Format and check code style with [prek](https://prek.j178.dev/). -->

.. <!-- ## Using the provided GitHub Actions configuration

.. When you push your changes to GitHub, the [unit test workflow](.github/workflows/unit-test.yaml)
.. compile your code on the CPU (with and without MPI) and on the GPU (with and without MPI). The
.. workflow also executes the unit tests on the CPU. You should run GPU unit tests locally, as GitHub
.. does not provide free GPU runners for GitHub Actions. As a one time step, you need to navigate to
.. the "Actions" tab of your repository and confirm that GitHub should execute actions for your fork.

.. When you push a new tag, the [release workflow](.github/workflows/release.yaml) will create a
.. new GitHub release with automatically generated release notes.

.. ## Maintaining your component

.. The HOOMD-blue developers will periodically update
.. [hoomd-component-template](https://github.com/glotzerlab/hoomd-component-template/), including
.. updates to the GitHub Actions workflow, pre-commit configuration, and CMake scripts. Merge these
.. changes into your fork to support the latest version of HOOMD-blue.

.. ## Documenting and releasing your component

.. When appropriate:

.. * Add a [conda-forge](https://conda-forge.org/) package.
.. * Announce your component on the [HOOMD-blue discussion board](https://github.com/glotzerlab/hoomd-blue/discussions). -->