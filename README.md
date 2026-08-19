# HOOMD-blue DEP component

`hoomd-DEP` provides external potentials and forces to simulate dielectrophoresis (DEP) using
either the HPMC or MD integrators within [**HOOMD-blue**](https://glotzerlab.engin.umich.edu/hoomd-blue/).
It includes template C++ and Python modules, an example unit test, CMake scripts to build the component,
and GitHub Actions workflows.

## Building the component

To build hoomd-DEP:

1. Build and install [**HOOMD-blue**](https://hoomd-blue.readthedocs.io/en/latest/building.html) from source.
2. Obtain the source code source.
    ```
    [~]$ git clone https://github.com/johnemmanuelbond/hoomd-DEP.git
    ```
3. Configure.
    ```
    [~/hoomd-DEP]$ cmake -B build -S . -GNinja -D ENABLE_MPI=ON
    ```
4. Build the component.
    ```
    [~/hoomd-DEP/build]$ ninja
    ```
5. Install the component to your current hoomd installation.
    ```
    [~/hoomd-DEP/build]$ ninja install
    ```
    
    Alternatively, you can link the component to your current hoomd build
    ```
    [~] ln -s ~/hoomd-DEP/build/src ~/hoomd-v7/build/hoomd/dep
    ```

Once installed, the template is available for import via:
```
import hoomd.dep
```

## TODO List

hoomd-DEP is still a work in progress. Before this package is ready the devs need to finish:

1. determine if the units module should be here or user-side
2. write and build docs for hosting on readthedocs
3. figure out which unit tests to include
4. deterine whether github actions are worth including for this repo

<!-- ## Creating a new component

To create a new component:

1. Fork [hoomd-component-template](https://github.com/glotzerlab/hoomd-component-template/).
2. Address all **TODO** comments (including those in `.github/`)
3. Add C++ and Python files to `src/`.
4. Add unit tests in `src/pytest`.
5. Format and check code style with [prek](https://prek.j178.dev/). -->

<!-- ## Using the provided GitHub Actions configuration

When you push your changes to GitHub, the [unit test workflow](.github/workflows/unit-test.yaml)
compile your code on the CPU (with and without MPI) and on the GPU (with and without MPI). The
workflow also executes the unit tests on the CPU. You should run GPU unit tests locally, as GitHub
does not provide free GPU runners for GitHub Actions. As a one time step, you need to navigate to
the "Actions" tab of your repository and confirm that GitHub should execute actions for your fork.

When you push a new tag, the [release workflow](.github/workflows/release.yaml) will create a
new GitHub release with automatically generated release notes.

## Maintaining your component

The HOOMD-blue developers will periodically update
[hoomd-component-template](https://github.com/glotzerlab/hoomd-component-template/), including
updates to the GitHub Actions workflow, pre-commit configuration, and CMake scripts. Merge these
changes into your fork to support the latest version of HOOMD-blue.

## Documenting and releasing your component

TODO: Document your component in `README.md` (this file) and remove documentation relevant to the
template.

When appropriate:

* Add [Sphinx](https://www.sphinx-doc.org) documentation and publish it on
[readthedocs](https://www.readthedocs.org).
* Add a [conda-forge](https://conda-forge.org/) package.
* Announce your component on the [HOOMD-blue discussion board](https://github.com/glotzerlab/hoomd-blue/discussions). -->
