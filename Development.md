# Development of HeatTrans package

If you would like to modify HeatTrans package yourself or contribute back to the original repository,
then the following instructions can help you.
First you need to install [Git](https://git-scm.com/) and
[clone](https://help.github.com/articles/cloning-a-repository/) the project
from its GitHub homepage to your local computer.

## Prerequisites

Essential:

* [Mathematica](https://www.wolfram.com/mathematica/) version 11.3 or later
* [AceFEM/AceGen](http://symech.fgg.uni-lj.si/) packages version 6.912 or later

Recommended:

* [Wolfram Workbench](https://www.wolfram.com/workbench/) for building documentation
* [WolframScript](https://www.wolfram.com/wolframscript/) for building the `.paclet` file from command line.
 On most systems it already comes bundled with Mathematica installation.

## Testing code

It is considered good practice that every (public) function in this package includes its own set of unit tests.
A bunch of them is collected in `Tests/Tests.wl` file, using the Mathematica testing
[framework](https://reference.wolfram.com/language/guide/SystematicTestingAndVerification.html).
It is recommended that you run them periodically during development and especially before every commit.
This can be done by calling script file `Tests/RunTests.wls` in command line
(first change directory to project root directory) or by evaluating whole notebook `Tests/RunTests.nb`.

### Integration of tests in Git hook

Unit test can be run automatically before every commit via Git client-side
[hooks](https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks).
File `pre-commit` should contain call to `Tests/RunTests.wls` script,
which exits with value 0 if all tests pass and aborts the commit otherwise.
Minimal example of `pre-commit` file content is:

    #!/bin/sh
    ./Tests/RunTests.wls

## How to build the package

There a 3 phases in building the package from source code. First you need to generate finite element subroutines
which are later used with AceFEM framework. Documentation notebooks have to be processed for proper integration
into Mathematica documentation center and finally the package is assembled in `.paclet` file.

### Generating element subroutines

Finite element subroutines for assembly of residual vector and tangent matrix are
written using AceGen package functions. AceGen code generates .C files,
which are then compiled to libraries for each operating system.
This procedure is happens by evaluating the whole notebook `Elements/GenerateElements.nb`.

### Building documentation

First, import HeatTrans in Wolfram Workbench:

* Select "File" -> "Open project from file system" ...
* Add path to project root directory and click "Finish"

Importing of the HeatTrans source needs to be done only once.
Next, build the documentation:
  
* Select "Window" -> "Show view" -> "Application Tools"
* In "Application Tools" panel choose the name of the project and click "Build" documentation

This will create a folder named _build_, which will contain
folder _HeatTrans_ that contains the build documentation of package.

### Packaging HeatTrans

Open terminal window (command line) in HeatTrans root directory and run file `Build.wls`.
This will leave you with a HeatTrans-X.Y.Z.paclet file in the _build_ folder.
