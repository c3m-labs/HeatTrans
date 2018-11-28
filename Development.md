## Testing

It is considered good practice that every (public) function in this package inclues its own set of unit tests. 
A bunch of them is collected in `Tests/Tests.wl` file, using the Mathematica 
testing [framework](https://reference.wolfram.com/language/guide/SystematicTestingAndVerification.html). 
It is strongly reccomended that tests are run periodically during development and especially before every commit. 

#### Integration of test in Git hook

Unit test can be run automatically before every commit via Git 
client-side [hooks](https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks). 
File `pre-commit` should contain call to `Tests/RunTests.wls` script, which exits with value 0 if all tests pass 
and aborts the commit otherwise. Minimal example of `pre-commit` file content is:

    #!/bin/sh
    ./Tests/RunTests.wls

## How to build the package

#### Prerequisites
* Version 11.3 or greater of Mathematica
* [Wolfram Workbench](https://www.wolfram.com/workbench/)

#### Building documentation
First, import HeatTrans in Workbench:

* Select File -> Import...
* Git - Projects from Git (Next)
* Existing local repository (Next)
* Add... (browse to HeatTrans, select, Next)
* Import as general project (Finish)

Importing of the HeatTrans source needs to be done only once.


Next, build the documentation:
  
* In the HeatTrans folder right click on docbuild.xml
* Choose Run As...
* Choose 2 Ant Build...
* Deselect all 
* Select *HeatTrans*
* Run


This will create a folder named _build_, which will contain a folder _HeatTrans_ that contains the build documentation 
of package.


#### Packaging HeatTrans

Open terminal window (command line) in HeatTrans root directory and run file _Build.wls_. 
This will leave you with a HeatTrans-X.Y.Z.paclet file in the _build_ folder.
