# HeatTrans
Package for non-stationary heat transfer simulation with AceFEM framework

![logo](FrontEnd/Logo.png)


## Installation

The HeatTrans release comes in the form of a `.paclet` file, which contains the entire package and its documentation. 
Download the latest release from the [Github repo's releases page](https://github.com/c3m-labs/HeatTrans/releases). 
To install, run the following command in the Wolfram Language:

    PacletInstall["/full/path/to/HeatTrans-X.Y.Z.paclet"]

This will permanently install the HeatTrans paclet. The Wolfram Language will always use the latest installed version of MeshTools. 
Installed versions can be enumerated using the command:

    PacletFind["HeatTrans"]

And all versions can be uninstalled using the command:

    PacletUninstall["HeatTrans"]

To update the documentation it may be necessary to restart Mathematica.


## Usage

After you have installed the paclet, load it to Mathematica session with `Get`. Then you can, for example, calculate temperature distribution over cross-section of steel bar with 2cm diameter after 60 seconds.

    Get["HeatTrans`"]
    
    result = HeatTransfer[Disk[{0,0},0.01],60,$DefaultMaterial]
    (* InterpolatingFunction[...] *)
    
To access the documentation, open the notebook interface help viewer and search for HeatTrans.


## Contributing and bug reports

Please use the repository [issues](https://github.com/c3m-labs/HeatTrans/issues) page to submit bugs or feature ideas. 

Contributions to this repository are very welcome. Guidelines on how to build paclet file from source code can be found in [Development.md]( Development.md ) file.
