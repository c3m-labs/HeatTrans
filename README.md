# HeatTrans
[Mathematica](http://www.wolfram.com/mathematica/) package for non-stationary heat transfer simulation with [AceFEM](http://symech.fgg.uni-lj.si/) framework

![logo](FrontEnd/Logo.png)


## Installation

To use "HeatTrans" package you need Mathematica version 11. or later and 
AceFEM package, version 6.912 or later ([trial](http://symech.fgg.uni-lj.si/Download.htm) version). 
Supported operating systems are 64-bit Windows, MacOS and Linux.
The HeatTrans release comes in the form of a `.paclet` file, which contains the entire package and its documentation. 
Download the latest release from the [Github repo's releases page](https://github.com/c3m-labs/HeatTrans/releases). 
To install, run the following command in the Mathematica:

    PacletInstall["/full/path/to/HeatTrans-X.Y.Z.paclet"]

This will permanently install the HeatTrans paclet. The Mathematica will always use the latest installed version of MeshTools. 
Installed versions can be enumerated using the command:

    PacletFind["HeatTrans"]

To update the documentation it may be necessary to restart Mathematica. 
All versions can be uninstalled with `PacletUninstall["HeatTrans"]`.


## Usage

After you have installed the paclet, load it to Mathematica session with `Get`. 
Then you can, for example, calculate temperature distribution over triangular steel bar after 60 seconds.

    Get["HeatTrans`"]
    
    result = HeatTransfer[Triangle[{{0,0},{1,1},{2,0}}/100],60,$DefaultMaterial]
    (* InterpolatingFunction[...] *)
    
    Plot3D[
        result[60, x, y],
        {x, y} \[Element] (result["ElementMesh"]),
        ColorFunction -> "TemperatureMap",
        AspectRatio -> Automatic,
        Boxed -> False,
        PlotRange -> All
    ]

![plot3D](https://i.imgur.com/F0pGWwV.png)
    
To access the documentation, open the notebook interface help viewer and search for HeatTrans.


## Contributing and bug reports

Please use the repository [issues](https://github.com/c3m-labs/HeatTrans/issues) page to submit bugs or feature ideas. 

Contributions to this repository are very welcome. Guidelines on how to build paclet file from source code can be found in [Development.md]( Development.md ) file.
