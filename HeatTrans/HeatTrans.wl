(* ::Package:: *)

(* ::Section::Closed:: *)
(*Begin package*)


(* ::Subsection::Closed:: *)
(*Header comments*)


(* :Title: HeatTrans *)
(* :Context: HeatTrans` *)
(* :Author: Matevz Pintar *)
(* :Summary: Package for non-stationary of heat transfer simulation with AceFEM framework. *)
(* :Copyright: C3M d.o.o., 2019 *)

(* :Summary: *)


(* ::Subsection::Closed:: *)
(*Begin package*)


(* This enables using AceFEM also without Notebook interface.
It has to be called before BeginPackage, because AceFEM doesn't follow standard conventions. *)
If[
	Not@$Notebooks,
	Get["AceFEM`Remote`"]
];

(* "AceFEM`" context should be called before "AceCommon`", because it loads the AceFEM package,
meanwhile "AceCommon`" doesn't load the package, but puts it on $ContextPath inside the package. *)
BeginPackage["HeatTrans`",{"NDSolve`FEM`","AceFEM`","AceCommon`"}];

(* Clear definitions from package symbols in public and private context. *)
ClearAll["`*","`*`*"];


(* ::Subsection::Closed:: *)
(*Public symbols*)


HeatTransfer;
MakeMesh;
$DefaultMaterial;


(* ::Section::Closed:: *)
(*Code*)


(* Begin private context *)
Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Utilities*)


(* ::Subsubsection::Closed:: *)
(*Paths*)


(* Get the directory of the package, even if the code is evaluated from notebook package editor. 
	RuleDelayed is crucial to avoid evaluation of NotebookDirectory[] when FrontEnd is not available. *)
$packageDirectory=ReplaceAll[
	DirectoryName[$InputFileName],
	("":>NotebookDirectory[])
];


(* This finds element library (.dll or similar) in appropriate directory.
If it doesn't exist it compiles it from source (.C file). This avoids problems
with compilation of elements for other operating systems, assuming users have appropriate 
compiler installed. *)

(* Messages are attached to the main symbol which uses the function to find element library.*)
HeatTransfer::unsupported="$SystemID `1` is not supported in this package.";
HeatTransfer::compile="Element library `1` was not found for this $SystemID.
It will be compiled from source (.C) and saved for reuse.";

getLibrary[name_]:=Module[
	{ext,srcDir,libDir,src,lib},
	ext=Switch[$SystemID,
		"Windows-x86-64",".W64.dll",
		"Linux-x86-64",".L64.dll",
		"MacOSX-x86-64",".M64.dll",
		_,Message[HeatTransfer::unsupported,$SystemID];Abort[]
	];
	
	srcDir=FileNameJoin[{$packageDirectory,"LibraryResources","Source"}];
	libDir=FileNameJoin[{$packageDirectory,"LibraryResources",$SystemID}];
	src=FileNameJoin[{srcDir,name<>".c"}];
	lib=FileNameJoin[{libDir,name<>ext}];
	
	(* Use message attached to General symbol instead of some package symbol. *)
	If[
		Not@FileExistsQ[src],
		Message[General::noopen,src];Abort[]
	];
	
	If[Not@DirectoryQ[libDir],CreateDirectory[libDir]];
	(* If element library doesn't exist, compile it from source (.C) 
	and copy it to appropriate directory. *)
	If[
		Not@FileExistsQ[lib],
		Message[HeatTransfer::compile,FileNameTake@lib];
		SMTMakeDll@src;
		CopyFile[FileNameJoin[{srcDir,FileNameTake@lib}],lib];
		DeleteFile[FileNameJoin[{srcDir,FileNameTake@lib}]]
	];
	lib
];


(* ::Subsection::Closed:: *)
(*HeatTransfer*)


(* ::Subsubsection::Closed:: *)
(*Mesh*)


(* This function improves the quality of 2D mesh. 
It is copied from FEMAddOns package ( https://github.com/WolframResearch/FEMAddOns ) *)

laplacianElementMeshSmoothing[mesh_]:=Module[
	{n, vec, mat, adjacencyMatrix, mass, 
	bndVertices, interiorVertices, stiffness, load, newCoords},

	n = Length[mesh["Coordinates"]];
	vec = mesh["VertexElementConnectivity"];
	mat = Unitize[vec.Transpose[vec]];
	vec = Null;
	adjacencyMatrix = mat - DiagonalMatrix[Diagonal[mat]];
	mass = DiagonalMatrix[SparseArray[Total[adjacencyMatrix, {2}]]];
	stiffness = N[mass - adjacencyMatrix];
	adjacencyMatrix = Null;
	mass = Null;

	bndVertices =  Flatten[Join @@ ElementIncidents[mesh["PointElements"]]];
	interiorVertices = Complement[Range[1, n], bndVertices];

	stiffness[[bndVertices]] = IdentityMatrix[n, SparseArray][[bndVertices]];

	load = ConstantArray[0., {n, mesh["EmbeddingDimension"]}];
	load[[bndVertices]] = mesh["Coordinates"][[bndVertices]];

	newCoords = LinearSolve[stiffness, load];

	ToElementMesh[
		"Coordinates" -> newCoords, 
		"MeshElements" -> mesh["MeshElements"], 
		"BoundaryElements" -> mesh["BoundaryElements"], 
		"PointElements" -> mesh["PointElements"], 
		"CheckIntersections" -> False, 
		"DeleteDuplicateCoordinates" -> False,
		"RegionHoles" -> mesh["RegionHoles"]
	]
];


MakeMesh::usage="MakeMesh[reg,ord] makes ElementMesh from 2D region reg with \"MeshOrder\"->ord.";
MakeMesh::badreg="Region should be a bounded region in 2D.";

MakeMesh//SyntaxInformation={"ArgumentsPattern"->{_,_}};

MakeMesh[region_,order:(1|2)]:=Module[
	{mesh,maxBound},
	
	If[
		Not[BoundedRegionQ[region]&&RegionDimension[region]==2],
		Message[MakeMesh::badreg];Return[$Failed,Module]
	];
	
	maxBound=Max[Differences/@RegionBounds[region]];
	mesh=ToElementMesh[
		region,
		"MeshOrder"->1,
		"MeshElementType"->TriangleElement,
		MaxCellMeasure->{"Length"->maxBound/10},
		AccuracyGoal->2
	];
	
	MeshOrderAlteration[
		laplacianElementMeshSmoothing[mesh],
		order
	]
];


(* ::Subsubsection::Closed:: *)
(*Setup*)


$DefaultMaterial::usage="$DefaultMaterial represents a set of default material properties.";

(* This is a utility symbol, useful for demonstration purpose. Units system is (kg,m,s). *)
$DefaultMaterial=<|"Conductivity"->55,"Density"->7800,"SpecificHeat"->470.|>;


solidDomainData[parameters_Association,other_:{}]:=Join[{
	"kt0 *"->parameters["Conductivity"],
	"rho0 *"->parameters["Density"],
	"cp0 *"->parameters["SpecificHeat"]
	},
	other
];


surfaceDomainData[parameters_Association,other_:{}]:=Join[{
	"h *"->parameters["ConvectionCoefficient"],
	"Tamb *"->parameters["AmbientTemperature"]
	},
	other
];


setupAceFEM//Options={"SaveResults"->False};

(* All parameters are set in this function, analysis function contains just the time stepping loop. *)
setupAceFEM[mesh_,parameters_,opts:OptionsPattern[]]:=Module[
	{order,type,solidTopology,surfaceTopology,solidElm,surfaceElm,resultsFile},
	
	(* Name of results file can be some given string or default name.*)
	resultsFile=ReplaceAll[
		OptionValue["SaveResults"],
		{Automatic|True->"HeatTransfer",x_/;Not@StringQ[x]->False}
	];
	
	order=mesh["MeshOrder"];
	type=Head@First@mesh["MeshElements"];
	solidTopology={type,order}/.{
		{TriangleElement,1}->"T1",
		{TriangleElement,2}->"T2",
		{QuadElement,1}->"Q1",
		{QuadElement,2}->"Q2S"
	};
	solidElm="HeatConductionD2"<>solidTopology;
	surfaceTopology=order/.{1->"L1",2->"L2"};
	surfaceElm="HeatConvectionD2"<>surfaceTopology;
	
	SMTInputData[];
	SMTAddDomain[{
		{"Solid",solidElm,solidDomainData[parameters],"Source"->getLibrary[solidElm]},
		{"Surface",surfaceElm,surfaceDomainData[parameters],"Source"->getLibrary[surfaceElm]}
	}];
	SMTAddMesh[
		mesh,
		{solidTopology->"Solid",surfaceTopology->"Surface"},
		"BoundaryElements"->True
	];
	SMTAddInitialBoundary["T",1->parameters["InitialTemperature"],"Type"->"InitialCondition"];
	SMTAnalysis["DumpInputTo"->resultsFile]
];


(* ::Subsubsection::Closed:: *)
(*Analysis*)


analysisAceFEM[mesh_,time_,parameters_,opts:OptionsPattern[]]:=Module[
	{method,subOpts,timeSteps,data,t0,\[CapitalDelta]tMin,\[CapitalDelta]tMax,step,reaped},
	
	method=OptionValue[HeatTransfer,{opts},Method];
	subOpts=If[ListQ@method,Rest@method,{}];
	
	setupAceFEM[mesh,parameters,FilterRules[Join[subOpts,{opts}],Options@setupAceFEM]];
	
	t0=(OptionValue[HeatTransfer,{opts},StartingStepSize]/.{Automatic->time/1000.});
	\[CapitalDelta]tMax=(OptionValue[HeatTransfer,{opts},MaxStepSize]/.{Automatic->time/10.});
	\[CapitalDelta]tMin=\[CapitalDelta]tMax/1000.;
	
	timeSteps={SMTRData["Time"]};
	data={SMTPostData["Temperature"]};
	
	SMTNextStep["\[CapitalDelta]t"->t0];
	reaped=Last@Last@Reap@While[
		While[
			(step=SMTConvergence[10^-8,10,{"Adaptive Time",7,\[CapitalDelta]tMin,\[CapitalDelta]tMax,time}]),
			SMTNewtonIteration[]
		];
		If[step[[4]]==="MinBound",SMTStatusReport["Analyze"];SMTStepBack[]];
		If[
			Not@step[[1]],
			(* Save AceFEM result files if appropriate option is given at setup.*)
			If[
				MatchQ[SMTSession[[8]],1|2],
				SMTPut[SMTIData["Step"],SMTRData["Time"]]
			];
			Sow@{SMTRData["Time"],SMTPostData["Temperature"]}
		];
		step[[3]],
		If[step[[1]],SMTStepBack[]];
		SMTNextStep["\[CapitalDelta]t"->step[[2]]]
	];
	
	timeSteps=Join[timeSteps,reaped[[All,1]]];
	data=Join[data,reaped[[All,2]]];
	(* Return InterpolatingFunction of temperature field, like NDSolve would.*)
	ElementMeshInterpolation[
		{timeSteps,mesh},
		Transpose[{data},{2,1,3}],
		InterpolationOrder->All,
		"ExtrapolationHandler"->{Function[Indeterminate],"WarningMessage"->False}
	]
];


analysisNDSolve[mesh_,time_,parameters_,opts:OptionsPattern[]]:=Module[
	{iniTemp,ambTemp,conCoeff,rho,cp,k,t0,\[CapitalDelta]tMax},
	(* This is repetitive, but connection with od strings with symbols is clear. *)
	rho=parameters["Density"];
	cp=parameters["SpecificHeat"];
	k=parameters["Conductivity"];
	iniTemp=parameters["InitialTemperature"];
	ambTemp=parameters["AmbientTemperature"];
	conCoeff=parameters["ConvectionCoefficient"];

	t0=(OptionValue[HeatTransfer,{opts},StartingStepSize]/.{Automatic->time/1000.});
	\[CapitalDelta]tMax=(OptionValue[HeatTransfer,{opts},MaxStepSize]/.{Automatic->time/10.});

	NDSolveValue[{
		rho*cp*D[u[t,x,y],t]-k*Laplacian[u[t,x,y],{x,y}]==NeumannValue[conCoeff*(ambTemp-u[t,x,y]),True],
		u[0,x,y]==iniTemp
		},
		u,
		{t,0,time},
		{x,y}\[Element]mesh,
		MaxStepFraction->1.,
		StartingStepSize->t0,
		MaxStepSize->\[CapitalDelta]tMax
	]
];


(* ::Subsubsection::Closed:: *)
(*Main*)


HeatTransfer::usage="HeatTransfer[reg, time, material] simulates heat transfer on 2D Region reg.
HeatTransfer[mesh, time, material] accepts ElementMesh mesh as description of geometry.";
HeatTransfer::eltyp="\"MeshElements\" of given ElementMesh should all either TriangleElement or QuadElement.";
HeatTransfer::bdmtd="Value of option Method->`1` should be \"AceFEM\", \"NDSolve\" or Automatic.";

HeatTransfer//Options={
	"InitialTemperature"->100.,
	"AmbientTemperature"->0.,
	"ConvectionCoefficient"->20.,
	"MeshOrder"->1,
	Method->Automatic,
	StartingStepSize->Automatic,
	MaxStepSize->Automatic
};

HeatTransfer//SyntaxInformation={"ArgumentsPattern"->{_,_,_,OptionsPattern[HeatTransfer]}};

HeatTransfer[region_?RegionQ,time_,material_,opts:OptionsPattern[]]:=Module[
	{mesh,order},
	
	order=OptionValue["MeshOrder"]/.Automatic->1;
	
	If[
		Not[BoundedRegionQ[region]&&RegionDimension[region]==2],
		Message[HeatTransfer::badreg];Return[$Failed,Module]
	];
	
	mesh=MakeMesh[region,order];
	If[
		MatchQ[mesh,_ElementMesh],
		HeatTransfer[mesh,time,material,opts],
		$Failed
	]
];

HeatTransfer[mesh_ElementMesh,time_,material_,opts:OptionsPattern[]]:=Module[
	{bcData,parameters,method},
	
	(* Check for non-mixed mesh and correct embedding dimension .*)
	If[
		Not@MatchQ[Head/@mesh["MeshElements"],{TriangleElement}|{QuadElement}],
		Message[HeatTransfer::eltyp];Return[$Failed,Module]
	];
	
	bcData=AssociationThread[
		{"InitialTemperature","AmbientTemperature","ConvectionCoefficient"},
		{OptionValue["InitialTemperature"],OptionValue["AmbientTemperature"],Clip[OptionValue["ConvectionCoefficient"],{0.,Infinity}]}
	];
	parameters=Merge[{material,bcData},First];
	
	(* Method should converted to one string, even if given with suboptions. *)
	method=OptionValue[Method]/.{{m_,___}:>m}/.{Automatic->"AceFEM"};
	Switch[method,
		"AceFEM",
		analysisAceFEM[mesh,time,parameters,opts],
		"NDSolve",
		analysisNDSolve[mesh,time,parameters,opts],
		_,
		Message[HeatTransfer::bdmtd,Style[method,ShowStringCharacters->True]];Return[$Failed,Module]
	]
];


(* ::Subsection:: *)
(*Postprocessing*)


(* ::Section::Closed:: *)
(*End package*)


End[]; (* "`Private`" *)


EndPackage[];


HeatTransfer::version="Recommended AceFEM package version is at least `1`.";

With[
	{ver=6.912},
	If[
		TrueQ@(ver>SMCSession[[16]]),
		Message[HeatTransfer::version,ToString@ver]
	]
];
