(* ::Package:: *)

(* ::Section::Closed:: *)
(*Header comments*)


(* :Title: HeatTrans *)
(* :Context: HeatTrans` *)
(* :Author: Matevz Pintar *)
(* :Summary: Package for non-stationary of heat transfer simulation with AceFEM framework. *)
(* :Copyright: C3M d.o.o., 2018 *)

(* :Summary:
	Add text here.
 *)


(* ::Section::Closed:: *)
(*Begin package*)


(* This enables using AceFEM also without Notebook interface. Explanation why this works
 is a bit complicated. *)
If[Not@$Notebooks,Get["AceFEM`Remote`"]];

(*  "AceFEM`" context should be called before "AceCommon`", because it loads the AceFEM package,
menawhile "AceCommon`" doesn't load the package, but puts this Context on $ContextPath inside the package.  *)
BeginPackage["HeatTrans`",{"NDSolve`FEM`","AceFEM`","AceCommon`"}];


(* ::Subsection::Closed:: *)
(*Public symbols*)


HeatTransfer;
MakeMesh;


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


(* Find element library (DLL or equivalent). *)
getLibrary[name_]:=With[{
	dir=FileNameJoin[{$packageDirectory,"LibraryResources",$SystemID}]
	},
	First[
		FileNames[name<>".*",dir],
		Message[getLibrary::noopen,name];Abort[]
	]
]


(* ::Subsection::Closed:: *)
(*HeatTransfer*)


(* ::Subsubsection::Closed:: *)
(*Mesh*)


(* This function is copied from FEMAddOns package ( https://github.com/WolframResearch/FEMAddOns ). *)
laplacianElementMeshSmoothing[mesh_] := 
Block[
	{n, vec, mat, adjacencymatrix2, mass2, laplacian2, typoOpt, 
	bndvertices2, interiorvertices, stiffness, load, newCoords}, 

	n = Length[mesh["Coordinates"]];
	vec = mesh["VertexElementConnectivity"];
	mat = Unitize[vec.Transpose[vec]];
	vec = Null;
	adjacencymatrix2 = mat - DiagonalMatrix[Diagonal[mat]];
	mass2 = DiagonalMatrix[SparseArray[Total[adjacencymatrix2, {2}]]];
	stiffness = N[mass2 - adjacencymatrix2];
	adjacencymatrix2 = Null;
	mass2 = Null;

	bndvertices2 =  Flatten[Join @@ ElementIncidents[mesh["PointElements"]]];
	interiorvertices = Complement[Range[1, n], bndvertices2];

	stiffness[[bndvertices2]] = IdentityMatrix[n, SparseArray][[bndvertices2]];

	load = ConstantArray[0., {n, mesh["EmbeddingDimension"]}];
	load[[bndvertices2]] = mesh["Coordinates"][[bndvertices2]];

	newCoords = LinearSolve[stiffness, load];

	typoOpt = If[$VersionNumber <= 11.3,
			"CheckIncidentsCompletness" -> False,
			"CheckIncidentsCompleteness" -> False
		];

	ToElementMesh["Coordinates" -> newCoords, 
		"MeshElements" -> mesh["MeshElements"], 
		"BoundaryElements" -> mesh["BoundaryElements"], 
		"PointElements" -> mesh["PointElements"], 
		typoOpt,
		"CheckIntersections" -> False, 
		"DeleteDuplicateCoordinates" -> False,
		"RegionHoles" -> mesh["RegionHoles"]
	]
]


MakeMesh//ClearAll

MakeMesh::usage="MakeMesh[reg,ord] makes ElementMesh from 2D region reg with \"MeshoOrder\"->ord.";

MakeMesh[region_,order:(1|2)]:=Module[
	{triMesh,maxBound},
	
	maxBound=Max[Differences/@RegionBounds[region]];
	triMesh=ToElementMesh[
		region,
		"MeshOrder"->1,
		MaxCellMeasure->{"Length"->maxBound/10},
		AccuracyGoal->2
	];
	MeshOrderAlteration[
		laplacianElementMeshSmoothing@SMTTriangularToQuad[triMesh],
		order
	]
]


(* ::Subsubsection::Closed:: *)
(*Setup*)


assembleDomainData//ClearAll

assembleDomainData[material_Association,other_:{}]:=Join[{
	"kt0 *"->material["Conductivity"],
	"rho0 *"->material["Density"],
	"cp0 *"->material["SpecificHeat"]
	},
	other
]


setup//ClearAll

setup//Options={"SaveResultsTo"->False};

(* We assume that profile mesh is always made of QuadElement (Q1 or Q2S topology) *)
setup[mesh_,material_,opts:OptionsPattern[]]:=Module[
	{order,solidElement,surfaceElement,resultsFile},
	(* Name of results file can be some given string or default name.*)
	resultsFile=ReplaceAll[
		OptionValue["SaveResultsTo"],
		{Automatic|True->"HeatTransfer",x_/;Not@StringQ[x]->False}
	];
	order=mesh["MeshOrder"];
	{solidElement,surfaceElement}=order/.{
		1-> {"HeatConductionD2Q1","HeatConvectionD2L1"},
		2-> {"HeatConductionD2Q2S","HeatConvectionD2L2"}
	};
	
	SMTInputData[];
	SMTAddDomain[{
		{"solid",solidElement,assembleDomainData[material],"Source"->getLibrary[solidElement]},
		{"surface",surfaceElement,{"h *"->10.,"Tamb *"->25.},"Source"->getLibrary[surfaceElement]}
	}];
	SMTAddMesh[mesh,{(order/.{1->"Q1",2->"Q2S"})->"solid",(order/.{1->"L1",2->"L2"})->"surface"},"BoundaryElements"->True];
	
	SMTAnalysis["DumpInputTo"->resultsFile]
]


(* ::Subsubsection::Closed:: *)
(*Analysis*)


analysis//ClearAll

analysis//Options={
	"InitialTemperature"->500,"AmbientTemperature"->25,"ConvectionCoefficient"->20.
	};

analysis[mesh_,time_,timeStep_,opts:OptionsPattern[]]:=Module[
	{initialTemp,convCoeff,allTimeSteps,reaped,data},
	
	initialTemp=OptionValue["InitialTemperature"];
	convCoeff=Clip[OptionValue["ConvectionCoefficient"],{0.,Infinity}];
	
	allTimeSteps=Range[0.,time,timeStep];
	
	SMTAddInitialBoundary["T",1->initialTemp,"Type"->"InitialCondition"];
	SMTDomainData["surface","Data","h*"->convCoeff];
	
	reaped=Reap@Do[
		SMTNextStep["t"->t];
		While[
			SMTConvergence[10^-8,15],
			SMTNewtonIteration[];
		];
		(* Save AceFEM result files if appropriate option is given at setup.*)
		If[
			MatchQ[SMTSession[[8]],1|2],
			SMTPut[SMTIData["Step"],SMTRData["Time"]]
		];
		Sow[SMTPostData["Temperature"]]
		,
		{t,allTimeSteps}
	];
	
	data=Transpose[Last@reaped,{2,1,3}];
	(* Return InterpolatingFunction of temperature field, like NDSolve would.*)
	ElementMeshInterpolation[
		{allTimeSteps,mesh},
		data,
		InterpolationOrder->All,
		"ExtrapolationHandler"->{Function[Indeterminate],"WarningMessage"->False}
	]
]


(* ::Subsubsection::Closed:: *)
(*Main*)


HeatTransfer//ClearAll

HeatTransfer::usage="HeatTransfer[reg, time, material] simulates 2D heat transfer.";
HeatTransfer::badreg="Profile cross-section should be bounded region in 2D.";

HeatTransfer//Options=Sort@Join[
	Flatten[Options/@{setup,analysis}],
	{"TimeStep"->Automatic,"MeshOrder"->1}
];

HeatTransfer//SyntaxInformation={"ArgumentsPattern"->{_,_,_,OptionsPattern[HeatTransfer]}};

HeatTransfer[region_,time_,material_,opts:OptionsPattern[]]:=Module[
	{mesh,timeStep,order},
	timeStep=OptionValue["TimeStep"]/.Automatic->1.;
	order=OptionValue["MeshOrder"]/.Automatic->1;
	
	If[
		Not[BoundedRegionQ[region]&&RegionDimension[region]==2],
		Message[HeatTransfer::badreg];Return[$Failed]
	];
	
	mesh=MakeMesh[region,order];
	
	setup[mesh,material,FilterRules[{opts},Options@setup]];
	analysis[mesh,time,timeStep,FilterRules[{opts},Options@analysis]]
]


(* ::Subsection:: *)
(*Postprocessing*)


(* ::Section::Closed:: *)
(*End package*)


End[]; (* "`Private`" *)


EndPackage[];
