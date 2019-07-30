(* ::Package:: *)

(* ::Section:: *)
(*Heat conduction element*)


(* ::Subsection:: *)
(*Initialization*)


Needs["AceGen`"];


If[
	Not@MatchQ[$topology,"L1"|"L2"],
	$topology="T1"
];

SMSInitialize[
	"HeatConvectionD2"<>$topology,
	"Environment"->"AceFEM",
	"Mode"->"Optimal"
];


SMSTemplate[
	"SMSTopology"->$topology,
	"SMSSymmetricTangent"->True,
	"SMSNodeID"->"T",
	"SMSDOFGlobal"->1
];


(* ::Subsection:: *)
(*Element definitions*)


ElementDefinitions[]:=Module[
	{},

	wgp\[DoubleRightTee]SMSIO["Integration weight",Ig];
	\[CapitalXi]={\[Xi],\[Eta],\[Zeta]}\[DoubleRightTee]SMSIO["Integration point",Ig];
	
	\[DoubleStruckCapitalN]h\[DoubleRightTee]Switch[
		$topology,
		"L1",
		{(1-\[Xi])/2, (1+\[Xi])/2},
		"L2",
		{\[Xi] (\[Xi]-1)/2, \[Xi] (1+\[Xi])/2, (1+\[Xi])(1-\[Xi])}
	];
	
	(* Nodal DOF - current *)
	\[DoubleStruckP]eIO\[DoubleRightTee]SMSIO["All DOFs"];
	\[DoubleStruckP]e=Flatten[\[DoubleStruckP]eIO];
	
	(* temperature field *)
	TIO=\[DoubleStruckP]eIO[[All,1]]; 
	T\[DoubleRightTee]\[DoubleStruckCapitalN]h.TIO;
	
	(* Discretization of coordinates *)
	SMSFreeze[\[DoubleStruckCapitalX],Append[\[DoubleStruckCapitalN]h.Take[SMSIO["All coordinates"],Length[\[DoubleStruckCapitalN]h]],\[Zeta]]];
	{X,Y,Z}=\[DoubleStruckCapitalX];
	\[DoubleStruckR]\[Xi]\[DoubleRightTee]SMSD[{X,Y},\[Xi]];
	r\[Xi]n\[DoubleRightTee]SMSSqrt[\[DoubleStruckR]\[Xi].\[DoubleStruckR]\[Xi]];
	fGauss\[DoubleRightTee]wgp*r\[Xi]n;
	
	(* Assumming kg/m/s unit system.*)
	{h\[DoubleStruckG],Tamb\[DoubleStruckG]}\[DoubleRightTee]SMSIO["Domain data",{
		"h"->{"h - convection coefficient",10.},
		"Tamb"->{"Tamb - ambient temperature",20.}
	}];
	(* variables that have to be considered constant for the evaluation of residual *)
	constant={};
	h\[CapitalDelta]T\[DoubleRightTee]SMSFreeze[h\[DoubleStruckG]*(T-Tamb\[DoubleStruckG])];
	AppendTo[constant,h\[CapitalDelta]T];
	
	W\[DoubleRightTee]h\[CapitalDelta]T*T ;
];


(* ::Subsection:: *)
(*Tangent and Residual*)


SMSStandardModule["Tangent and residual"];
	
SMSDo[Ig,1,SMSIO["No. integration points"]];
	
	ElementDefinitions[];

	SMSDo[i,1,SMSNoDOFGlobal];
		Rg\[DoubleRightTee]fGauss*SMSD[W, \[DoubleStruckP]e, i, "Constant"->SMSVariables[constant]];
		SMSExport[SMSResidualSign*Rg,p$$[i],"AddIn"->True];
		SMSDo[j,If[SMSSymmetricTangent,i,1],SMSNoDOFGlobal];		
			Kt=SMSD[Rg, \[DoubleStruckP]e, j];
			SMSExport[Kt,s$$[i,j],"AddIn"->True];
		SMSEndDo[];
	SMSEndDo[];
SMSEndDo[];(* end integration point loop *)


(* ::Subsection:: *)
(*Post-processing*)


SMSStandardModule["Postprocessing"];
(* Nodal post-processing is redundant.
Visualisation of line elements would only interfere with visualisation of 2D elements.  *)


(* ::Subsection:: *)
(*Code generation*)


SMSMainTitle="Heat convection element.";


SMSWrite[];
