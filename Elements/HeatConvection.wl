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
	"SMSDOFGlobal"->1,
	"SMSDomainDataNames"->{"h - convective coefficient","Tamb - ambient temperature"},
	"SMSDefaultData"->{10., 25.} 
	(* Assumming kg/m/s unit system.*)
];


(* ::Subsection:: *)
(*Element definitions*)


ElementDefinitions[]:=Module[{},
	XYZ\[RightTee]SMSReal[Table[nd$$[i,"X",j],{i,SMSNoNodes},{j,SMSNoDimensions}]];
	{Xi,Yi}\[DoubleRightTee]Transpose[XYZ];
	
	dof\[RightTee]SMSReal[Table[nd$$[i,"at",j],{i,SMSNoNodes},{j,SMSDOFGlobal[[i]]}]];
	
	{t,\[CapitalDelta]t}\[RightTee]SMSReal[{rdata$$["Time"],rdata$$["TimeIncrement"]}];
	
	{h\[DoubleStruckG],Tamb\[DoubleStruckG]}\[RightTee]SMSReal[Table[es$$["Data",i],{i,Length[SMSDomainDataNames]}]];
	
	{\[Xi],\[Eta],\[Zeta],wGauss}\[RightTee]Array[SMSReal[es$$["IntPoints",#1,IpIndex]]&,4];
	
	Ni\[DoubleRightTee]Switch[
		$topology,
		"L1",
		{(1-\[Xi])/2, (1+\[Xi])/2},
		"L2",
		{\[Xi] (\[Xi]-1)/2, \[Xi] (1+\[Xi])/2, (1+\[Xi])(1-\[Xi])}
	];
	
	X\[RightTee]SMSFreeze[Ni.Xi];
	Y\[RightTee]SMSFreeze[Ni.Yi];
	Z\[RightTee]SMSFreeze[\[Zeta]];
	
	r\[Xi]\[DoubleRightTee]SMSD[{X,Y},\[Xi]];
	\[Xi]\[Eta]\[Zeta]ToXYZ={};

	fGauss\[DoubleRightTee]wGauss*SMSSqrt[r\[Xi].r\[Xi]];
	
	Ti\[DoubleRightTee]Flatten[dof];
	T\[DoubleRightTee]Ni.Ti;

	h\[CapitalDelta]T\[DoubleRightTee]SMSFreeze[h\[DoubleStruckG]*(T-Tamb\[DoubleStruckG])];
	constant={h\[CapitalDelta]T};
	
	\[CapitalPi]\[DoubleRightTee]h\[CapitalDelta]T*T ;
];


(* ::Subsection:: *)
(*Tangent and Residual*)


SMSStandardModule["Tangent and residual"];
	
	NoIp\[RightTee]SMSInteger[es$$["id","NoIntPoints"]];
	
	SMSDo[IpIndex,1,NoIp];
		
		ElementDefinitions[];
		
		SMSDo[i,1,Length[Ti]];
			Rg\[DoubleRightTee]fGauss*SMSD[\[CapitalPi], Ti, i, "Constant"->constant];
			SMSExport[SMSResidualSign*Rg,p$$[i],"AddIn"->True];
			SMSDo[j,If[SMSSymmetricTangent,i,1],Length[Ti]];		
				Kt=SMSD[Rg, Ti, j];
				SMSExport[Kt,s$$[i,j],"AddIn"->True];
			SMSEndDo[];
		SMSEndDo[];
	SMSEndDo[];


(* ::Subsection:: *)
(*Post-processing*)


SMSStandardModule["Postprocessing"];
	
	NoIp\[RightTee]SMSInteger[es$$["id","NoIntPoints"]];
	
	GPlotQuantities={};
	NPlotQuantities={};
	
	SMSDo[IpIndex,1,NoIp];
		
		ElementDefinitions[];
		
		If[
			GPlotQuantities=!={},
			SMSGPostNames=GPlotQuantities[[All,1]];
			SMSExport[GPlotQuantities[[All,2]],gpost$$[IpIndex,#1]&];
		];
	SMSEndDo[];
	
	NPlotQuantities={{"Temperature",Ti}};
	
	If[
		NPlotQuantities=!={},
		SMSNPostNames=NPlotQuantities[[All,1]];
		SMSExport[NPlotQuantities[[All,2]],npost$$[#2,#1]&];
	];


(* ::Subsection:: *)
(*Code generation*)


SMSMainTitle="Heat convection element.";


SMSWrite[];
