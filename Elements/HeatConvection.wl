(* ::Package:: *)

(* ::Subsection::Closed:: *)
(*Begin package*)


BeginPackage["HeatTrans`Element`Convection`",{"AceFEM`","AceCommon`","AceGen`","AceEnvironment`"}];


makeHeatConvectionElement::usage="makeHeatConvectionElement[model, topology] generates heat convection element for AceFEM.";


(* ::Subsection::Closed:: *)
(*Element code*)


(*Begin["`Private`"];*)


(* ::Subsubsection::Closed:: *)
(*Phases and modules*)


inputOutput//Clear;
inputOutput[]:=(
	XYZ\[RightTee]SMSReal[Table[nd$$[i,"X",j],{i,SMSNoNodes},{j,SMSNoDimensions}]];
	{Xi,Yi}\[DoubleRightTee]Transpose[XYZ];
	
	dof\[RightTee]SMSReal[Table[nd$$[i,"at",j],{i,SMSNoNodes},{j,SMSDOFGlobal[[i]]}]];
	
	{t,\[CapitalDelta]t}\[RightTee]SMSReal[{rdata$$["Time"],rdata$$["TimeIncrement"]}];
	
	{h\[DoubleStruckG],Tamb\[DoubleStruckG]}\[RightTee]SMSReal[Table[es$$["Data",i],{i,Length[SMSDomainDataNames]}]];
)


discretization//Clear;
discretization[model_String,topology_String]:=(
	{\[Xi],\[Eta],\[Zeta],wGauss}\[RightTee]Array[SMSReal[es$$["IntPoints",#1,IpIndex]]&,4];
	Ni\[DoubleRightTee]Switch[
		topology,
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

	fGauss\[DoubleRightTee]Switch[model,
		"D2",
		wGauss*SMSSqrt[r\[Xi].r\[Xi]],
		"AX",
		wGauss*2Pi*Y*SMSSqrt[r\[Xi].r\[Xi]]
	];
)


constitutiveEquations//Clear;
constitutiveEquations[task_String]:=(
	
	Ti\[DoubleRightTee]Flatten[dof];
	T\[DoubleRightTee]Ni.Ti;
	
	If[task==="Post",
		NPlotQuantities={{"Temperature",Ti}};
		GPlotQuantities={};
		Return[]
	];
	
	h\[CapitalDelta]T\[DoubleRightTee]SMSFreeze[h\[DoubleStruckG]*(T-Tamb\[DoubleStruckG])];
	constant={h\[CapitalDelta]T};
	
	\[CapitalPi]\[DoubleRightTee]h\[CapitalDelta]T*T ;
);


(* ::Subsubsection::Closed:: *)
(*Main function*)


makeHeatConvectionElement[model_String,topology_String]:=Block[
	{i,j},
	(* TODO: Add checks which model types and topologies are accepted. *)
	SMSInitialize[
		"HeatConvection"<>model<>topology,
		"Environment"->"AceFEM",
		"Mode"->"Optimal"
	];
	SMSTemplate[
		"SMSTopology"->topology,
		"SMSSymmetricTangent"->True,
		"SMSNodeID"->"T",
		"SMSDOFGlobal"->1,
		"SMSDomainDataNames"->{"h - convective coefficient","Tamb - ambient temperature"},
		"SMSDefaultData"->{10., 25.}
	];
	(* ============================================================= *)
	SMSStandardModule["Tangent and residual"];
	
	NoIp\[RightTee]SMSInteger[es$$["id","NoIntPoints"]];
	
	SMSDo[IpIndex,1,NoIp];
		
		inputOutput[];
		discretization[model,topology];
		constitutiveEquations["\[CapitalPi]"];
		
		SMSDo[i,1,Length[Ti]];
			Rg\[DoubleRightTee]fGauss*SMSD[\[CapitalPi], Ti, i, "Constant"->constant];
			SMSExport[SMSResidualSign*Rg,p$$[i],"AddIn"->True];
			SMSDo[j,If[SMSSymmetricTangent,i,1],Length[Ti]];		
				Kt=SMSD[Rg, Ti, j];
				SMSExport[Kt,s$$[i,j],"AddIn"->True];
			SMSEndDo[];
		SMSEndDo[];
	SMSEndDo[];
	
	(* ============================================================= *)
	SMSStandardModule["Postprocessing"];
	
	NoIp\[RightTee]SMSInteger[es$$["id","NoIntPoints"]];
	
	GPlotQuantities={};
	NPlotQuantities={};
	
	SMSDo[IpIndex,1,NoIp];
		
		inputOutput[];
		discretization[model,topology];
		constitutiveEquations["Post"];
		
		If[
			GPlotQuantities=!={},
			SMSGPostNames=GPlotQuantities[[All,1]];
			SMSExport[GPlotQuantities[[All,2]],gpost$$[IpIndex,#1]&];
		];
	SMSEndDo[];
	
	If[
		NPlotQuantities=!={},
		SMSNPostNames=NPlotQuantities[[All,1]];
		SMSExport[NPlotQuantities[[All,2]],npost$$[#2,#1]&];
	];
	
	SMSMainTitle="Heat convection element.";
	SMSSubTitle=model/.{"AX"->"Axisymmetric model","D2"->"2D continuum model"};
	
	SMSWrite[];
]


(* ::Subsection::Closed:: *)
(*End package*)


(*End[];*)


EndPackage[];
