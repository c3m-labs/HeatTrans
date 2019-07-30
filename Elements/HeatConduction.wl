(* ::Package:: *)

(* ::Section:: *)
(*Heat conduction element*)


(* ::Subsection:: *)
(*Initialization*)


Needs["AceGen`"];


If[
	Not@MatchQ[$topology,"T1"|"T2"|"Q1"|"Q2S"],
	$topology="T1"
];

SMSInitialize[
	"HeatConductionD2"<>$topology,
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


ElementDefinitions[]:=Module[{},

	wgp\[DoubleRightTee]SMSIO["Integration weight",Ig];
	\[CapitalXi]={\[Xi],\[Eta],\[Zeta]}\[DoubleRightTee]SMSIO["Integration point",Ig];

	\[DoubleStruckCapitalN]h\[DoubleRightTee]Switch[
		$topology,
		"T1",
		{\[Xi],\[Eta],1-\[Xi]-\[Eta]},
		"T2",
		\[Kappa]=1-\[Xi]-\[Eta];
		{(2 \[Xi]-1) \[Xi],(2 \[Eta]-1) \[Eta],(2 \[Kappa]-1) \[Kappa],4 \[Xi] \[Eta],4 \[Eta] \[Kappa],4 \[Kappa] \[Xi]},
		"Q1",
		1/4 {(1-\[Xi]) (1-\[Eta]),(1+\[Xi]) (1-\[Eta]),(1+\[Xi]) (1+\[Eta]),(1-\[Xi]) (1+\[Eta])},
		"Q2S",
		{
		1/4(1-\[Xi])(1-\[Eta])(-\[Xi]-\[Eta]-1), 1/4(1+\[Xi])(1-\[Eta])(\[Xi]-\[Eta]-1), 1/4(1+\[Xi])(1+\[Eta])(\[Xi]+\[Eta]-1), 1/4(1-\[Xi])(1+\[Eta])(-\[Xi]+\[Eta]-1),
		1/2(1-\[Xi]^2)(1-\[Eta]), 1/2(1+\[Xi])(1-\[Eta]^2), 1/2(1-\[Xi]^2)(1+\[Eta]), 1/2(1-\[Xi])(1-\[Eta]^2)
		}
	];

	(* Nodal coordinates in real space *)
	\[DoubleStruckCapitalX]IO\[DoubleRightTee]SMSIO["All coordinates"];
	(* Gauss point coordinate in real space*)
	\[DoubleStruckCapitalX]\[RightTee]SMSFreeze[Append[\[DoubleStruckCapitalN]h.\[DoubleStruckCapitalX]IO,\[Zeta]]];
	{X,Y,Z}=\[DoubleStruckCapitalX];

	{t,\[CapitalDelta]t}\[DoubleRightTee]SMSIO[{"Time","TimeIncrement"}];

	(* Nodal DOF - current *)
	\[DoubleStruckP]eIO\[DoubleRightTee]SMSIO["All DOFs"];
	\[DoubleStruckP]e=Flatten[\[DoubleStruckP]eIO];
	(* Nodal DOF - previous *)
	\[DoubleStruckP]enIO\[DoubleRightTee]SMSIO["All DOFs n"];
	\[DoubleStruckP]en=Flatten[\[DoubleStruckP]enIO];

	(* temperature field *)
	TIO=\[DoubleStruckP]eIO[[All,1]]; 
	T\[DoubleRightTee]\[DoubleStruckCapitalN]h.TIO;
	TnIO=\[DoubleStruckP]enIO[[All,1]]; 
	Tn\[DoubleRightTee]\[DoubleStruckCapitalN]h.TnIO;

	(* Approximate data for aluminium. Assumming kg/m/s unit system.*)
	{kt0,kt1,kt2}\[DoubleRightTee]SMSIO["Domain data",{
		"kt0"-> {"kt0 - conductivity (constant term)",220.},
		"kt1"-> {"kt1 - conductivity (linear term)",0.},
		"kt2"-> {"kt2 - conductivity (quadratic term)",0.}
	}];
	{rho0,rho1,rho2}\[DoubleRightTee]SMSIO["Domain data",{
		"rho0"-> {"rho0 - density (constant term)",2700.},
		"rho1"-> {"rho1 - density (linear term)",0.},
		"rho2"-> {"rho2 - density (quadratic term)",0.}
	}];
	{cp0,cp1,cp2}\[DoubleRightTee]SMSIO["Domain data",{
		"cp0"-> {"cp0 - specific heat (constant term)",900.},
		"cp1"-> {"cp1 - specific heat (linear term)",0.},
		"cp2"-> {"cp2 - specific heat (quadratic term)",0.}
	}];
	{Q}\[DoubleRightTee]SMSIO["Domain data",{
		"Q"-> {"Q - heat source",0.}
	}];

	kt \[DoubleRightTee] kt0 + kt1*T + kt2*T^2;
	rho \[DoubleRightTee] rho0 + rho1*T + rho2*T^2;
	cp \[DoubleRightTee] cp0 + cp1*T + cp2*T^2;
	
	(* Kinematics *)
	\[DoubleStruckCapitalJ]e\[DoubleRightTee]SMSD[\[DoubleStruckCapitalX],\[CapitalXi]]; 
	Jed\[DoubleRightTee]SMSDet[\[DoubleStruckCapitalJ]e];
	\[Xi]\[Eta]\[Zeta]ToXYZ={\[CapitalXi],\[DoubleStruckCapitalX],SMSSimplify@SMSInverse[\[DoubleStruckCapitalJ]e]};

	fGauss\[DoubleRightTee]Jed*wgp;

	\[DoubleStruckCapitalD]T\[DoubleRightTee]SMSD[T,\[DoubleStruckCapitalX],"Dependency"->\[Xi]\[Eta]\[Zeta]ToXYZ];
	DTDt\[DoubleRightTee](T-Tn)/\[CapitalDelta]t;
	bT\[DoubleRightTee]SMSFreeze[cp*rho*DTDt+Q];
	constant={bT};
	
	W\[DoubleRightTee](kt*\[DoubleStruckCapitalD]T.\[DoubleStruckCapitalD]T)/2+bT*T;
];


(* ::Subsection:: *)
(*Tangent and Residual*)


SMSStandardModule["Tangent and residual"];

SMSDo[Ig,1,SMSIO["No. integration points"]];

	ElementDefinitions[];
	
	SMSDo[i,1,SMSNoDOFGlobal];
		(* Assuming unit thickness. *)
		Rg\[DoubleRightTee]fGauss*SMSD[W,\[DoubleStruckP]e,i,"Constant"->SMSVariables[constant]];
		SMSExport[SMSResidualSign*Rg,p$$[i],"AddIn"->True];
		SMSDo[j,If[SMSSymmetricTangent,i,1],SMSNoDOFGlobal];
			Kg\[DoubleRightTee]SMSD[Rg,\[DoubleStruckP]e,j];
			SMSExport[Kg,s$$[i,j],"AddIn"-> True];
		SMSEndDo[];
	SMSEndDo[];
SMSEndDo[];(* end integration point loop *)


(* ::Subsection:: *)
(*Post-processing*)


SMSStandardModule["Postprocessing"];

(* nodal post-processing *)
{temperature}\[DoubleRightTee]Transpose[SMSIO["All DOFs"]];
SMSIO[
	{"Temperature"->temperature},
	"Export to",
	"Nodal point"
];


(* ::Subsection:: *)
(*Code generation*)


SMSMainTitle="Transient heat conduction element.";


SMSWrite[];
