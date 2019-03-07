(* ::Package:: *)

(* ::Subsection:: *)
(*Begin package*)


(* We will keep element subroutine generating functions in a Global` context to avoid
complications with symbol shadowing. 
One reason for this might be non-conventional AceGen context structure. *)

(*BeginPackage["HeatTrans`Element`Conduction`",{"AceFEM`","AceCommon`","AceGen`","AceEnvironment`"}];*)


makeHeatConductionElement::usage="makeHeatConductionElement[model, topology] generates heat conduction element for AceFEM .";


(* ::Subsection:: *)
(*Element code*)


(*Begin["`Private`"];*)


(* ::Subsubsection:: *)
(*Phases and modules*)


inputOutput//ClearAll;

inputOutput[]:=(
	XYZ\[RightTee]SMSReal[Table[nd$$[i,"X",j],{i,SMSNoNodes},{j,SMSNoDimensions}]];
	{Xi,Yi}\[DoubleRightTee]Transpose[XYZ];
	
	dof\[RightTee]SMSReal[Table[nd$$[i,"at",j],{i,SMSNoNodes},{j,SMSDOFGlobal[[i]]}]];
	dofp\[RightTee]SMSReal[Table[nd$$[i,"ap",j],{i,SMSNoNodes},{j,SMSDOFGlobal[[i]]}]];
	
	{t,\[CapitalDelta]t}\[RightTee]SMSReal[{rdata$$["Time"],rdata$$["TimeIncrement"]}];
	
	{kt0,kt1,kt2,rho0,rho1,rho2,cp0,cp1,cp2,Q}\[RightTee]SMSReal[Table[es$$["Data",i],{i,Length[SMSDomainDataNames]}]];
);


discretization//ClearAll;

discretization[model_String,topology_String]:=(
	{\[Xi],\[Eta],\[Zeta],wGauss}\[RightTee]Array[SMSReal[es$$["IntPoints",#1,IpIndex]]&,4];
	Ni\[DoubleRightTee]Switch[
		topology,
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
	
	X\[RightTee]SMSFreeze[Ni.Xi];
	Y\[RightTee]SMSFreeze[Ni.Yi];
	Z\[RightTee]SMSFreeze[\[Zeta]];
	
	Jm\[DoubleRightTee]SMSD[{X,Y},{\[Xi],\[Eta]}];
	Jmi\[DoubleRightTee]SMSSimplify[SMSInverse[Jm]];
	Jd\[DoubleRightTee]Det[Jm];
	\[Xi]\[Eta]\[Zeta]ToXYZ={{\[Xi],\[Eta]},{X,Y},Jmi};

	fGauss\[DoubleRightTee]Switch[model,
		"D2",
		Jd*wGauss,
		"AX",
		2Pi*Y*Jd*wGauss
	];
);


constitutiveEquations//ClearAll;

constitutiveEquations[task_String]:=(
	
	Ti\[DoubleRightTee]Flatten[dof];
	T\[DoubleRightTee]Ni.Ti;
	Tpi\[DoubleRightTee]Flatten[dofp];
	Tp\[DoubleRightTee]Ni.Tpi;
	
	kt \[DoubleRightTee] kt0 + kt1*T + kt2*T^2;
	rho \[DoubleRightTee] rho0 + rho1*T + rho2*T^2;
	cp \[DoubleRightTee] cp0 + cp1*T + cp2*T^2;
	
	If[task==="Post",
		NPlotQuantities={{"Temperature",Ti}};
		GPlotQuantities={};
		(* Probably we do not need to show these quantities. *)
		(*GPlotQuantities={{"Conductivity",kt},{"Density",rho},{"SpecificHeat",cp}};*)
		Return[]
	];
	
	\[CapitalDelta]T\[DoubleRightTee]SMSD[T,{X,Y,Z},"Dependency"->\[Xi]\[Eta]\[Zeta]ToXYZ];
	(* This use of SMSMax avoids division by 0 in the first analysis step if
	we want to record initial state at time==0. Small value insteadof \[CapitalDelta]t is chosen arbitrarily. *)
	\[Rho]\[CapitalDelta]T\[DoubleRightTee]SMSFreeze[1/SMSMax[\[CapitalDelta]t,10^-6]*rho*cp*(T-Tp)];
	constant={\[Rho]\[CapitalDelta]T};
	
	\[CapitalPi]\[DoubleRightTee](kt*\[CapitalDelta]T.\[CapitalDelta]T)/2+\[Rho]\[CapitalDelta]T*T-T*Q;
);


(* ::Subsubsection:: *)
(*Main function*)


makeHeatConductionElement[model_String,topology_String]:=Block[
	{i,j},
	(* TODO: Add checks which model types and topologies are accepted. *)
	SMSInitialize[
		"HeatConduction"<>model<>topology,
		"Environment"->"AceFEM",
		"Mode"->"Optimal"
	];
	SMSTemplate[
		"SMSTopology"->topology,
		"SMSSymmetricTangent"->True,
		"SMSNodeID"->"T",
		"SMSDOFGlobal"->1,
		"SMSDomainDataNames"->{
			"kt0 - conductivity (constant term)","kt1 - conductivity (linear term)","kt2 - conductivity (quadratic term)",
			"rho0 - density (constant term)","rho1 - density (linear term)","rho2 - density (quadratic term)",
			"cp0 - specific heat (constant term)","cp1 - specific heat (linear term)","cp2 - specific heat (quadratic term)",
			"Q - heat source"
		},
		"SMSDefaultData"->{
			220., 0., 0.,
			2700., 0., 0.,
			900., 0., 0.,
			0.
		} (* Assumming kg/m/s unit system.*)
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
	
	SMSMainTitle="Heat conduction (non-stationary) element. \n Temperature dependent material properties.";
	SMSSubTitle=model/.{"AX"->"Axisymmetric model","D2"->"2D continuum model"};
	
	SMSWrite[];
];


(* ::Subsection:: *)
(*End package*)


(*End[];*)


(*EndPackage[];*)
