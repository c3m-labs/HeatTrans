(* ::Package:: *)

(* ::Subsection:: *)
(*Description*)


(* This file contains various (unit) tests for "HeatTrans" package. *)

(* "HeatTrans" package must be already loaded before running these tests, 
	otherwise testing is aborted. *)
If[
	Not@MemberQ[$Packages,"HeatTrans`"],
	Print["Error: Package is not loaded!"];Abort[];
];


BeginTestSection["Tests"];


$miniTestMesh=ToElementMesh[
	"Coordinates"->{{1.,0.},{1.,1.},{1.5,0.},{1.5,1.},{2.,0.},{2.,1.},{1.,0.5},{1.4,0.4},{2.,0.5}}, 
	"MeshElements"->{QuadElement[{{1,3,8,7},{7,8,4,2},{3,5,9,8},{8,9,6,4}}]}
];


(* ::Subsection:: *)
(*MakeMesh*)


VerificationTest[
	With[{
		mesh=MakeMesh[Disk[],1]
		},
		(* Test should return True *)
		And[
			Head@First@mesh["MeshElements"]==QuadElement,
			mesh["MeshOrder"]==1
		]
	],
	TestID->"MakeMesh_success_order=1"
];


VerificationTest[
	With[{
		mesh=MakeMesh[Disk[],2]
		},
		(* Test should return True *)
		And[
			Head@First@mesh["MeshElements"]==QuadElement,
			mesh["MeshOrder"]==2
		]
	],
	TestID->"MakeMesh_success_order=2"
];


VerificationTest[
	MakeMesh[Disk[],"badValue"],
	MakeMesh[Disk[],"badValue"],
	TestID->"MakeMesh_badMeshOrder"
];


VerificationTest[
	MakeMesh[Ball[],1],
	$Failed,
	{MakeMesh::badreg},
	TestID->"MakeMesh_fail_3D_region"
];


VerificationTest[
	MakeMesh[HalfPlane[{{0,0},{1,0}},{0,1}],1],
	$Failed,
	{MakeMesh::badreg},
	TestID->"MakeMesh_fail_unbounded_region"
];


(* ::Subsection:: *)
(*HeatTransfer*)


VerificationTest[
	HeatTransfer[$miniTestMesh,1,$DefaultMaterial,Method->"AceFEM"]["Domain"]//First,
	{0.,1.},
	TestID->"HeatTransfer_AceFEM"
];


VerificationTest[
	HeatTransfer[$miniTestMesh,1,$DefaultMaterial,Method->"NDSolve"]["Domain"]//First,
	{0.,1.},
	TestID->"HeatTransfer_NDSolve"
];


VerificationTest[
	HeatTransfer[$miniTestMesh,1,$DefaultMaterial,
		Method->"AceFEM",
		MaxStepSize->0.2,
		StartingStepSize->0.2
	]["Coordinates"]//First,
	{0.,0.2,0.4,0.6,0.8,1.},
	TestID->"HeatTransfer_AceFEM-step-size"
];


VerificationTest[
	HeatTransfer[$miniTestMesh,1,$DefaultMaterial,
		Method->"NDSolve",
		MaxStepSize->0.2,
		StartingStepSize->0.2
	]["Coordinates"]//First,
	{0.,0.2,0.4,0.6,0.8,1.},
	TestID->"HeatTransfer_NDSolve-step-size"
];


VerificationTest[
	HeatTransfer[$miniTestMesh,1,$DefaultMaterial,Method->"BadValue"],
	$Failed,
	{HeatTransfer::bdmtd},
	TestID->"HeatTransfer_wrong-Method"
];


VerificationTest[
	HeatTransfer[ToElementMesh@Disk[],1,$DefaultMaterial],
	$Failed,
	{HeatTransfer::quadElms},
	TestID->"HeatTransfer_fail_triangleElements"
];


(* Remove global symbols used only for testing. *)
Remove[$miniTestMesh];


(* ::Subsection:: *)
(*EndTestSection*)


EndTestSection[];
