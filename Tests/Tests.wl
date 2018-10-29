(* ::Package:: *)

(* ::Subsection::Closed:: *)
(*Description*)


(* ::Text:: *)
(*These are various test for "HeatTrans" paclet.*)


(* "HeatTrans.wl" must be loaded before running these tests, otherwise testing is aborted. *)
If[
	Not@MemberQ[$Packages,"HeatTrans`"],
	Print["Error: Package is not loaded!"];Abort[];
];


(* Currently it is unclear what this line does, it is automatically gnerated during conversion to .wlt *)
BeginTestSection["Tests"]


(* ::Subsection::Closed:: *)
(*MakeMesh*)


VerificationTest[
	With[{mesh=MakeMesh[Disk[],1]},
		(* Test should return True *)
		And[
			Head@First@mesh["MeshElements"]==QuadElement,
			mesh["MeshOrder"]==1
		]
	],
	TestID->"MakeMesh_success_order=1"
]


VerificationTest[
	With[{mesh=MakeMesh[Disk[],2]},
		(* Test should return True *)
		And[
			Head@First@mesh["MeshElements"]==QuadElement,
			mesh["MeshOrder"]==2
		]
	],
	TestID->"MakeMesh_success_order=2"
]


(* ::Subsection::Closed:: *)
(*HeatTransfer*)


VerificationTest[
	HeatTransfer[Disk[],10,$DefaultMaterial],
	_InterpolatingFunction,
	SameTest->MatchQ,
	TestID->"HeatTransfer_success"
]


(* ::Subsection::Closed:: *)
(*EndTestSection*)


EndTestSection[]
