#!/usr/bin/env wolframscript
(* ::Package:: *)

Print["  Running \"HeatTrans\" tests..."];

(* Modify paths depending from which directory the script is started. *)
PacletDirectoryAdd[ParentDirectory@Directory[]];
Get["HeatTrans`"];

With[{
	tr=TestReport[FileNameJoin[{Directory[],"Tests.wl"}]]
	},
	Exit@If[
		TrueQ[tr["AllTestsSucceeded"]],
		Print["  All tests passed in ",ToString@QuantityMagnitude@tr["TimeElapsed"]," seconds."];
		0, (* exit code if all test pass *)
		Print["  Tests failed! IDs: ",tr["TestsFailedIndices"]];
		1 (* exit code if any test fails *)
	]
]
