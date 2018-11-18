#!/usr/bin/env wolframscript
(* ::Package:: *)

Print["  Running \"HeatTrans\" tests..."];

(* Modify paths depending from which directory the script is started.
By default we assume script is started from git root directory. *)
PacletDirectoryAdd[Directory[]];
Get["HeatTrans`"];

With[{
	tr=TestReport[FileNameJoin[{Directory[],"Tests","Tests.wl"}]]
	},
	If[
		TrueQ@tr["AllTestsSucceeded"]
		,
		Print["\n","  All tests passed in ",ToString@QuantityMagnitude@tr["TimeElapsed"]," seconds."];
		Quit[0] (* exit code if all test pass *)
		,
		Print["\n","  Tests failed! IDs: ",tr["TestsFailedIndices"]];
		Quit[1] (* exit code if any test fails *)
	]
]
