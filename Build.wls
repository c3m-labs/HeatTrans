#!/usr/bin/env wolframscript
(* ::Package:: *)

getGitRevision[dir_]:=Module[{description},
	SetDirectory[dir];
	(* For meaning of this see: https://git-scm.com/docs/git-describe *)
	description=First[
		ReadList["!git describe --long --tags --dirty --always",String],
		Return[$Failed]
	];
	ResetDirectory[];
	description
]


getGitCommitCount[dir_]:=Module[{value},
	SetDirectory[dir];
	(* For meaning of this see: https://git-scm.com/docs/git-describe *)
	value=First[
		ReadList["!git rev-list --count HEAD",Number],
		Return[$Failed]
	];
	ResetDirectory[];
	value
]


$name="HeatTrans";
(* Get the distribution directory by using this notebook\[CloseCurlyQuote]s directory. *)
$source=If[$Notebooks,NotebookDirectory[],Directory[]];
(* Construct the target directory from this notebooks base directory. *)
$target=FileNameJoin[{$source,"build",$name}];


Module[
	{directories,original,modified},
	If[
		Not@DirectoryQ[$target],
		CreateDirectory[$target,CreateIntermediateDirectories->True]
	];
	(* Copy some files from the source to the target directory. *)
	Map[
		CopyFile[
			FileNameJoin[{$source,#}],
			FileNameJoin[{$target,#}],
			OverwriteTarget->True
		]&,
		{"HeatTrans.wl","PacletInfo.m"}
	];
	(* Delete old and copy the new source directories. 
	"Documentation" directory is already created by documentation build procedure. *)
	directories={"FronEnd","Kernel","LibraryResources"};
	Map[
		If[
			FileExistsQ[FileNameJoin[{$target,#}]],
			DeleteDirectory[FileNameJoin[{$target,#}],DeleteContents->True]
		]&,
		directories
	];
	Map[
		CopyDirectory[FileNameJoin[{$source,#}],FileNameJoin[{$target,#}]]&,
		directories
	];
	
	original=Import[FileNameJoin[{$source,"PacletInfo.m"}]];
	(* Description of git repository (git describe ...) returns a string that cannot 
	be used as BuildNumber because only integers are accepted. 
	Otherwise PackPaclet returns $Failed and procedure fails. Instead now commit count is used. *)
	modified=With[{
		noCommits=getGitCommitCount[$source],
		revision=getGitRevision[$source]
		},
		Fold[
			Insert[#1,#2,4]&,
			original,
			{Internal->True,BuildNumber->noCommits,"Hash"->revision}
		]
	];
	
	Export[
		FileNameJoin[{$target,"PacletInfo.m"}],
		modified,
		"Comments"->{"Paclet Info File","Created on "<>DateString[]}
	];
]


Module[
	{pacletPath,newPaclet},
	Needs["PacletManager`"];
	(* Change into the build directory. *)	
	SetDirectory[FileNameJoin[{$source,"build"}]];
	(*Create the paclet*)
	pacletPath=PackPaclet[$name];
	ResetDirectory[];
	(*Uninstall a possibly old version. *)
	If[PacletFind[$name]=!={},PacletUninstall[$name]];
	(*Install the paclet*)
	newPaclet=PacletInstall[pacletPath];
	RebuildPacletData[];
	If[
		Length@PacletFind["HeatTrans"]===1,
		Print[0],
		Print[1]
	]
]
