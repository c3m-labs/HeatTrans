(* ::Package:: *)

(* Paclet Info File *)
(* BuildNumber and Internal values should be inserted during build procedure. *)
Paclet[
	Name -> "HeatTrans",
	Version -> "0.1.0",
	WolframVersion -> "11.+",
	Description -> "Package for non-stationary heat transfer simulation.",
	Creator -> "Matevz Pintar",
	Publisher -> "C3M d.o.o.",
	URL->"http://www.c3m.si",
	Thumbnail->"FrontEnd/icon.png",
	Extensions -> {
		{"Kernel",
			Root -> ".",
			Context ->{"HeatTrans`"}
		},
		{"Documentation",
			Language -> "English",
			MainPage -> "Guides/HeatTrans"
		}
	}
]
