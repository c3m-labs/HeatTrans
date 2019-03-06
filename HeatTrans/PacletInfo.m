(* ::Package:: *)

(* Paclet Info File *)
(* BuildNumber and Internal values should be inserted during build procedure. *)
Paclet[
	Name -> "HeatTrans",
	Version -> "0.3.1",
	WolframVersion -> "11.1+",
	Description -> "Package for non-stationary heat transfer simulation.",
	Creator -> "Matevz Pintar",
	Publisher -> "C3M d.o.o.",
	URL->"https://github.com/c3m-labs/HeatTrans",
	Thumbnail->"FrontEnd/Icon.png",
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
