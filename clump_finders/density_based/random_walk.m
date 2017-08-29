readVTK[file_String, label_String, type_String] := 
  Module[{str, n, data, dim, processLine, 
    readTime},(*helper functions:*)(*-parse lines in the vtk \
header.*)(*e.g."DIMENSIONS 321 161 161"\[Rule]{321,161,161}*)
   processLine[line_String] := 
    Map[Read[StringToStream[#], Number] &, 
     Drop[StringSplit[line, " "], 1]];
   (*-parse the time out of the comment field*)(*eg "PRIMITIVE vars \
at time= 2.000000e+01, level= 0, domain= 0"\[Rule]20*)
   readTime[line_String] := 
    Block[{regex, timestring}, 
     regex = RegularExpression["time=\\s*([0-9e\\.+\\-]+)"];
     timestring = First[StringCases[line, regex -> "$1"]];
     Read[StringToStream[timestring], Number]];
   (*open the file for reading:*)
   str = OpenRead[file, BinaryFormat -> True];
   (*read the header:*)time = readTime[Find[str, "time"]];
   dim = Map[If[# > 1, # - 1, #] &, 
     processLine[Find[str, "DIMENSION"]]];
   n = Apply[Times, dim];
   origin = processLine[Find[str, "ORIGIN"]];
   spacing = processLine[Find[str, "SPACING"]];
   (*find the data:*)Find[str, label];
   BinaryRead[str, "Character8"];
   If[type == "scalar",(*LOOKUP_TABLE*)
    Block[{}, Read[str, Record]; BinaryRead[str, "Character8"]]];
   (*read the data and close the file:*)
   If[type == "vector", n = 3*n];
   data = 
    BinaryReadList[str, "Real32", n, ByteOrdering -> +1];(*vtk is big-
   endian*)Close[str];
   (*store vectors in a 3D array so that data[[k,j,i]]={vx,xy,vz}*)
   If[type == "vector", data = Partition[data, 3]];
   (*Partition along the x,
   then y axes*)(*Output will always be a 3D array,but Nz and Ny may=
   1*)data = Partition[data, dim[[1]]];
   data = Partition[data, dim[[2]]];
   data];

readdensity[file_] := 
 Module[{densitydata = 
    readVTK[file, "density", "scalar"][[1, All, All]]}, 
  Reverse[densitydata]]

pickpoint[array_] := Module[{point = 0.0, s},
  While[point == 0.0,
   s = {RandomInteger[{1, Length[array]}], 
     RandomInteger[{1, Length[array]}]};
   If[array[[s[[1]], s[[2]]]] == 1.0, point = s]
   
   ];
  s
  ]

edgewalk[pos_, array_] := 
 Module[{angle = RandomReal[{0 , 2 Pi}], out, deltax, deltay, 
   samplepos},
  samplepos = pos;
  deltax = Cos[angle]; deltay = Sin[angle];
  While[array[[Round[samplepos[[1]]], Round[samplepos[[2]]]]] == 1.0,
   samplepos = samplepos + {deltax, deltay};
   ];
  out = Round[samplepos];
  
  If[out[[1]] > Length[array], out[[1]] = Length[array]];
  If[out[[2]] > Length[array], out[[2]] = Length[array]];
  If[out[[1]] < 1, out[[1]] = 1];
  If[out[[2]] < 1, out[[2]] = 1];
  
  out
  
  ]

walklength[initialpos_, array_] := 
 Module[{final}, final = edgewalk[initialpos, array];
  Norm[initialpos - final]]



walklist[array_] := 
 Module[{walks = {}, i = 0, padded}, padded = ArrayPad[array, 1];
  For[i = 0, i < 100, i++,
   walks = AppendTo[walks, walklength[pickpoint[padded], padded]];
   ];
  walks
  ]

mediansize[file_] := 
 Module[{data = readdensity[file], t}, 
  t = Map[If[# < 0.5, 0.0, 1.0] &, data, {2}];
  Median[walklist[t]]
  ]

meansize[file_] := 
 Module[{data = readdensity[file], t}, 
  t = Map[If[# < 0.5, 0.0, 1.0] &, data, {2}];
  Mean[walklist[t]]
  ]


maxsize[file_] := 
 Module[{data = readdensity[file], t}, 
  t = Map[If[# < 0.5, 0.0, 1.0] &, data, {2}];
  Max[walklist[t]]
  ]

historgram[file_] := 
 Module[{data = readdensity[file], t}, 
  t = Map[If[# < 0.5, 0.0, 1.0] &, data, {2}];
  Histogram[N[walklist[t]], {"Log", {0.3}}]
  ]
