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

readdensity[fname_] := 
  readVTK[fname, "density", "scalar"][[1, All, All]];

(*produces a histogram of number of clumps of each sizes (sizes are \
given in units of cstcool)     *)

clumphistogram[file_] := 
 Module[{rawdensity = readdensity[file], d, p1, comp, clumpIDs, sizes,
    physicalsize, deltax, cstcool}, 
  d = Map[Map[If[# > 7.0, 0.0, 1.0] &, #] &, rawdensity];
  p1 = ArrayPlot[d, Frame -> False, PlotRangePadding -> None, 
    ImageSize -> Last[Dimensions[d]]];
  comp = MorphologicalComponents[p1];
  clumpIDs = Tally[Flatten[comp]];
  
  (*clump id's is a list like {{1,40},{2,70},{3,2887},{4,10}}. 
  To be read as clum number 1 is 40 pixels in area. 
  clump #2 is 70 pixels in area. 
  clump number 4 is 10 pixels in area so on... *)
  
  sizes = Sqrt[clumpIDs[[All, 2]]];
  deltax = (origin[[1]]*-2)/Dimensions[rawdensity][[1]];
  cstcool = 0.1^3.5 (* 
  REMEMBER TO CHANGE THIS FOR DIFFERENT SIMULATIONS*);
  
  Histogram[sizes*deltax/cstcool, {"Log", {0.3}}, 
   AxesLabel -> {"size (\!\(\*SubscriptBox[\(c\), \
\(s\)]\)\!\(\*SubscriptBox[\(t\), \(cool\)]\))", "no of clumps"}, 
   ImageSize -> Medium]
  
  ]

masssize[array_, threshold_] := 
 Module[{d, p1, comp, clumpIDs, sizes, pos, totalclumps, 
   clumpmass = {}, clump, sum, i},
  d = array;
  d = Map[Map[If[# > threshold, 0.0, 1.0] &, #] &, array];
  p1 = ArrayPlot[d, Frame -> False, PlotRangePadding -> None];
  comp = MorphologicalComponents[p1];
  clumpIDs = Tally[Flatten[comp]];
  
  (*clump id's is a list like {{1,40},{2,70},{3,2887},{4,10}}. 
  To be read as clum number 1 is 40 pixels in area. 
  clump #2 is 70 pixels in area. 
  clump number 4 is 10 pixels in area so on... *)
  
  sizes = Sqrt[clumpIDs[[All, 2]]];
  totalclumps = Max[comp] + 1;
  
  For[clump = 0, clump < totalclumps, clump++,
   sum = 0;
   pos = Position[comp, clump]; 
   For[i = 1, i < Length[pos] + 1, i++, 
    sum = sum + array[[pos[[i, 1]], pos[[i, 2]]]] ];
   AppendTo[clumpmass, sum]
   ];
  Transpose[{sizes, clumpmass}]
  
  ]


showclumps[file_, threshold_] := 
 Module[{rawdensity = readdensity[file], d, p1, comp},
  d = Map[Map[If[# > threshold, 0.0, 1.0] &, #] &, rawdensity];
  p1 = ArrayPlot[d, Frame -> False, PlotRangePadding -> None, 
    ImageSize -> Last[Dimensions[d]]];
  comp = MorphologicalComponents[p1];
  Colorize[comp]
  ]
