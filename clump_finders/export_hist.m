(*function to read VTK files into mathematica available \
here:https://github.com/mkmcc/vtk-utils/blob/master/read-vtk.m \
apologies for the terrible formatting...*)(*readVTK:returns a named \
data block from a vtk file as a mathematica array.stores \
`time',`origin',and `spacing' as global variables.example \
usage:ListContourPlot[readVTK["merged/out.0020.vtk","pressure",\
"scalar"][[1,All,All]]] \
ListStreamPlot[readVTK["merged/out.0020.vtk","velocity","vector"][[1,\
All,All,{1,2}]]]*)
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


readpressure[file_] := 
 Module[{pressuredata = 
    readVTK[file, "pressure", "scalar"][[1, All, All]]}, 
  Reverse[pressuredata]]

readtemp[file_] := readpressure[file]/readdensity[file]

readvelocity[file_] := 
 Module[{vel = 
    readVTK[file, "velocity", "vector"][[1, All, All, {1, 2}]]}, 
  Transpose[vel]]

readvx[file_] := 
 Module[{vel = readvelocity[file], vx}, vx = vel[[All, All, 1]];
  Reverse[Transpose[vx]]]

readvy[file_] := 
 Module[{vel = readvelocity[file], vy}, vy = vel[[All, All, 2]];
  Reverse[Transpose[vy]]]

grad[data_] := 
 Module[{out, x = dx[data], y = dy[data]}, 
  out = Transpose[{x, y}, {3, 1, 2}];
  out = out[[2 ;; -2, 2 ;; -2]];
  out(*=ArrayPad[out,1]*)]

baroclinic[fname_] := 
 Module[{gp = 
    grad[readVTK[fname, "pressure", "scalar"][[1, All, All]] // 
      Transpose], 
   gd = grad[
     readVTK[fname, "density", "scalar"][[1, All, All]] // Transpose],
    out}, out = -(gp[[All, All, 1]]*gd[[All, All, 2]] - 
      gp[[All, All, 2]]*gd[[All, All, 1]]);
  Transpose[ArrayPad[out, 1]]]

vorticity[file_] := 
 Module[{v = 
    readVTK[file, "velocity", "vector"][[1, All, All, {1, 2}]] // 
     Transpose, vx, vy, out}, vx = v[[All, All, 1]];
  vy = v[[All, All, 2]];
  out = dx[vy] - dy[vx];
  out = out[[2 ;; -2, 2 ;; -2]];
  Transpose[ArrayPad[out, 1]]]

sliceclumps[slice_] := 
 If[Max[slice] == 1, 
  Map[Length, Select[Split[slice], First[#] == 1 &]] // N, 0]

mscale[file_] := Module[{data, t, tmp, bins, b, bins2},
  data = readdensity[file];
  t = Map[If[# < 0.5, 0, 1] &, data, {2}];
  tmp = Flatten[Map[sliceclumps, t]];
  bins = 10^Range[0, 3, 0.3];
  b = BinCounts[tmp, {bins}];
  bins2 = ListConvolve[{0.5, 0.5}, bins];
  ListLogLogPlot[Transpose[{bins2, b*bins2^2}]]]

meadianmass[file_] := Module[{data, t, tmp, bins, b, bins2},
  data = readdensity[file];
  t = Map[If[# < 0.5, 0, 1] &, data, {2}];
  tmp = Flatten[Map[sliceclumps, t]];
  tmp = Select[tmp, # > 0 &];
  Quantile[tmp, 0.95]
  ]


charsize[file_] := Module[{data = readdensity[file], t, tmp},
  t = Map[If[# < 0.5, 0, 1] &, data, {2}]; 
  tmp = Flatten[Map[sliceclumps, t]];
  tmp = Select[tmp, # > 0 &];
  Mean[N[tmp]]
  ]


maxslice[file_] := Module[{data = readdensity[file], t, tmp},
  t = Map[If[# < 0.5, 0, 1] &, data, {2}]; 
  tmp = Flatten[Map[sliceclumps, t]];
  tmp = Select[tmp, # > 0 &];
  Max[N[tmp]]
  ]

hist2[file_] := Module[{data = readdensity[file], t, list},
  t = Map[If[# < 0.5, 0, 1] &, data, {2}]; 
  list = Flatten[Map[sliceclumps, t]];
  list = Select[list, # > 0 &];
  Histogram[N[list], {"Log", {0.1}}, {"Log", "PDF"}]
  ]

hist2data[file_] := Module[{data = readdensity[file], t, list},
  t = Map[If[# < 0.5, 0, 1] &, data, {2}]; 
  list = Flatten[Map[sliceclumps, t]];
  list = Select[list, # > 0 &];
  N[list]]

myhistogram[data_] :=
 Block[{bins, y, x},
  bins = Exp[Range[Log[10], Log[500], 0.2]];
  y = BinCounts[data, {bins}];
  y = y/Total[y] // N;
  x = ListConvolve[{0.5, 0.5}, bins];
  Transpose[{x, y}]]

avgclumps[files_] := Module[{data, tmp, x, y},
  data = Map[hist2data, files];
  tmp = Map[myhistogram, data];
  x = tmp[[1, All, 1]];
  tmp = Map[#[[All, 2]] &, tmp];
  y = Map[Mean, Transpose[tmp]];
  Transpose[{x, y}]
  ]

clumpcurves[filenames_] := 
 Module[{step0, step1, step2, step3, step4, data0, data1, data2, 
   data3, data4, grouped},
  step0 = filenames[[1 ;; 21]];
  step1 = filenames[[22 ;; 41]];
  step2 = filenames[[42 ;; 61]];
  step3 = filenames[[62 ;; 81]];
  step4 = filenames[[82 ;; 101]];
  grouped = List[step0, step1, step2, step3, step4];
  Map[avgclumps, grouped]
  ]

myfmt[num_] := ToString[PaddedForm[num, {10, 6}]]

exporthist[files_] := Module[{data, x, arranged, fmtd, header},
  data = clumpcurves[files];
  x = data[[1, All, 1]];
  arranged = Transpose[RotateRight[Append[data[[All, All, 2]], x]]];
  fmtd = Map[myfmt, arranged, {2}];
  header = {"# [1] = size (cells)", " [n+1] = PDF[n]"};
  Export["/Users/NeerajAir/clumpfind_repo/data/frag2048.dat",
   Append[Prepend[fmtd, header], {}]]
  ]
