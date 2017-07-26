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
  Module[{str, n, data, processLine, 
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
   Reverse[densitydata]];

sliceclumps[slice_] :=
  If[Max[slice] == 1,
   Map[Length, Select[Split[slice], First[#] == 1 &]] // N,
   0];

(* return a list of all clumps along the x-direction in a single file \
*)
getAllClumpSizesHelper[t_] :=
  
  With[{clumplist = Flatten[Map[sliceclumps, t]]},
   N[Select[clumplist, # > 0 &]]];

(* return a list of all clumps in a single file *)
(* returns a \
2-element list, first item is simulation time; second item is a list \
of all clump sizes in the file *)

getAllClumpSizes[file_] :=
  
  Module[{data = readdensity[file], t, allsizes},
   t = Map[If[# < 0.5, 0, 1] &, data, {2}];
   allsizes = Join[getAllClumpSizesHelper[t],
     getAllClumpSizesHelper[Transpose[t]]];
   
   {time, allsizes}];

(* return a histogram of a list of clump sizes in data with log-x \
bins *)
myhistogram[data_] := Block[{bins, y, x},
   bins = Exp[Range[Log[10], Log[2048], 0.2]];
   y = BinCounts[data, {bins}];
   y = y/Total[y] // N;
   x = ListConvolve[{0.5, 0.5}, bins];
   Transpose[{x, y}]];

(* given a simulation direcotry, return a clump size histogram for \
each step in the simulation *)
(* we assume 5 equally-spaced steps in \
simulation time. *)
(* returns a list of 5 histograms, one for each \
step *)
makeClumpHistograms[dir_] :=
  
  Module[{fnames = Drop[FileNames[dir <> "*.vtk"], 1],
    clumpdata, tlim,
    seg1, seg2, seg3, seg4, seg5},
   clumpdata = Map[getAllClumpSizes, fnames];
   tlim = Last[clumpdata[[All, 1]]];
   
   (* simulation is always split into 5 intervals of equal time.  \
(this may not line up with VTK file numbers if there are restarts!) *)

      seg1 = Select[clumpdata, 0.0 tlim < #[[1]] <= 0.2 tlim &];
   seg2 = Select[clumpdata, 0.2 tlim < #[[1]] <= 0.4 tlim &];
   seg3 = Select[clumpdata, 0.4 tlim < #[[1]] <= 0.6 tlim &];
   seg4 = Select[clumpdata, 0.6 tlim < #[[1]] <= 0.8 tlim &];
   seg5 = Select[clumpdata, 0.8 tlim < #[[1]] <= 1.0 tlim &];
   
   seg1 = Flatten[seg1[[All, 2]]];
   seg2 = Flatten[seg2[[All, 2]]];
   seg3 = Flatten[seg3[[All, 2]]];
   seg4 = Flatten[seg4[[All, 2]]];
   seg5 = Flatten[seg5[[All, 2]]];
   
   Map[myhistogram, {seg1, seg2, seg3, seg4, seg5}]];

myfmt[num_] := ToString[FortranForm[num]];

exporthist[dir_, outfile_] :=
  
  Module[{data, x, table, fmtd, header, nres, line},
   data = makeClumpHistograms[dir];
   (* all 5 histograms have the same x-column; 
   split off and put it as the first column *)
   
   x = data[[1, All, 1]];
   data = data[[All, All, 2]];
   
   table = Prepend[data, x];
   table = Transpose[table];
   
   (* format the numbers and add a header *)
   
   fmtd = Map[myfmt, table, {2}];
   header = {"# [1] = size (cells),", " [n+1] = PDF[n]"};
   fmtd = Prepend[fmtd, header];
   
   (* add information about simulation resolution and cstcool sizes \
for later analysis *)
   nres = Max[dim];
   cstcool = nres/(4^(Range[5] - 1));
   line = Prepend[cstcool, nres];
   line = Map[ToString[PaddedForm[#, 10]] &, line];
   header = {"# [1] = simulation resolution,", " [n+1] = cstcool"};
   fmtd = Join[{header, line}, fmtd];
   
   (* write to file *)
   Export[outfile, Append[fmtd, {}]]];

dir = "/Volumes/Neerajstoarage/simulationdata/bigstepfrag/1024/norst/";
exporthist[dir, "test.dat"];
