(*read file parameters*)

tlim = 15.0;
amrtimes = {10, 13};
alpha = -1.0;
gm = N[5/3];

f = gm^1.5/((gm - 1)*(2.0 - alpha));

drat = 10.0;

lambda0 = 1.0/f;
lambda0 *= 1/(drat^(2.5 - alpha));
lambda0 *= 4.0;

cstcool[x_] := 
 Module[{droplist, t, sorted, pos}, t = Table[2^-i, {i, 1, 5}];
  t = t/Total[t];
  t = N[Accumulate[t]];
  sorted = Sort[Append[t, x]];
  pos = Position[sorted, x][[1, 1]];
  N[4^-pos]
  ]


tfloor[x_] := (f*lambda0*cstcool[x])^(2.0/(5.0 - 2.0*alpha))
  
  (*start with functions to read athhdf files*)

  rhoplus[file_] :=
 Module[{b, loc, pos, data, rho, part},
  b = Import[file, "Data"];
  
  (*read in logical ordering of blocks (assuming 2D data!)*)
  
  loc = b[[2]][[All, {1, 2}]];
  loc = Map[Reverse, loc];
  
  (*sort blocks*)
  pos = Ordering[loc];
  data = b[[3, All, All, 1, All, All]];
  rho = data[[1, pos, All, All]];
  
  (*partition meshblocks into a 2D array*)
  
  part = Length[Select[loc, First[#] == 0 &]];
  rho = Partition[rho, part];
  
  (*use arrayflatten to compress into one big array*)
  
  rho = ArrayFlatten[rho];
  rho]
presplus[file_] := Module[{b, loc, pos, data, P, part},
  b = Import[file, "Data"];
  
  (*read in logical ordering of blocks (assuming 2D data!)*)
  
  loc = b[[2]][[All, {1, 2}]];
  loc = Map[Reverse, loc];
  
  (*sort blocks*)pos = Ordering[loc];
  data = b[[3, All, All, 1, All, All]];
  P = data[[2, pos, All, All]];
  
  (*partition meshblocks into a 2D array*)
  
  part = Length[Select[loc, First[#] == 0 &]];
  P = Partition[P, part];
  
  (*use arrayflatten to compress into one big array*)
  
  P = ArrayFlatten[P];
  P]
  
tempplus[file_] := presplus[file]/rhoplus[file]
(*threshold set at 2.0*tfloor*)

mybinarize[data_, time_] := Module[{tflr, x},
  x = time/tlim;
  tflr = tfloor[x]; 
  Map[If[# < 2.0*tflr, 1.0, 0.0] &, 
   data, {2}]](*thresholds data based on a temperature floor*)

readhdf5time[file_] := 
 Block[{str = OpenRead[file, BinaryFormat -> True], line}, 
  line = Find[str, "Time"];
  SetStreamPosition[str, StreamPosition[str] - StringLength[line]];
  line = BinaryReadList[str, "Character8", 256] // StringJoin;
  Close[str];
  
  line = StringSplit[line, "Time"][[-1]];
  str = StringToStream[line];
  SetStreamPosition[str, 36];
  BinaryRead[str, "Real64", ByteOrdering -> +1]]
  
  (*functions to analyze a single athdf file*)

   (*return a list of clump sizes along the 1D skewer slice*)
(*slice is \
presumed to contain only 1s and 0s... we return a list of all the \
lengths of contiguous patches of 1s*)

sliceclumps[slice_] :=
  If[Max[slice] == 1,
    Map[Length, Select[Split[slice], First[#] == 1 &]],
    0] // N;
(*return a list of all clumps along the x-direction in a single file*)

getAllClumpSizesHelper[t_] :=
  
  With[{clumplist = Flatten[Map[sliceclumps, t]]},
   N[Select[clumplist, # > 0 &]]];
myhistogram[clumpdata_] :=
 Module[{bins, y, x, clumps, time, dim, dx},
  dim = clumpdata[[1]];
  time = clumpdata[[2]];
  clumps = clumpdata[[3]];
  dx = N[1.0/dim];(*assuming a domain size of 1.0*)
  
  clumps = clumps*dx; (*clump size in physical units*)
  
  bins = Exp[Range[Log[10^-3], Log[1.0], N[1/29*Log[10^3]]]];
  
  (*bin clump data and normalize*)
  y = BinCounts[clumps, {bins}];
  
  (*average bins*)
  x = Exp[ListConvolve[{0.5, 0.5}, Log[bins]]];
  Transpose[{x, x*y*dx}]
  ]
(* function to check whether the thresholded data t is isotropic. *)
(* \
isotropic if integral[(delta y)/y d ln x] \[LessEqual] 1 *)

isotropicQ[dim_, time_, t_] := 
 Module[{c, ct, diff, Delta, clumpdata, clumpdatat},
  clumpdata = Join[{dim}, {time}, List[getAllClumpSizesHelper[t]]];
  clumpdatat = 
   Join[{dim}, {time}, List[getAllClumpSizesHelper[Transpose[t]]]];
  c = myhistogram[clumpdata];
  ct = myhistogram[clumpdatat];
  
  (* minor hack to throw out early data since our integral test \
doesnt work for the initial condition *)
  
  If[ct[[5, 2]] == c[[5, 2]] == 0, Return[False]];
  
  diff = Block[{x, y, xt, yt},
    {x, y} = Transpose[c];
    {xt, yt} = Transpose[ct];
    Transpose[{x, Abs[yt - y]/(1 + yt + y)}]];
  
  Delta = With[{x1 = Min[diff[[All, 1]]], x2 = Max[diff[[All, 1]]],
     int = Interpolation[diff, InterpolationOrder -> 1]},
    NIntegrate[int[x]/x, {x, x1, x2}]];
  
  Return[Delta < 1.0]]
(*function to check if cold has had time to cool to the actual floor \
value.*)
mintempQ[file_] := Module[{temp, time, x, ret},
  time = readhdf5time[file];
  x = time/tlim;
  temp = tempplus[file];
  ret = If[Min[temp] < 1.1*tfloor[x], True, False];
  ret
  ]
findclumpshelper[file_] :=
  Module[{t, ret, dim, time},
   time = readhdf5time[file];
   t = mybinarize[tempplus[file], time];
   dim = Max[Dimensions[t]];
   ret = If[(isotropicQ[dim, time, t] && mintempQ[file]),
     Sort[
      Join[getAllClumpSizesHelper[t], 
       getAllClumpSizesHelper[Transpose[t]]]],
     -1];
   {dim, time, ret}];
findclumps[file_] :=
 Block[{outfile = file <> ".clumps", ret},
  If[FileExistsQ[outfile], Return[Get[outfile]]];
  
  ret = findclumpshelper[file];
  Save[file <> ".clumps", ret];
  Return[ret]]

  (*function to analyze a set of athdf files*)

   SetDirectory["/Volumes/LaCie/simdata"]

   mkfilename[num_] := 
 With[{numstr = 
    ToString[PaddedForm[num, 4, NumberPadding -> {"0", "X"}]]},
  dir <> "/shatter.out3." <> numstr <> ".athdf"]

analyze[files_] :=
 Module[{clumpdata, x, y},
  clumpdata = Map[findclumps[#] &, files];
  clumpdata = Select[clumpdata, Length[#[[-1]]] > 0 &];
  Print["keeping " <> ToString[Length[clumpdata]] <> " / " <> 
    ToString[Length[files]] <> " files."];
  (*bring clump sizes to physical coordinates*)
  
  If[Length[Last[clumpdata]] == 0,
   Return[{{1, 1}}]];
  
  tmp = Map[myhistogram, clumpdata];
  x = tmp[[1, All, 1]];
  y = tmp[[All, All, 2]];
  y = Map[Mean, Transpose[y], {1}];
  Return[Transpose[{x, y}]]]


  (*now, analyze simulations*)

  dir = "test_run";

data = Map[analyze,
   Partition[Map[mkfilename, Range[100]], 20]];

ListLogLogPlot[data, Joined -> True]
