rhoplus[file_] := Module[{b, loc, pos, data, rho, part},
  
			 b = Import[file, "Data"];
  
			 (* read in logical ordering of blocks (assuming 2D data!) *)
  
			 loc = b[[2]][[All, {1, 2}]];
			 loc = Map[Reverse, loc];
  
			 (* sort blocks *)
			 pos = Ordering[loc];
			 data = b[[3, All, All, 1, All, All]];
			 
			 rho = data[[1, pos, All, All]];
  
			 (* partition meshblocks into a 2D array *)
  
			 part = Length[Select[loc, First[#] == 0 &]];
			 rho = Partition[rho, part];
  
			 (* use arrayflatten to compress into one big array *)
  
			 rho = ArrayFlatten[rho];
  rho
  ]


  presplus[file_] := Module[{b, loc, pos, data, P, part},
  
			    b = Import[file, "Data"];
  
			    (* read in logical ordering of blocks (assuming 2D data!) *)
  
			    loc = b[[2]][[All, {1, 2}]];
			    loc = Map[Reverse, loc];
  
			    (* sort blocks *)
			    pos = Ordering[loc];
			    data = b[[3, All, All, 1, All, All]];
			    
			    P = data[[2, pos, All, All]];
  
			    (* partition meshblocks into a 2D array *)
  
			    part = Length[Select[loc, First[#] == 0 &]];
			    P = Partition[P, part];
  
			    (* use arrayflatten to compress into one big array *)
  
			    P = ArrayFlatten[P];
  P
  ]

  tempplus[file_] := presplus[file]/rhoplus[file]

  (*return a list of all clumps along the x-direction in a single file*)

  getAllClumpSizesHelper[t_] := 
  With[{clumplist = Flatten[Map[sliceclumps, t]]}, 
       N[Select[clumplist, # > 0 &]]];

(*return a list of clump sizes along the 1D skewer slice*)
(*slice is \
presumed to contain only 1s and 0s... we return a list of all the \
 lengths of contiguous patches of 1s*)

sliceclumps[slice_] := 
  If[Max[slice] == 1, 
     Map[Length, Select[Split[slice], First[#] == 1 &]], 0] // N;

  (*function to read time for athena++ athdf files*)

  readhdf5time[file_] :=
 
  Block[{str = OpenRead[file, BinaryFormat -> True], line},
	line = Find[str, "Time"];
	Close[str];
	line = StringSplit[line, "Time"][[-1]];
	str = StringToStream[line];
	SetStreamPosition[str, 36];
	BinaryRead[str, "Real64", ByteOrdering -> +1]]

  (*return a list of all clumps in a single file*)
  (*returns a \
   2-element list,first item is simulation time;second item is a list of \
   all clump sizes in the file*)

  getAllClumpSizesplus[file_] := 
  Module[{data = tempplus[file], t, allsizes, 
	  tfloor},(*first threshold the data on temperature*)
   
    tfloor = Min[data];
    t = Map[If[# < 1.5*tfloor, 1.0, 0.0] &, 
	    data, {2}]; (*seperate int hot,cold gas based on 1.5* 
			 floor temp*)
    (*count clumps by slicing along both x-and y-
     directions*)
   allsizes = 
    Join[getAllClumpSizesHelper[t], 
	 getAllClumpSizesHelper[Transpose[t]]];
    {readhdf5time[file], allsizes}];

(*return a histogram of a list of clump sizes in data with log-x bins*)

myhistogram[data_] := 
  Block[{bins, y, x}, 
   bins = Exp[
	      Range[Log[10], Log[Max[data]], N[1/29*Log[Max[data]/10.0]]]];
	(*bin clump data and normalize*)y = BinCounts[data, {bins}];
	(*NORMALIZE IF NEEDED    y=y/Total[y]//N;   *)
	 (*average bins*)
	 x = Exp[ListConvolve[{0.5, 0.5}, Log[bins]]];
	 Transpose[{x, y}]];

	(*given a simulation direcotry,return a clump size histogram for each \
	 step in the simulation*)
	(*we assume 5 equally-spaced steps in \
	 simulation time.*)
	(*returns a list of 5 histograms,one for each step*)

	makeClumpHistograms[dir_] := 
	Module[{fnames = Drop[FileNames[dir <> "*.athdf"], 1], clumpdata, 
		tlim, segments}, clumpdata = Map[getAllClumpSizesplus, fnames];
	  tlim = Last[clumpdata[[All, 1]]];
	  (*simulation is always split into 5 intervals of equal time.(this \
may not line up with VTK file numbers if there are \
								       restarts!)*)(*start by splitting into five segments based on \
										    simulation time*)
   segments = 
    Table[Select[
		 clumpdata, ((i - 1)/5.0) tlim < #[[1]] <= (i/5.0) tlim &], {i, 
									     5}];
	  (*remove the time elements and flatten each segment into one big \
	   list of clump sizes*)
	  segments = Map[Flatten[#[[All, 2]]] &, segments];
	  (*now,calculate the histogram for each segment and return*)
	  Map[myhistogram, segments]];



	(*format numbers for printing to a text file.*)
	(*whyyyyyy is this so \
	 hard in mathematica???*)

	myfmt[num_] := NumberForm[num, {10, 4}, NumberPadding -> {" ", "0"}];

	(*clean up the data write it to a file*)

	exporthist[dir_, outfile_] := 
	Module[{data, x, table, fmtd, header, nres, line}, 
	       data = makeClumpHistograms[dir];
	       (*all 5 histograms have the same x-column;
		split off and put it as the first column*)x = data[[1, All, 1]];
	       data = data[[All, All, 2]];
	       table = Prepend[data, x];
	       table = Transpose[table];
	       (*format the numbers and add a header*)
	       fmtd = Map[myfmt, table, {2}];
	       header = {"# [1] = size (cells),", " [n+1] = PDF[n]"};
	       fmtd = Prepend[fmtd, header];
	       (*add information about simulation resolution and cstcool sizes \
		for later analysis*)nres = Max[dim];
	       cstcool = nres/(4^(Range[5] - 1));
	       line = Prepend[cstcool, nres];
	       line = Map[ToString[PaddedForm[#, 10]] &, line];
	       header = {"# [1] = simulation resolution,", " [n+1] = cstcool"};
	       fmtd = Join[{header, line}, fmtd];
	       (*write to file*)Export[outfile, Append[fmtd, {}]]];

	(*read the directory from the command-line options*)
	If[Length[ARGV] < 3, 
	   Block[{}, 
		 Print["Usage: mash export_hist.m path-to-vtk-files/ outfile"];
		 Exit[0]]];

	dir = ARGV[[2]];
	outfile = ARGV[[3]];

	exporthist[dir, outfile];


