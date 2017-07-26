data = Import["../data/frag2048.dat"];

(* drop headers *)
data = Select[data, VectorQ[#, NumberQ] &];

(* first line contains resolution, followed by predicted values for cstcool *)
nres = First[First[data]];
cstcool = Rest[First[data]];
data = Rest[data];

(* convert back to {x,y} pairs and drop zeros *)
data = Transpose[data];

data = With[{x = First[data], d = Rest[data]},
   Map[Transpose[{x, #}] &, d]];

data = Map[Select[#, #[[2]] > 0 &] &, data];


(* make sure we've parsed the file correctly *)
Export["histogram-check.pdf",
       ListLogLogPlot[data, Joined -> True]];

Print[""];

Print["nres:"];
Print[nres];

Print[""];
Print["cstcool:"];
Print[cstcool];

Print[""];


(* taken from:
https://mathematica.stackexchange.com/questions/91784/how-to-find-numerically-all-roots-of-a-function-in-a-given-range
*)

(* NDSolve has an adaptive step size... use it to systematically find
   all the roots of the function.  This is much more reliable than
   using, eg, FindRoot! *)

rootSearchD[f_, x1_, x2_, ops : OptionsPattern[]] :=
    Block[{},
          Last[Last[Reap[
              NDSolve[{y'[x$] == f'[x$], 
                       y[x1] == f[x1]}, 
                      y, {x$, x1, x2},
                      Method -> {"EventLocator", "Event" -> y[x$],
                                 "EventAction" :> Sow[x$]},
                      ops]]]]];

Options[rootSearchD] = Options[NDSolve];


myplot[d_, prediction_] :=
    Module[{nlm, rng, tpad, bpad},
  
           (* have to use a different fitting function before cooling
              kicks in *)
           nlm = If[prediction < nres,
                    NonlinearModelFit[Log[d],
                                      { amp + a x - b Log[(1 + (Exp[x]/x0)^n)^(1/n)],
                                        0 < n < 10,
                                        1 < x0 < 2048,
                                        -2 < a < 0,
                                        0 < b < 5},
                                      {amp, a, b, n, x0}, x],
                    NonlinearModelFit[Log[d],
                                      { amp + a x,
                                        -2 < a < 0},
                                      {amp, a, b, n, x0}, x]];
  
           (* knee is defined where g''(x) > 1 *)
           rng = If[prediction < nres,
                    rootSearchD[-nlm''[Log[#]] - 1 &, 1, 2048],
                    {2 nres, 3 nres}];
  
           (* choose padding so that top and bottom plots line up *)
           tpad = {{30, 5}, {0, 5}};
           bpad = {{30, 5}, {20, 0}};
  
           Column[{
               (* top shows a plot with the data and fitting function *)
               LogLogPlot[Exp[nlm[Log[x]]], {x, 1, 2048},

                          (* show measured points, along with the
                             predicted value cstcool *)
                          Epilog -> {Point[Log[d]],
                                     {Thick, Dashed, Lighter[Blue, 0.6], 
                                      Line[{Log[{prediction, 3 10^-5}], Log[{prediction, 10}]}]}},

                          (* mark the knee with a gray rectangle *)
                          Prolog -> {Opacity[0.1], 
                                     Rectangle[Log[{Min[rng], 3 10^-5}], Log[{Max[rng], 10}]]},

                          PlotRange        -> {3 10^-5, 10^1}, 
                          Frame            -> True, 
                          Axes             -> False,
                          AspectRatio      -> 1,
                          PlotRangePadding -> None,
                          FrameTicks       -> {{Automatic, Automatic}, {Automatic, Automatic}},
                          ImageSize        -> 250, 
                          ImagePadding     -> tpad],
    
               (* bottom shows a plot with second derivative of fitting function *)
               LogLinearPlot[{-nlm''[Log[x]], 1}, {x, 1, 2048},

                             (* mark knee with a gray rectangle *)
                             Prolog -> {Opacity[0.1], 
                                        Rectangle[{Log[Min[rng]], 0}, {Log[Max[rng]], 100}]},

                             (* mark predicted value cstcool *)
                             Epilog -> {Thick, Dashed, Lighter[Blue, 0.6], 
                                        Line[{{Log[prediction], 0}, {Log[prediction], 100}}]},

                             PlotRange        -> All, 
                             Frame            -> True, 
                             Axes             -> False,
                             AspectRatio      -> 1/5,
                             PlotRangePadding -> None,
                             FrameTicks       -> {{None, None}, {Automatic, Automatic}},
                             ImageSize        -> 250, 
                             ImagePadding     -> bpad]},

                  Center, Scaled[0.0], ItemSize -> Full]]

p1 = Row[MapThread[myplot, {data, cstcool}]];

Export["scale.pdf", p1, "PDF"];
