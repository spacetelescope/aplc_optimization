#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Chapter:: *)
(*Definitions*)


(* ::Section:: *)
(*Initializations*)


(* ::Subsection:: *)
(*Package loading*)


(*PackageVersion="5.0.9";
Init[PackageVersion]*)
GetWolframOpticsLibrary["current"]


Needs["Developer`"]


Off[SFT::"NewPlan"]
Off[ToLyotPlane::"NewPlan"]
SetOptions[ListLogPlot,Joined->True,Frame->True,Axes->None,PlotStyle->AbsoluteThickness[1.5]];


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Padding for File Name *)


(* ::Text:: *)
(*Function to pad string numbers and complete with zeros in front *)


PaddingForFileNames[x_,ndigit_]:=PaddedForm[x,ndigit,NumberPadding->{"0",""},NumberSigns->{"",""}]//ToString


(* ::Section::Closed:: *)
(*Array Size Increase*)


(* ::Subsection:: *)
(*Size increase*)


(* ::Text:: *)
(*Function to oversize a given array by an integer factor and inclusion of the values from the original array*)


SizeIncrease=Compile[{{image,_Real,2},{coeff,_Integer}},Module[{i,j,i0,j0,n,res},
n=Length[image]*coeff;
res=Table[0.,{n},{n}];
Do[
Do[
i0=Quotient[i+coeff-1,coeff];
j0=Quotient[j+coeff-1,coeff];
res[[i,j]]=image[[i0, j0]];
,{i,1,n}];
,{j,1,n}];
res
]
];


(* ::Section::Closed:: *)
(*Size increase and error diffusion algorithm*)


(* ::Text:: *)
(*Function to oversize a given array by an integer factor and inclusion of values from the original array after error diffusion algorithm (e.g. Dorrer & Zuegel 2007)*)


SizeIncreaseAndErrorDiffusion=Compile[{{image,_Real,2},{coeff,_Integer}},Module[{i,j,n,res,err, SuperImage},
n=Length[image]*coeff;
res=Table[0.,{n},{n}];
err=Table[0.,{n},{n}];
SuperImage=SizeIncrease[image, coeff];
(*Print[Dimensions[SuperImage]];*)
Do[
Do[
If[SuperImage[[i,j]] > 0.5,res[[i,j]]=1., res[[i, j]]=0.];
err[[i,j]]=res[[i,j]]-SuperImage[[i,j]];
If[ i+1 <= n && j-1 <= n  && j-1 >= 1 , SuperImage[[i+1, j-1]]=SuperImage[[i+1, j-1]]+-3./16. *err[[i,j]]];
If[ i+1 <= n && j+0 <= n , SuperImage[[i+1, j+0]]=SuperImage[[i+1, j+0]]+-5./16. *err[[i,j]]];
If[ i+1 <= n && j+1 <= n , SuperImage[[i+1, j+1]]=SuperImage[[i+1, j+1]]+-1./16. *err[[i,j]]];
If[ i +0<= n && j+1 <= n , SuperImage[[i+0, j+1]]=SuperImage[[i+0, j+1]]+-7./16. *err[[i,j]]];
,{j,1,n}];
,{i,1,n}];
res
]
];


(* ::Chapter:: *)
(*Parameters*)


(* ::Text:: *)
(*Pupil size*)


nPup=1000;


ResizingFactor=2;
nPupSuper=ResizingFactor*nPup;
RoundApodKeyword = 0;


(* ::Section:: *)
(*Import*)


(* ::Subsection:: *)
(*Aperture*)


FilepathAperture = "/user/kstlaurent/git/progressive_refinement_coronagraphy/masks/LUVOIR/TelAp_LUVOIR_gap_pad01_bw_ovsamp04_N1000.fits"
Aperture=Import[FilepathAperture,{"Data",1}];
ArrayPlot[Aperture, PlotRange-> All,ColorFunction ->GrayLevel]


(* ::Subsection:: *)
(*Lyot Stop*)


FilepathLyotStop = "/user/kstlaurent/git/progressive_refinement_coronagraphy/masks/LUVOIR/LS_LUVOIR_ID0120_OD0982_no_struts_gy_ovsamp4_N1000.fits"
LyotStop=Import[FilepathLyotStop,{"Data",1}];
ArrayPlot[LyotStop, PlotRange-> All,ColorFunction ->GrayLevel]


(* ::Subsection:: *)
(*Shaped pupil*)


FilepathShapedPupil = "/user/kstlaurent/git/progressive_refinement_coronagraphy/surveys/luvoir_BW10_small_N1000_telserv3/solutions/LUVOIR_N1000_FPM350M0150_IWA0340_OWA01200_C10_BW10_Nlam5_LS_IDD0120_OD0982_no_ls_struts.fits"
ShapedPupil=Import[FilepathShapedPupil,{"Data",1}];
ArrayPlot[ShapedPupil, PlotRange-> All,ColorFunction ->GrayLevel]


(* ::Chapter:: *)
(*Resize with error diffusion algorithm*)


(* ::Subsection:: *)
(*Increase size array*)


(* ::Subsubsection:: *)
(*Pupil size increase for apertures and Lyot stop*)


(* ::Text:: *)
(*Resizing of the different pupil arrays*)


If[ResizingFactor==1,
SuperAperture=Aperture;
SuperLyotStop=LyotStop;,
SuperAperture=Ceiling[SizeIncrease[Aperture, ResizingFactor]];
SuperLyotStop=SizeIncrease[LyotStop, ResizingFactor];
]


(* ::Subsubsection:: *)
(*Pupil size increase for the apodizer*)


(* ::Text:: *)
(*Resizing of the shaped pupil array and possible use of error diffusion algorithm*)


SuperShapedPupil=SizeIncreaseAndErrorDiffusion[ShapedPupil, ResizingFactor];


(* ::Chapter:: *)
(*Export*)


(* ::Section:: *)
(*Fits file*)


(* ::Subsection:: *)
(*Oversize*)


(* ::Subsubsection:: *)
(*Filenames*)


FilenameSuperAperture = StringReplace[FilepathAperture,{"N1000" -> "N"<>PaddingForFileNames[nPupSuper, 4]<>"_Ceiling_EDF"<>PaddingForFileNames[ResizingFactor,2]}]
FilenameSuperLyotStop = StringReplace[FilepathLyotStop,{"N1000" -> "N"<>PaddingForFileNames[nPupSuper, 4]<>"_EDF"<>PaddingForFileNames[ResizingFactor,2]}]
FilenameSuperShapedPupil = StringReplace[FilepathShapedPupil,{"LUVOIR_N1000" -> "LUVOIR_N"<>PaddingForFileNames[nPupSuper, 4]<>"_EDF"<>PaddingForFileNames[ResizingFactor,2]}]


(* ::Subsubsection:: *)
(*Export*)


Export[FilenameSuperAperture, SuperAperture];
Export[FilenameSuperShapedPupil, SuperShapedPupil*SuperAperture];
Export[FilenameSuperLyotStop, SuperLyotStop];
