//###############################################
// Include this macro (setSPRACEcolors.C) in 
//	gROOT->ProcessLine(".x setSPRACEColors.C+");
// and in your
//	.rootlogon.C
// The colors can be accessed by name or number:
//	histo->SetFillColor(SPcolors[1]); // OK
//	histo->SetFillColors(SPcolors[SPdred]); //OK
//###############################################

   enum SPcolor {
      SPred,
      SPdred,
      SPgreen,
      SPdgreen,
      SPblue,
      SPdblue,
      SPyellow,
      SPdyellow,
      SPbeige,
      SPdbeige,
      SPbrown,
      SPdbrown,
      SPgrey,
      SPdgrey,
      SPsgrey,
      SPpurple,
   };

   int SPcolors[16] = {
      TColor::GetColor(177, 4, 14),    TColor::GetColor(140, 21, 21),  TColor::GetColor(0, 155, 118),
      TColor::GetColor(23, 94, 84),    TColor::GetColor(0, 152, 219),  TColor::GetColor(0, 124, 146),
      TColor::GetColor(234, 171, 0),   TColor::GetColor(178, 111, 22), TColor::GetColor(210, 194, 149),
      TColor::GetColor(179, 153, 93),  TColor::GetColor(141, 60, 30),  TColor::GetColor(94, 48, 50),
      TColor::GetColor(182, 177, 169), TColor::GetColor(77, 79, 83), TColor::GetColor(232, 237, 229),  TColor::GetColor(83, 40, 79)
   };

