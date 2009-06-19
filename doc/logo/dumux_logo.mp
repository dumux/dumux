%  The DuMuX Logo              Metapost version 1.000 or higher needed!
%  design by Jochen Fritz 2009
%  
%  This is the metapost file for the colored logo. To compile a b/w logo, simply comment the drawoptions
%  in lines 37 and 42. To scale the logo to different sizes, alter u = ... in line 19 but don't forget to
%  alter the fontsizes in line 16 and 35 by the same factor!! 

prologues:=3;
verbatimtex
%&latex
\documentclass[10pt]{article}
\usepackage[T1]{fontenc}
\usepackage{arev}
\usepackage[dvips]{color}
\begin{document}
\fontsize{10mm}{10mm}\selectfont
etex
beginfig(1);
u = 5mm;
z11 = (0u,0u);
z0 = ( 1u,0u);
z1 = ( 3.25u, 1.95u);
z2 = ( 5.5u, 3.75u);
z3 = ( 7u, 1.5u);
z4 = ( 10u, 0u);
z45 = (9.5u, -0.1u);
z5 = ( 11u, 0.5u);
z6 = ( 12u, 0u);
z7 = ( 13u, 0.25u);
z8 = ( 14u, 0u);
z9 = ( 15u, 0.1u);
z10 = ( 16u, 0u);
z12 = ( 17u, 0u);
labeloffset:=0;
label.llft(btex $\textsf{DuMu}^\textsf{\fontsize{6mm}{6mm}\selectfont X} $ etex , (17.2u,3.8u)) withcolor black;
pickup pencircle xscaled 0.3u yscaled 0.05u rotated 10;
drawoptions(withcolor (0.88, 0.58, 0.09));
draw z11{right}..z0..z1..z2;
pickup pencircle xscaled 0.3u yscaled 0.05u rotated 20;
draw z2..z3..{right}z45;
drawoptions(withcolor(0.05,0.45,0.81));
draw z45{right}..z4..z5..z6..z7..z8..z9..z10..z12; 
endfig;
end;