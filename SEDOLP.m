(* ::Package:: *)

Off[CreateDirectory::ioerr];
Print[StyleForm["=====================================================","Section",FontSize->14]]
Print[StyleForm["PACKAGE:","Section",FontSize->14],StyleForm[" SEDOLP ","Section",FontSize->14] ]
Print[StyleForm["Por: Mat. \[CapitalOAcute]scar Iv\[AAcute]n de Jes\[UAcute]s Mungu\[IAcute]a y Dr. Jorge Ch\[AAcute]vez Carlos, (2019)","Section",FontSize->12]]
Print[StyleForm["=====================================================","Section",FontSize->14]]
Print[StyleForm["Link de Notas y descarga:","Section",FontSize->10]]
Print["https://github.com/NuclearGeorge/Notas_EDO_Lineales"]
(*Print["Descarga: https://raw.githubusercontent.com/NuclearGeorge/Notas_EDO_Lineales/\
master/SEDOLP.m"]*)
Print[StyleForm["Este paquete adquiere resuelve: Sistemas de Ecuaciones Diferenciales Ordinarias Lineales Planas, de la forma:","Section",FontSize->12]]
Print[StyleForm["\!\(\*SubscriptBox[\(x\), \(1\)]\)'= a \!\(\*SubscriptBox[\(x\), \(1\)]\) + b \!\(\*SubscriptBox[\(x\), \(2\)]\),  \!\(\*SubscriptBox[\(x\), \(2\)]\)' = c \!\(\*SubscriptBox[\(x\), \(1\)]\) + d \!\(\*SubscriptBox[\(x\), \(2\)]\),  o escrita en forma matricial:","Section",FontSize->12]] 
Print[StyleForm["\!\(\*OverscriptBox[\(x\), \(_\)]\)' =  \!\(\*TagBox[\((\[NoBreak]\*GridBox[{
{a, b},
{c, d}
},\nGridBoxAlignment->{\"Columns\" -> {{Center}}, \"ColumnsIndexed\" -> {}, \"Rows\" -> {{Baseline}}, \"RowsIndexed\" -> {}},\nGridBoxSpacings->{\"Columns\" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, \"ColumnsIndexed\" -> {}, \"Rows\" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, \"RowsIndexed\" -> {}}]\[NoBreak])\),
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\) \!\(\*OverscriptBox[\(x\), \(_\)]\) ","Section",FontSize->12]]
Print[StyleForm["donde {a,b,c,d} son par\[AAcute]metros reales seleccionados por el usuario.","Section",FontSize->12]]
Print["-----------------------------------------------------"];
Print["El paquete fu\[EAcute] cargado exitosamente"];
Print[StyleForm["=====================================================","Section",FontSize->14]]

INPI:=
{a=Input["Introduce el par\[AAcute]metro a = "];
b=Input["Introduce el par\[AAcute]metro b = "];
c1=Input["Introduce el par\[AAcute]metro c = "];
d=Input["Introduce el par\[AAcute]metro d = "];
A={{a,b},{c1,d}};
If[Det[A]== 0,Print["\[DownExclamation] Det A = 0 ! Por favor introduzca unos nuevos par\[AAcute]metros para poder continuar."]];};

INP[ai_,bi_,c1i_,di_]:=
{a=ai;b=bi;c1=c1i;d=di;
A={{a,b},{c1,d}};
If[Det[A]== 0,Print["\[DownExclamation] Det A = 0 ! Por favor introduzca unos nuevos par\[AAcute]metros para poder continuar."]];};


SIS:={val=Eigenvalues[A];
vect=Eigenvectors[A];
\[CapitalDelta]=Tr[A]^2-4Det[A];
sis[x_,y_]:={a x+b y,c1 x+d y};
Print["===================================================================="];
Print["El Sistema de Ecuaciones Diferenciales es: \!\(\*OverscriptBox[\(x\), \(_\)]\)'="MatrixForm[A]," \!\(\*OverscriptBox[\(x\), \(_\)]\)"];
Print["El punto cr\[IAcute]tico es:"];

(*----------Clasificaci\[OAcute]n de estabilidad----------*)
If[\[CapitalDelta]>= 0,{If[val[[1]]== val[[2]],{If[val[[1]]<0,Print["Nodo Degenerado Atractor"],Print["Nodo Degenerado Repulsor"]],vec1=vect[[1]],vec2=Transpose[Eigenvectors[A-val[[1]]{{1,0},{0,1}}]][[2]],L={{val[[1]],1},{0,val[[1]]}},P=Transpose[{-vec1,vec2}],Y ={{Exp[val[[1]]t],t Exp[val[[1]]t]},{0, Exp[val[[1]]t]}}},{If[val[[1]]val[[2]]<0,Print["Nodo Hiperb\[OAcute]lico"],If[val[[1]]<0,Print["Nodo Atractor"],Print["Nodo Repulsor"]]],vec1=vect[[1]],vec2=vect[[2]],L={{val[[1]],0},{0,val[[2]]}},P=Transpose[{vec1,vec2}],Y ={{Exp[val[[1]]t],0},{0, Exp[val[[2]]t]}}}]},{If[Re[val[[1]]]<0,Print["Nodo Espiral Atractor"],If[Re[val[[1]]]==0,Print["Nodo centro"],Print["Nodo Espiral Repulsor"]]],vec1=Re[vect[[1]]],vec2=Im[vect[[2]]],L={{Re[val[[1]]],-Im[val[[1]]]},{Im[val[[1]]],Re[val[[1]]]}},P=Transpose[{vec1,vec2}],Y =Exp[Re[val[[1]]]t]{{Cos[Im[val[[1]]]t],-Sin[Im[val[[1]]] t]},{Sin[Im[val[[1]]] t], Cos[Im[val[[1]]] t]}}}];
Print["Los valores propios del sistema son: ", val];
Print["Forma can\[OAcute]nica de la matriz A: \[CapitalLambda] ="MatrixForm[L]];};


(*----------Soluciones ----------*)
SOL:={Clear[c];,Print["Soluci\[OAcute]n en la base can\[OAcute]nica: \!\(\*OverscriptBox[\(y\), \(_\)]\) =", Y//MatrixForm,"\!\(\*TagBox[TagBox[
RowBox[{\"(\", GridBox[{
{SubscriptBox[\"c\", \"1\"]},
{SubscriptBox[\"c\", \"2\"]}
},\nGridBoxAlignment->{\"Columns\" -> {{Center}}, \"ColumnsIndexed\" -> {}, \"Rows\" -> {{Baseline}}, \"RowsIndexed\" -> {}},\nGridBoxSpacings->{\"Columns\" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, \"ColumnsIndexed\" -> {}, \"Rows\" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, \"RowsIndexed\" -> {}}], \")\"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\)"],Print["Soluci\[OAcute]n en la base \!\(\*OverscriptBox[\(x\), \(_\)]\) = ",X= P.Y.{Subscript[c, 1],Subscript[c, 2]}],SOLY[t_]=Y.{Subscript[c, 1],Subscript[c, 2]};,SOLX[t_]=X};


(*SOLCIINP:={
t0=Input["Introduce el valor del tiempo inicial Subscript[t, 0] = "];
x1i=Input["Introduce la condici\[OAcute]n inicial x1(t0) = "];
x2i=Input["Introduce la condici\[OAcute]n inicial Subscript[x, 2](Subscript[t, 0]) = "];};


SOLCII:={SOLCIINP,solC=Solve[SOLX[t0]\[Equal]{x1i,x2i},{Subscript[c, 1],Subscript[c, 2]}][[1]],Print[solC],SOLCIY[t_]=SOLY[t]/.solC,SOLCIX[t_]=SOLX[t]/.solC};
*)
SOLCI[t0i_,x1ii_,x2ii_]:={{t0,x1i,x2i}={t0i,x1ii,x2ii},solC=Solve[SOLX[t0]=={x1i,x2i},{Subscript[c, 1],Subscript[c, 2]}][[1]],Print[solC],SOLCIY[t_]=SOLY[t]/.solC,SOLCIX[t_]=SOLX[t]/.solC};


(*----------------Elementos gr\[AAcute]ficos--------------*)
EF[xl_,xr_,yl_,yr_]:=StreamPlot[sis[x,y],{x,xl,xr},{y,yl,yr},ImageSize->700,StreamStyle-> Black,StreamPoints-> Fine,Axes->True,AxesLabel->{Subscript[x, 1],Subscript[x, 2]},LabelStyle->Directive[Black,Bold ,50,FontFamily->"Courier New"],PlotRange-> All,Prolog-> {White,Rectangle[Scaled[{-2,-2}],Scaled[{2,2}]]}];
EFO[xl_,xr_,yl_,yr_,ti_,tf_]:=Show[EF[xl,xr,yl,yr],ParametricPlot[SOLCIX[t],{t,ti,tf},PlotStyle->Red,ImageSize->700,Axes->True,AxesLabel->{Subscript[x, 1],Subscript[x, 2]},LabelStyle->Directive[Black,Bold ,50,FontFamily->"Courier New"],PlotRange-> Automatic,Prolog-> {White,Rectangle[Scaled[{-2,-2}],Scaled[{2,2}]]}]];
