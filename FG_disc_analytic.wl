(* ::Package:: *)

(*Uncomment first part and comment the list tu use from terminal*)
args = (*$ScriptCommandLine[[2 ;;]]*){"/Users/baptisteguilleminot/Documents/M2/Stage_cleaned/data/SDPB_d=6_ci=0_lmax=16_Nmax=20_minmax=max_improv={TH, subTH, SPC}.txt","6","/Users/baptisteguilleminot/Documents/M2/Stage_cleaned/data/"};
coeffAmpPath = args[[1]];
d = ToExpression[args[[2]]];
outPath = args[[3]];


(* ::Section:: *)
(*Tools FG*)


PolQ[J_, d_, z_]:=(4/Pi) (Sqrt[Pi]Gamma[J+1]Gamma[(d-2)/2])/(2^(J+1) Gamma[J+(d-1)/2]z^(J+1)) Hypergeometric2F1[(J+1)/2,(J+2)/2, J+(d-1)/2,1/z^2];
PolQDer[J_,d_,z_]:=(D[PolQ[J,d,zz],{zz,1}]/.zz->z);
z0t[s_]:= (s+4)/(s-4);
z0u[s_]:= (s+4)/(4-s);
tt[s_,z_]:=((4-s)/2) (1-z);
uu[s_,z_]:=((4-s)/2 )(1+z);


MyNIntegrateFaster[expr_,{x_,a_,b_},p_]:=Module[{temp=expr},(*Print["integrate..."];*)temp=SetPrecision[temp,8*p];
temp=NIntegrate[temp,{x,a,b},Method->"GaussKronrodRule",WorkingPrecision->4*p,PrecisionGoal->p,AccuracyGoal->p,MaxRecursion->100];
temp=SetPrecision[temp,p];
Return[temp];];
HalfCompactify[s0_,\[Phi]_]:=(2 s0)/(1+Cos[\[Phi]]);
CompactIntegrand[x_,var_,s0_]:=D[HalfCompactify[s0,\[Phi]],\[Phi]] (x/.var->HalfCompactify[s0,\[Phi]]);


(* ::Section:: *)
(*Setting up amplitude*)


(* ::Subsection:: *)
(*Rho ansatz*)


ClearAll[gp,wp,mp,nq,I\[Epsilon]];
(* shorthands: *)
gp=256;(* global prec, probably 256 *)
wp=512;(* working prec, probably 512 *)
mp=32;(* mid prec *)
nq=NumericQ;
I\[Epsilon]=SetPrecision[I 10^-wp,gp];

ClearAll[\[Alpha],c];
(* some notations *)
\[Alpha][i___]:=ToExpression[ToString[\[Alpha]]<>StringRiffle[ToString/@{i},"x"]];(* rename \[Alpha][1,2,..]=\[Alpha]1x2x... *)
c[i___]:=ToExpression[ToString[c]<>StringRiffle[ToString/@{i},"x"]];(* rename c[1,2,...]=c1x2x... *)

ClearAll[\[Rho],disc\[Rho],s,t,u,csp,motherpoint];
(* Basic function definitions *)
\[Rho][s_,\[Sigma]_]:=(Sqrt[\[Sigma]-4]-Sqrt[4-s])/(Sqrt[\[Sigma]-4]+Sqrt[4-s]);(* conformal map *)
disc\[Rho][s_,\[Sigma]_]:=(2 Sqrt[s-4] Sqrt[\[Sigma]-4])/(s+\[Sigma]-8)
s[\[Rho]_,\[Sigma]_]=SolveValues[\[Rho][s,\[Sigma]]==\[Rho],s][[1]]//Quiet;(* inverse conf map *)
{t[s_,x_],u[s_,x_]}=SolveValues[{x==1+(2t)/(s-4),s+t+u==4},{t,u}][[1]];(* t and u as func of s and x *)
csp=First@SolveValues[Eliminate[{s+t+u==4,s==t,s==u},{u,t}],s];(* cross sym point *)
motherpoint=SolveValues[\[Rho][csp,\[Sigma]]==0,\[Sigma]][[1]];(* mother point *)

ClearAll[SetPrecGrid,\[Phi]g,\[Rho]g,sg,tg,pg,lg,\[Sigma]g];
(* All grids definitions *)
SetPrecGrid[expr_]:=SetPrecision[expr//Chop[#,10^-gp]&,wp]//Union;
\[Phi]g[n_?nq]:=Table[(\[Pi]/2)*(1+Cos[\[Pi]*i/(n+1)]),{i,n,1,-1}]//SetPrecGrid; (* \[Phi] Chebychev grid *)
\[Rho]g[n_?nq]:=Exp[I*\[Phi]g[n]]//SetPrecGrid;(* corresponding \[Rho] (conformal map) grid *)
sg[n_?nq]:=(s[\[Rho]g[n],motherpoint]//SetPrecGrid)+I\[Epsilon];(* corresponding s (energy squared) grid *)
tg[n_?nq]:=4/\[Pi] {\[Phi]g[Round[n/2]],\[Phi]g[2000][[-Round[n/2];;-1]]}//Flatten//SetPrecGrid;(* t grid for positivity constraint *)
pg[n_?nq]:=Table[1,{i,Length@\[Sigma]g[n]}];(* list of max \[Rho] powers for each \[Sigma]i point. (simplest is all 1, but can be customed to insist on chosen wavelet centers)*)
lg[n_?nq]:=Table[i,{i,0,n,2}];(* even spin list *)
\[Sigma]g[1]:=sg[1]//Re(* wavelet centers grid *)
\[Sigma]g[n_]:=Module[{i=1,j=1,temp,exit,previous,out},
previous=\[Sigma]g[n-1];
For[i=1,i<=100,i++,
For[j=1,j<=Length[sg[i]],j++,
temp=sg[i][[j]]//Re;
If[!MemberQ[previous,temp],Goto[exit]];
];
];
Return["\[Sigma]g error ..."];
Abort[];
Label[exit];
out={previous,temp}//Flatten//Sort;
Return[out//SetPrecGrid]
];


ClearAll[\[Phi],Nd,nl,nd];
(* dim dependant math functions *)
\[Phi][d_,s_]:=(s-4)^((d-3)/4) s^(-1/4);(* partial-wave normalization *)
Nd[d_]:=(16\[Pi])^((2-d)/2)/Gamma[(d-2)/2];(* useful \[Pi]-like normalization term *)
nl[d_,0]:=2^(-3+2 d) \[Pi]^(1/2 (-3+d)) Gamma[1/2 (-1+d)];(* useful normalization term *)
nl[d_,l_]:=((4\[Pi])^(d/2) (d+2l-3)Gamma[d+l-3])/(\[Pi] Gamma[(d-2)/2]Gamma[l+1]);(* useful normalization term *)
nd[d_]:=2^(5-d) nl[d,0];(* threshold inspired term normalization term *)


ClearAll[\[Rho]ansatz,discT\[Rho]ansatz,discU\[Rho]ansatz,cond,Mraw,discTMraw,discUMraw];
(* \[Rho] ansatz: *)
\[Rho]ansatz[s_,t_,u_,n_,m_,\[Sigma]i_,\[Sigma]j_]:=Evaluate[1/6 (Plus@@((\[Rho][#1,\[Sigma]i]^n*\[Rho][#2,\[Sigma]j]^m)&@@@Permutations[{s,t,u},{2}]))];(* cross sym \[Rho] ansatz with arbitrary powers *)
(*This implementation of the different discontinuities works only for n,m\in {0,1}...*)
discT\[Rho]ansatz[s_,t_,u_,n_,m_,\[Sigma]i_,\[Sigma]j_]:=(\[Rho]ansatz[S,T,U,n,m,\[Sigma]ii,\[Sigma]jj]/. {\[Rho][T,aa_]:>disc\[Rho][T,aa],\[Rho][U,aa_]->0,\[Rho][S,aa_]->0,1->0})/. {S->s,T->t,U->u,\[Sigma]ii->\[Sigma]i,\[Sigma]jj->\[Sigma]j}
discU\[Rho]ansatz[s_,t_,u_,n_,m_,\[Sigma]i_,\[Sigma]j_]:=(\[Rho]ansatz[S,T,U,n,m,\[Sigma]ii,\[Sigma]jj]/. {\[Rho][U,aa_]:>disc\[Rho][U,aa],\[Rho][T,aa_]->0,\[Rho][S,aa_]->0,1->0})/. {S->s,T->t,U->u,\[Sigma]ii->\[Sigma]i,\[Sigma]jj->\[Sigma]j}
cond[i_,j_,n_,m_,\[Sigma]i_,\[Sigma]j_,pgi_,pgj_]:=Boole[And@@{Abs[\[Sigma]i]<=Abs[\[Sigma]j],n<=m,(m!=0||j<=1),(n!=0||i<=1),n<=Min[pgi,pgj],m<=Max[pgi,pgj]}];(* condition to avoid double counting *)
Mraw[d_,s_,t_,u_,\[Sigma]g_,pg_]:=
Table[{\[Alpha][i,j,n,m],nd[d]*\[Rho]ansatz[s,t,u,n,m,\[Sigma]g[[i]],\[Sigma]g[[j]]]*(* \[Rho] ansatz *)
cond[i,j,n,m,\[Sigma]g[[i]],\[Sigma]g[[j]],pg[[i]],pg[[j]]]},(* condition to avoid double counting *)
{i,Length@\[Sigma]g},{j,Length@\[Sigma]g},{n,0,pg[[i]]},{m,0,pg[[j]]}]//Flatten[#,2]&;
discTMraw[d_,s_,t_,u_,\[Sigma]g_,pg_]:=
Table[{\[Alpha][i,j,n,m],nd[d]*discT\[Rho]ansatz[s,t,u,n,m,\[Sigma]g[[i]],\[Sigma]g[[j]]]*(* \[Rho] ansatz *)
cond[i,j,n,m,\[Sigma]g[[i]],\[Sigma]g[[j]],pg[[i]],pg[[j]]]},(* condition to avoid double counting *)
{i,Length@\[Sigma]g},{j,Length@\[Sigma]g},{n,0,pg[[i]]},{m,0,pg[[j]]}]//Flatten[#,2]&;
discUMraw[d_,s_,t_,u_,\[Sigma]g_,pg_]:=
Table[{\[Alpha][i,j,n,m],nd[d]*discU\[Rho]ansatz[s,t,u,n,m,\[Sigma]g[[i]],\[Sigma]g[[j]]]*(* \[Rho] ansatz *)
cond[i,j,n,m,\[Sigma]g[[i]],\[Sigma]g[[j]],pg[[i]],pg[[j]]]},(* condition to avoid double counting *)
{i,Length@\[Sigma]g},{j,Length@\[Sigma]g},{n,0,pg[[i]]},{m,0,pg[[j]]}]//Flatten[#,2]&;


(*N[discT\[Rho]ansatz[2,1,1,0,1,1,1]]
discT\[Rho]ansatz[s,t,u,1,0,\[Sigma]i,\[Sigma]j]
discT\[Rho]ansatz[s,t,u,0,1,\[Sigma]i,\[Sigma]j]
discT\[Rho]ansatz[s,t,u,0,0,\[Sigma]i,\[Sigma]j]*)
N[nd[6]]


(* ::Subsection:: *)
(*Divergent piece*)


(* ::Text:: *)
(*It was easier to implement the FG projection on this term without using the amplitude, see below. The fine threshold tuning is not implemented so only even dimension works.*)


ClearAll[\[CapitalDelta],\[CapitalDelta]M];
(* leading order TH term: *)
\[CapitalDelta][d_?nq,s_]:=((4-s)/4)^(-((d-3)/2))*Which[OddQ[d],(\[Pi] (-1)^((d-3)/2))/Log[(4-s)/4] ((Sqrt[\[Sigma]-4]-Sqrt[4-s])/(Sqrt[\[Sigma]-4]+Sqrt[4-s]))^2/.\[Sigma]->8,True,Sin[\[Pi] (d-3)/2]];

\[CapitalDelta][d_,s_,t_,u_]:=\[CapitalDelta][d,s]+\[CapitalDelta][d,t]+\[CapitalDelta][d,u];(* crossing sym threshold function *)

\[CapitalDelta]M[d_,s_,t_,u_]:={{{\[Alpha][0,0,0,0],nd[d]*\[CapitalDelta][d,s,t,u]}}};

(* give all corrections powers *)
ClearAll[p,pstep,\[Delta]\[CapitalDelta],\[Delta]M];
pstep[d_]:=Ceiling[(d-3)/2]-(d-3)/2
p[d_?nq]:=Rationalize@Table[n,{n,-((d-3)/2)+pstep[d],-10^-9,pstep[d]}]/;!OddQ[d];
p[d_?nq]:=Table[i,{i,1,(d-3)/2}]/;OddQ[d]

(* subleading TH terms: *)
\[Delta]\[CapitalDelta][d_?nq,s_,p_]:=Which[OddQ[d],((4-s)/4)^(-((d-3)/2)+p) (\[Pi] (-1)^((d-3)/2))/Log[(4-s)/4] ((Sqrt[\[Sigma]-4]-Sqrt[4-s])/(Sqrt[\[Sigma]-4]+Sqrt[4-s]))^2/.\[Sigma]->8,True,(4-s)^p];
\[Delta]\[CapitalDelta][d_,s_,t_,u_,p_]:=\[Delta]\[CapitalDelta][d,s,p]+\[Delta]\[CapitalDelta][d,t,p]+\[Delta]\[CapitalDelta][d,u,p];

\[Delta]M[d_,s_,t_,u_]:=Table[{\[Alpha][0,0,0,i],nd[d]*\[Delta]\[CapitalDelta][d,s,t,u,p[d][[i]]]},{i,1,Length@p[d]}]


(* ::Subsection:: *)
(*Putting it together*)


ClearAll[M,discTM,discUM];
(*M[d_,s_,t_,u_,\[Sigma]g_,pg_]:=(Mraw[d,s,t,u,\[Sigma]g,pg]+\[CapitalDelta]M[d,s,t,u]+\[Delta]M[d,s,t,u]+\[Epsilon]M[d,s,t,u])*)
M[d_,s_,t_,u_,\[Sigma]g_,pg_]:=Join[Mraw[d,s,t,u,\[Sigma]g,pg],\[CapitalDelta]M[d,s,t,u],{\[Delta]M[d,s,t,u]}(*,\[Epsilon]M[d,s,t,u]*)]//Flatten[#,1]&
discTM[d_,s_,t_,u_,\[Sigma]g_,pg_]:=Join[discTMraw[d,s,t,u,\[Sigma]g,pg](*,discT\[CapitalDelta]M[d,s,t,u]*)(*,\[Delta]M[d,s,t,u],\[Epsilon]M[d,s,t,u]*)]//Flatten[#,1]&
discUM[d_,s_,t_,u_,\[Sigma]g_,pg_]:=Join[discUMraw[d,s,t,u,\[Sigma]g,pg](*,discU\[CapitalDelta]M[d,s,t,u]*)(*,\[Delta]M[d,s,t,u],\[Epsilon]M[d,s,t,u]*)]//Flatten[#,1]&


ClearAll[myM,myDiscTM,myDiscTM];
myM[d_?NumberQ,Nmax_]:=M[Rationalize@d,s,t,u,\[Sigma]g[Nmax],pg[Nmax]];
myDiscTM[d_?NumberQ,Nmax_]:=discTM[Rationalize@d,s,t,u,\[Sigma]g[Nmax],pg[Nmax]];
myDiscUM[d_?NumberQ,Nmax_]:=discUM[Rationalize@d,s,t,u,\[Sigma]g[Nmax],pg[Nmax]];


ClearAll[mysdpbout];
mysdpbout=ToExpression@Import[coeffAmpPath];


(* ::Text:: *)
(*Nmax is hard coded, may want to ease that later. If no comparison Pj, Qj necessary, possible to comment exprSTU.*)


Clear[exprSTU,discTexprSTU,discUexprSTU]
exprSTU[s_,t_,u_]=1/(32\[Pi]) myM[d,20][[All,1]] . myM[d,20][[All,2]]/.mysdpbout[[2]];
discTexprSTU[s_,t_,u_]=1/(32\[Pi]) myDiscTM[d,20][[All,1]] . myDiscTM[d,20][[All,2]]/.mysdpbout[[2]];
discUexprSTU[s_,t_,u_]=1/(32\[Pi]) myDiscUM[d,20][[All,1]] . myDiscUM[d,20][[All,2]]/.mysdpbout[[2]];


Clear[expr,discTExpr,exprSymb, discTExprSymb]
exprSymb[s_,z_]=exprSTU[s, tt[s, z], uu[s, z]];
discTExprSymb[s_,z_]=discTexprSTU[s, tt[s, z], uu[s, z]];
expr[s_,z_]:=exprSymb[Rationalize[s],Rationalize[z]];
discTExpr[s_,z_]:=discTExprSymb[Rationalize[s],Rationalize[z]];
(*discUExpr[s_,z_]=discUexprSTU[s,tt[s,z],uu[s,z]];*)


(* ::Text:: *)
(*Different tests to compare with Julia implementation*)


(*exprSTUrebuilt[ss_, tt_, uu_] := Module[{terms},
  terms = M[d, ss, tt, uu, \[Sigma]g[20], pg[20]];
  (terms[[All, 1]] . terms[[All, 2]] /. mysdpbout[[2]]) / (32 Pi)
]
discTexprSTUrebuilt[ss_, tt_, uu_] := Module[{terms},
  terms = discTM[d, ss, tt, uu, \[Sigma]g[20], pg[20]];
  (terms[[All, 1]] . terms[[All, 2]] /. mysdpbout[[2]]) / (32 Pi)
]
exprRe[s_, z_] := exprSTUrebuilt[Rationalize[s], tt[Rationalize[s], Rationalize[z]], uu[Rationalize[s], Rationalize[z]]]
discTExprRe[s_, z_] := discTexprSTUrebuilt[Rationalize[s], tt[Rationalize[s], Rationalize[z]], uu[Rationalize[s], Rationalize[z]]]
N[exprRe[-20.1-20.1 I,-20.1-20.1 I]]
expr[-20.1-20.1 I,-20.1-20.1 I]
N[discTExprRe[20.1+20.1 I,20.1+20.1 I]]
discTExpr[Rationalize[20.1+20.1 I],Rationalize[20.1+20.1 I]]*)


(*testM=ParallelTable[{x,y,xx,yy,Re[exprRe[x+I y, xx+I yy]],Im[exprRe[x+I y, xx+I yy]]},{x,-20.1, 19.9, 10},{y,-20.1, 19.9, 10},{xx,-20.1, 19.9, 10},{yy,-20.1, 19.9, 10}];
testDiscM=ParallelTable[{x,y,xx,yy,Re[discTExprRe[x+I y, xx+I yy]],Im[discTExprRe[x+I y, xx+I yy]]},{x,-20.1, 19.9, 10},{y,-20.1, 19.9, 10},{xx,-20.1, 19.9, 10},{yy,-20.1, 19.9, 10}];
Export[outPath<>"testM.txt", testM];
Export[outPath<>"testDiscM.txt", testDiscM];*)


(*testM=ParallelTable[{x,y,xx,yy,Re[expr[x+I y, xx+I yy]],Im[expr[x+I y, xx+I yy]]},{x,-20.1, 19.9, 10},{y,-20.1, 19.9, 10},{xx,-20.1, 19.9, 10},{yy,-20.1, 19.9, 10}];
testDiscM=ParallelTable[{x,y,xx,yy,Re[discTExpr[x+I y, xx+I yy]],Im[discTExpr[x+I y, xx+I yy]]},{x,-20.1, 19.9, 10},{y,-20.1, 19.9, 10},{xx,-20.1, 19.9, 10},{yy,-20.1, 19.9, 10}];
Export[outPath<>"testM.txt", testM];
Export[outPath<>"testDiscM.txt", testDiscM];*)


Print["Amplitude loaded"]


(*discTNumexprSTU[s_,t_,u_]=(exprSTU[s,t+10^-12I,u]-exprSTU[s,t-10^-12I,u])/(2I);
discTNumexprSTU[20,30,4-50]-discTexprSTU[20,30,4-50]*)


(* ::Section:: *)
(*Setting up FG*)


Clear[exprJmain1,exprJkeyhole,exprJthresh,exprJ]
p0=16;
(*Integration of the regular part (with rho ansatz)*)
exprJmain1[J_,s_]:=MyNIntegrateFaster[(PolQ[J,d,z]discTExpr[s,z])//CompactIntegrand[#,z,z0t[s]]&,{\[Phi],0,\[Pi]},p0];
(*exprJmain2[J_,s_]:=MyNIntegrateFaster[(PolQ[J,d,z]discUExpr[s,z])//CompactIntegrand[#,z,z0u[s]]&,{\[Phi],0,\[Pi]},p0];*)
Clear[exprJkeyhole,exprJthresh]
(*Definition of the function computing the divergent threshold terms when n half integer*)
exprJkeyhole[J_,s_,n_]:=Module[{div,divT,exprSTU,expr,DiscExprT,exprJmain,exprJkeyhole1},
div[ss_]:=(4-ss)^(-n-1/2);
divT[ss_]:=(-1)^n (ss-4)^(-n-1/2);
exprSTU[ss_,t_,u_]:=div[ss]+div[u]+div[t];
expr[ss_,z_]:=exprSTU[ss, tt[ss,z],uu[ss,z]];
DiscExprT[ss_,z_]:=divT[tt[ss,z]];
exprJmain[JJ_,ss_]:=MyNIntegrateFaster[(PolQ[JJ,d,z]DiscExprT[ss,z])//CompactIntegrand[#,z,z0t[ss]+Exp[I Arg[-1+z0t[ss]]]10^-5]&,{\[Phi],0,\[Pi]},p0];
exprJkeyhole1[JJ_, ss_, ord_: 8] := Module[
	  {z0, dir, eps, \[Alpha], c0, qk},
	  z0  = z0t[ss];
	  dir = Exp[I Arg[z0 - 1]];   (* direction of the cut / keyhole ray *)
	  eps = 10^-5;
	  \[Alpha]   = n +1/2;
	  (* regular prefactor of the discontinuity at the branch point *)
	  c0 = Limit[(z - z0)^\[Alpha] DiscExprT[ss, z],z -> z0,Direction -> dir];
	  qk[k_] := SeriesCoefficient[PolQ[JJ, d, z], {z, z0, k}];
	  c0*Sum[qk[k]*dir^(k+1- \[Alpha])*eps^(k + 1 - \[Alpha])/(k + 1 - \[Alpha]),{k, 0, ord}]
	];
exprJmain[J,s]+exprJkeyhole1[J,s]
];
(*Computation of both integer and half integer divergent terms*)
exprJthresh[J_,s_]:=Module[{Residue, Keyhole},
Residue=Sum[If[Boole[EvenQ[d]]-2(n-Floor[(d-3)/2])==0, Sin[\[Pi] (d-3)/2],1]\[Alpha][0,0,0,Boole[EvenQ[d]]-2(n-Floor[(d-3)/2])] (Pi (-1)^(n+1))/(n-1)! (D[PolQ[J,d,z],{z,n-1}]/.z->z0t[s])((s-4)/2)^-n,{n,Floor[(d-3)/2],1,-1}];
Keyhole=Sum[If[Boole[OddQ[d]]-2(n-Floor[(d-4)/2])==0, Sin[\[Pi] (d-3)/2],1]\[Alpha][0,0,0,Boole[OddQ[d]]-2(n-Floor[(d-4)/2])]exprJkeyhole[J,s,n],{n,Floor[(d-4)/2],0,-1}];
nd[d]/(32 Pi) (Residue+Keyhole/.mysdpbout[[2]])
];
(*Full expression of the projection.*)
(*exprJ[J_,s_]:=(*If[Re[s]>4,*)exprJmain1[Rationalize[J,10^-p0],Rationalize[s,10^-p0]]+exprJthresh[Rationalize[J,10^-p0],Rationalize[s,10^-p0]](*,exprJmain2[J,s]+exprJthresh[J,s]]*)*)
exprJ[J_, s_] := Module[{Jrat, srat, prec = p0},
  Jrat = SetPrecision[Rationalize[J, 10^-prec], prec];
  srat = SetPrecision[Rationalize[s, 10^-prec], prec];
  N[exprJmain1[Jrat, srat], prec] + N[exprJthresh[Jrat, srat], prec]
];


(*Definition of the usual projection with Pj*)
Clear[Pgen,IPt]
Pgen[J_,d_,z_]:=Hypergeometric2F1[-J,J+d-3,(d-2)/2,(1-z)/2];
IPt[J_,s_]:=MyNIntegrateFaster[Sin[\[Theta]]*(1-Cos[\[Theta]]^2)^((d-4)/2) Pgen[J,d,Cos[\[Theta]]]*expr[s,Cos[\[Theta]]],{\[Theta],0,Pi},p0]


(* ::Text:: *)
(*Testing the case J even.*)


J=2;
testS=50.301+60.201 I ;
Timing[a=exprJmain1[J,s]+exprJthresh[J,s]]
Timing[b=IPt[J,Rationalize[testS]]]
Abs[(a-b)/(a+b)]


(* ::Text:: *)
(*Building comparison for Julia code*)


(*testJ=ParallelTable[{x,y,xx,yy,Re[exprJ[x+I y, xx+I yy]],Im[exprJ[x+I y, xx+I yy]]},{x,0, 20, 5},{y,-10, 10, 10},{xx,-20.1, 19.9, 10},{yy,-20.1, 19.9, 10}];
Export[outPath<>"testJ.txt", testJ];*)


(* ::Section:: *)
(*Setting the grid computation*)


(*Function to compute the J-plane (adaptiveGridS) and s-plane (adaptiveGridJ) plots efficiently (adaptive meshing).*)
Clear[adaptiveGridS,adaptiveGridJ]
adaptiveGridS[exprJNum_,s_,largeRes_,smallRes_,xRange_, yRange_,coarseIn_:Automatic]:=Module[{coarse,vals,dx,dy,grad,thresh,mask,inds,cells,refined},(*coarse grid*)
coarse=If[coarseIn===Automatic,
DistributeDefinitions[exprJNum];ParallelTable[{x,y,exprJNum[Rationalize[x+I y,10^-p0],Rationalize[s,10^-p0]]},{x,xRange[[1]],xRange[[2]],largeRes},{y,yRange[[1]],yRange[[2]],largeRes}],
coarseIn];vals=Log[Abs[coarse[[All,All,3]]]];
dx=Abs[Differences[vals,{1}]][[;;,1;;-2]];
dy=Abs[Differences[vals,{0,1}]][[1;;-2,;;]];
grad=Sqrt[dx^2+dy^2]/(1+Abs[(vals[[1;;-2,1;;-2]]+vals[[2;;,1;;-2]]+vals[[1;;-2,2;;]]+vals[[2;;,2;;]])/4]);
thresh=Quantile[DeleteCases[Flatten[grad],Indeterminate],0.9];
mask=Map[#>thresh&,grad,{2}];
inds=Position[mask,True];
cells=({xRange[[1]]+largeRes (#[[1]]-1),yRange[[1]]+largeRes (#[[2]]-1)}&)/@inds;
Print["Refining ",Length[cells]," cells"];
DistributeDefinitions[exprJNum];refined=ParallelMap[Function[{pt},With[{x=pt[[1]],y=pt[[2]]},Table[{xx,yy,exprJNum[Rationalize[xx+I yy,10^-p0],Rationalize[s,10^-p0]]},{xx,x,x+(largeRes-smallRes),smallRes},{yy,y,y+(largeRes-smallRes),smallRes}]]],cells];
DeleteDuplicatesBy[Join[Flatten[coarse,1],Flatten[refined,2]],#[[1;;2]]&]
]
adaptiveGridJ[exprJNum_,J_,largeRes_,smallRes_,xRange_, yRange_,coarseIn_:Automatic]:=Module[{coarse,vals,dx,dy,grad,thresh,mask,inds,cells,refined},(*coarse grid*)
coarse=If[coarseIn===Automatic,
DistributeDefinitions[exprJNum];ParallelTable[{x,y,exprJNum[Rationalize[J,10^-p0],Rationalize[x+I y,10^-p0]]},{x,xRange[[1]],xRange[[2]],largeRes},{y,yRange[[1]],yRange[[2]],largeRes}],
coarseIn];vals=Log[Abs[coarse[[All,All,3]]]];
dx=Abs[Differences[vals,{1}]][[;;,1;;-2]];
dy=Abs[Differences[vals,{0,1}]][[1;;-2,;;]];
grad=Sqrt[dx^2+dy^2]/(1+Abs[(vals[[1;;-2,1;;-2]]+vals[[2;;,1;;-2]]+vals[[1;;-2,2;;]]+vals[[2;;,2;;]])/4]);
thresh=Quantile[DeleteCases[Flatten[grad],Indeterminate],0.9];
mask=Map[#>thresh&,grad,{2}];
inds=Position[mask,True];
cells=({xRange[[1]]+largeRes (#[[1]]-1),yRange[[1]]+largeRes (#[[2]]-1)}&)/@inds;
Print["Refining ",Length[cells]," cells"];
DistributeDefinitions[exprJNum];refined=ParallelMap[Function[{pt},With[{x=pt[[1]],y=pt[[2]]},Table[{xx,yy,exprJNum[Rationalize[J,10^-p0],Rationalize[xx+I yy,10^-p0]]},{xx,x,x+(largeRes-smallRes),smallRes},{yy,y,y+(largeRes-smallRes),smallRes}]]],cells];
DeleteDuplicatesBy[Join[Flatten[coarse,1],Flatten[refined,2]],#[[1;;2]]&]
]


CloseKernels[];
LaunchKernels[12];
(*Timing[dataAdapt=Table[Print["Start coarse grid."];
						adaptiveGrid[exprJ,s+0.01I, 0.2, 0.05,{0.2,8},{-4,4}],{s,-5,50,5}]]*)
dataAdapt=Table[adaptiveGridJ[exprJ,JJ, 0.5, 0.25,{-10.001,100.001},{0.001,100.001}],{JJ,2,4,0.1}]

CloseKernels[];




Export[outPath<>"6dMinN20test.mx", dataAdapt];


ListDensityPlot[dataIPt[[1]]/. {x_,y_,z_}:>{x,y,Log[Abs[z]]},MeshFunctions->{#3&},Mesh->20,PlotLegends->Automatic,ImageSize->Large]


ListDensityPlot[dataFG[[1]]/. {x_,y_,z_}:>{x,y,Log[Abs[z]]},MeshFunctions->{#3&},Mesh->20,PlotLegends->Automatic,ImageSize->Large]


Manipulate[ListDensityPlot[dataAdapt[[k]]/. {x_,y_,z_}:>{x,y,Log[Abs[z]]},MeshFunctions->{#3&},Mesh->20,PlotLegends->Automatic,ImageSize->Large],{k,1,5,1}]
