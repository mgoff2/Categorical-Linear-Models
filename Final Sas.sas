data Dusty;
	input Condition $ Years $ Smoking $ RAN RAY;
	Groupsize=RAN+RAY;
		datalines;
Dusty lt10 Smoke 203 30
Dusty lt10 Nsmke 120 8
Dusty ge10 Smoke 160 57
Dusty ge10 Nsmke 85 11
Ndust lt10 Smoke 1340 14
Ndust lt10 Nsmke 1004 12
Ndust ge10 Smoke 1360 24
Ndust ge10 Nsmke 985 10
;
run; 

*Full Model; 

proc logistic data=Dusty order=data;
	class Condition (ref="Dusty") Years(ref="lt10") Smoking (ref="Nsmke")/param=ref;
	model RAY/Groupsize = condition years smoking condition*years condition*smoking years*smoking / alpha=.1;
	run; 

*Reduced Model; 

	proc logistic data=Dusty order=data;
	class Condition(ref="Ndust") Years(ref="lt10") Smoking (ref="Nsmke")/param=ref;
	model RAY/Groupsize = condition years smoking/ alpha=.1;
	run; 


*Variations on the reduced model; 

	proc logistic data=Dusty order=data;
	class Condition(ref="Ndust") Years(ref="lt10") Smoking (ref="Nsmke")/param=ref;
	model RAY/Groupsize = condition years smoking / alpha=.1;
	oddsratios condition;
	run; 

*To the reduced model, testing: 
*+beta4 = years*smoking =42.744
*+beta5 condition*years =41.643
*+beta6 = condition*smoking =41.137
*+beta4,5 =40.135
*+beta4,6 =38.991
*+beta5,6 =38.394;

*P-Value for ReducedvsFull;

	data pchi; 
		pchiqs=1-probchi(44.464-40.135, 2);
		proc print;
		run;


*Goodness of Fit;


	proc logistic data=Dusty order=data;
	class Condition(ref="Ndust") Years(ref="lt10") Smoking (ref="Nsmke")/param=ref;
	model RAY/Groupsize = condition years smoking/ scale=none aggregate=(Condition Years Smoking);
	output out=temp p=pred;
	run; 
	proc print data=temp;
	run;

*Estimated Probability;

proc logistic data=Dusty order=data;
	class Condition(ref="Ndust") Years(ref="lt10") Smoking (ref="Nsmke")/param=ref;
	model RAY/Groupsize = condition years smoking;
	output out=disprob p=phat upper=ucl lower=lcl  ;
	run; 
	proc print data=disprob;
	run;

	proc logistic data=Dusty order=data;
	class Condition(ref="Ndust") Years(ref="lt10") Smoking (ref="Nsmke")/param=ref;
	model RAY/Groupsize = condition years smoking;
	output out=disprob2 p=phat upper=ucl lower=lcl predprobs=(individual crossvalidate);
	run; 
	proc print data=disprob2;
	run;


	*############################;

*Question 3; 



data Arthritis;
	input Sex $ Treatment $ Outcome $ count;
	datalines;
Female Active Marked 50
Female Active Some 15
Female Active None 20
Female Placebo Marked 18
Female Placebo Some 21
Female Placebo None 56
Male Active Marked 15
Male Active Some 7
Male Active None 21
Male Placebo Marked 3
Male Placebo Some 2
Male Placebo None 30
;
run; 

*Saturated Model;

proc logistic data=Arthritis;
	weight count;
	class Sex (ref="Male") Treatment (ref="Placebo")/ param=ref;
	model Outcome (ref="None") =Treatment Sex Treatment*Sex/link=glogit scale=none aggregate=(Treatment Sex); 
	output out=prob predprob=I;
	run;
	proc print data=prob;
	run;

*Reduced Model;

proc logistic data=Arthritis;
	weight count;
	class Sex (ref="Male") Treatment (ref="Placebo")/ param=ref;
	model Outcome (ref="None") =Treatment Sex/link=glogit scale=none aggregate=(Treatment Sex); 
	output out=prob predprob=I;
	run;
	proc print data=prob;
	run;

	*LRT = 1.103, 2 chi-sq df. P=;

data pvalue2;
	pval=1-probchi(1.103, 2);
run;
proc print;
run;

*Conditional Independence of sex and response, given treatment;

proc logistic data=Arthritis;
	weight count;
	class Treatment (ref="Placebo")/ param=ref;
	model Outcome (ref="None") =Treatment/link=glogit; 
	run;



data pvalue2;
	pval=1-probchi(17.137, 2);
run;
proc print;
run;



*Ordinal model;

proc logistic data=Arthritis order=data;
	weight count;
	Class Treatment Sex; 
	model Outcome = Treatment Sex/link=clogit;
	run; 

	proc logistic data=Arthritis descending;
	weight count;
	Class Treatment Sex; 
	model Outcome = Treatment Sex/link=clogit;
	run; 




	data pvalue2;
	pval=1-probchi(17.137, 2);
run;
proc print;
run;



*############################;
*Section 4;

Data Danish;
	input Xmgm $ Yspv $ Zwkr $ count;
	datalines;
Bad Low Low 103
Bad Low High 87
Bad High Low 35
Bad High High 42
Good Low Low 59
Good Low High 109
Good High Low 78
Good High High 206
;
run;

*Model X, Y, Z;

proc genmod data=Danish order=data;
	class Xmgm Yspv Zwkr; 
	model count= Xmgm Yspv Zwkr /dist=poisson;
	output out=temp1 p=phat;
	run; 

*Model XY, Z;

proc genmod data=Danish order=data;
	class Xmgm Yspv Zwkr Xmgm*Yspv; 
	model count= Xmgm Yspv Zwkr /dist=poisson;
	output out=temp2 p=phat;
	run;

*Model XY, YZ;

	proc genmod data=Danish order=data;
	class Xmgm Yspv Zwkr; 
	model count= Xmgm Yspv Zwkr Xmgm*Yspv Yspv*Zwkr /dist=poisson;
	output out=temp3 p=phat;
	run;

*Model XY, YZ, XZ;

proc genmod data=Danish order=data;
	class Xmgm Yspv Zwkr; 
	model count= Xmgm Yspv Zwkr Xmgm*Yspv Xmgm*Zwkr Yspv*Zwkr/dist=poisson;
	output out=temp4 p=phat;
	run; 

*Model XYZ;

proc genmod data=Danish order=data;
	class Xmgm Yspv Zwkr; 
	model count= Xmgm Yspv Zwkr Xmgm*Yspv Xmgm*Zwkr Yspv*Zwkr Xmgm*Yspv*Zwkr /dist=poisson;
	output out=temp1 p=phat;
	run;

*Part c;

proc print data=temp4;
run;


data pvalue;
	pv=1-probchi(0.0003,1);
	proc print;
	run; 



*Part e;

proc genmod data=Danish order=data;
	class Xmgm Yspv Zwkr; 
	model count= Xmgm Yspv Zwkr /dist=poisson;
	output out=temp1 p=phat;
	run; 

data pvalue2;
pv=1-probchi(116.2776, 4);
proc print;
run;

data pvalue2;
pv=1-probchi(125.8666, 4);
proc print;
run;
