
/*************************************************************************************************************/
/*                                       MACRO SPM  version 1.1 (May. 2015)                                  */ 
/*************************************************************************************************************/

**************************************************************************************************************;
*	The Joint Modelling of Longitudinal Outcomes and Multiple Competing Risk Dropouts		                  ;
*	Code written by Wei Wang and Michael E. Griswold(2015)													  ;
* 																					                          ;
*	Reference:																		                          ;
* 	Wang W, Wang W, Mosley T.H. and Griswold M.E. (2015). A SAS macro for the joint modelling                 ;
*   of longitudinal outcomes and multiple competing risk dropouts. Computer Methods and Programs 			  ;	
* 	in Biomedicine. Submitted.																				  ;
* 																					                          ;
**************************************************************************************************************;
*																					                          ;
*	PROGRAM DESCRIPTION:															                          ;
*   --------------------------------------------------------------------------------------------------------- ;
*	This program implement a joint model for the longitudinal outcome and cause-specific dropouts by assuming ;
*   that the association between the survival and the longitudinal submodels are linked through the shared    ;
*   random effects. The submodel for the longitudinal outcome takes the form of a linear mixed effects,       ;
*   allowing flexibility of random intercept and/or random slope. The shared parameter model encoded in this  ; 
*   SAS macro can accommodate up to three different causes for dropout, and two parametric baseline hazard    ;
*   functions are available for the survival component, the exponential or Weibull distribution. In addition, ; 
*   information criterion fit statistics AIC and BIC are provided for parametric baseline hazard function     ;
*   selection.                                                                                                ;
*																					                          ;
*																					                          ;
*	DATASET REQUIREMENTS:															                          ;
*   --------------------------------------------------------------------------------------------------------- ;
*	The macro requires a SAS data set with the following variables:                                           ;
*   (1) subject identification variable (ID)                                                                  ; 
*   (2) longitudinal outcome variable (LDAOUTCOME)                                                            ; 
*   (3) list of predictors for the longitudinal outcome (LDAPRED)                                             ; 
*   (4) dropout indicator (RISKIND)                                                                           ; 
*   (5) survival/censoring time (SURVTIME)                                                                    ; 
*   (6) list of predictors for the survival submodels (SURVPRED).                                             ; 
*	The dataset is assumed to be entered with ni records for individual i, and of note, our SPM macro only    ; 
*   allows time-independent covariates for survival component prediction, so the variables RISKIND, SURVTIME  ;
*   and SURVPRED should be consistent across all ni observations for the same individual i. However,          ;
*   predictors for the longitudinal outcome can be time-varying (see LDAPRED2 in                              ;
*   the data set example).                                                                                    ;
*                                                                                                             ;
*   Format of data set exmaple:                                                                               ;
*                                                                                                             ;
*	ID	TIME LDAOUTCOME LDAPRED1 LDAPRED2 RISKIND SURVTIME SURPRED1 SURPRED2                    			  ;
*	1	0	 13.2		1		 21		  0		   22.5    55		0			                              ;
*	1	6	 12.5		1		 28		  0		   22.5	   55       0					                      ;
*	1	21	 8.3		1		 25		  0		   22.5	   55	    0   				                      ;
*	2	0	 16.9		0		 23		  1		   10.2	   48		1				                          ;
*	2	5	 18.2		0		 15		  1		   10.2	   48		1				                          ;
*	3	0	 12.0		0		 23		  2		   13.5	   32		0					                      ;
*	3	7	 10.2		0		 35		  2		   13.5	   32		0						                  ;
*	4	0	 23.5		1		 12		  0		   21.3	   65		1					                      ;
*	4	4.5	 21.0		1		 13		  0		   21.3	   65		1				                          ;
*	4	20	 12.7		1		 45		  0		   21.3	   65		1						                  ;
*																					                          ;
*																					                          ;
*	MODEL SPECIFICATION (SEE REFERENCE FOR DETAILS)         						                          ;
*   --------------------------------------------------------------------------------------------------------- ;
*	longitudinal submodel: linear mixed model with random intercept and/or random slope   					  ;
*																					                          ;
*	survival submodel: parametric baseline hazard function as exponential or weibull with time-independent    ;
*	                   covariates  	   			    	                                                      ;
*																					                          ;
*	Shared random effect in both longitudinal submodel and survival submodel   								  ;
*													                 				                          ;
*													                 				                          ;
*	MACRO VARIABLES:																                          ;
*   --------------------------------------------------------------------------------------------------------- ;
*	DATA: 		    SAS data set to fit the shared parameter model                                  	      ;
*																					                          ;
*	ID:             Subject identification																	  ;
*																					                          ;
*	NRISK:		    Number of competing risk dropouts, K = 1, 2 or 3	                                      ;
*																					                          ;
*	LDAOUTCOME:		Longitudinal outcome variable                             		                          ;
*																					                          ;
*	LDAPRED:		List of predictors for the longitudinal outcome including time variable             	  ;
*																					                          ;
*	LDATIME:		Time variable for the longitudinal outcome		                                          ;
*																					                          ;
*	RISKIND: 		Dropout indicator, Ki = 1, 2, …, K                                      	              ;
*																					                          ;
*	SURVTIME:       Survival/censoring time																	  ;
*																					                          ;
*	SURVPRED:		List of predictors for the survival submodels (time-independent covariates required)      ;
*																					                          ;
*	RANDOM:		    Random effect indicator                             		                              ;
*					1 = random intercept only in longitudinal submodel						 				  ;
*					2 = random intercept and random slope in longitudinal submodel 							  ;
*																					                          ;
*	NLMOPTS:		User specified PROC NLMIXED procedure option, if missing then default option for NLMIXED  ;
*                   procedure is used for model fitting														  ;
*																					                          ;
*	PRINTALL:		Outpout data set print out or not	                                                      ;
*					F: output data set is not printed out								                      ;
*					T: output data set is printed out, default												  ;
*																					                          ;
*	ALLFIT:			Output dataset 1 includes information criterion information for all models with different ;
*                   combinations of baseline hazard function (exponential or Weibull)			              ;
*																					                          ;
*	ALLPAR:		    Output dataset 2 includes parameter estimates from the shared parameter model and the     ;
*                   separate model for all models with different combinations of baseline hazard function     ;
*                   (exponential or Weibull)                             		                              ;
*																					                          ;
*	LEASTAICPAR:	Output dataset 3 includes parameter estimates from the AIC selected shared parameter      ; 
*                   model and corresponding separate model                                       		      ;
*																					                          ;
*																					                          ;
*	EXAMPLE CODE:																	                          ;
*																					                          ;
*   %include 'SPM.sas'														                                  ;
*																					                          ;
*   %SPM (DATA        = aricdata,                                                                             ;
*         ID          = idnum,                                                                                ;
*         NRISK       = 2,                                                                                    ;
*         LDAOUTCOME  = total_z,                                                                              ; 
*         LDAPRED     = time agev2c agev2ctime female femaletime educ1 educ2 hypert25 hypertime,              ; 
*         LDATIME     = time,                                                                                 ;
*         RISKIND     = cat,                                                                                  ;
*         SURVTIME    = _t,                                                                                   ;
*         SURVPRED    = agev2c female educ1 educ2 hyper25,                                                    ;
*         RANDOM      = 2,                                                                                    ;
*         NLMOPTS     = maxfunc = 50000 maxiter = 5000,                                                       ;
*         PRINTALL    = T,                                                                                    ; 
*         ALLFIT      = out01,                                                                                ; 
*         ALLPAR      = out02,                                                                                ;
*         LEASTAICPAR = out03                                                                                 ;
*        )                                                                                                    ;
*																					                          ;
*																					                          ;
**************************************************************************************************************;

**************************************************************************************************************;
**************************************************Final Macro*************************************************;
**************************************************************************************************************;

** Define Macro;

** SPM1 is Macro to model the longiudinal and survial data jointly using specifying parametric baseline hazard function
SURVdist1, SURVdist2, SURVdist3 as Exponential or Weibull distribution;

%macro SPM1(data       = , 
            id         = ,
            nrisk      = ,
            LDAoutcome = ,
            LDApred    = , 
            LDAtime    = , 
            riskind    = ,
            SURVtime   = ,
            SURVpred   = ,
            Random     = ,
            nlmopts    =,
            outaic =,
            outpar =
            );

*******************************************************************************************************
Step 1: Parameterize the dropout indicator and survival time, and specify the random effects
Re-organize survival time and drop out status, if a subject drops out for one recorded reason then this 
subject is right censored for all other reasons;
*******************************************************************************************************;

data _dat; set &data; 
  _id = &id;
  _Y  = &LDAoutcome;
  if &nrisk  = 3 then do;
      _time = &SURVtime;
      _time1 = &SURVtime;
      _time2 = &SURVtime;
    if &riskind = 0 then do; 
    _event = 0;
      _event1 = 0;
      _event2 = 0;
    end;
    else if &riskind = 1 then do;
      _event = 1;
      _event1 = 0;
      _event2 = 0;
    end;
    else if &riskind = 2 then do;
      _event = 0;
      _event1 = 1;
      _event2 = 0;
    end;
    else if &riskind = 3 then do;
      _event = 0;
      _event1 =0;
      _event2 = 1;
    end;
  end;
  if &nrisk  = 2 then do;
      _time = &SURVtime;
      _time1 = &SURVtime;
    if &riskind = 0 then do; 
      _event = 0;
      _event1 = 0;
    end;
    else if &riskind = 1 then do;
      _event = 1;
      _event1 = 0;
    end;
    else if &riskind = 2 then do;
      _event = 0;
      _event1 = 1;
    end;
  end;
  if &nrisk  = 1 then do;
    _time = &SURVtime;
    _event = &riskind;
  end;
run;

data _dat;
set _dat;
if &ldaoutcome ne .;
run;

proc sort data= _dat; by _id _time; run;

data _dat; set _dat; by _id;
last = last._id; *for non time-varying survival models, just use last obs;
run;

/* Specify Random Effect used */ 

data _null_;
format parameter parameter1 parameter2 parameter3 parameter4 parameter5 parameter6 $50.;
if &random = 1 then do;
parameter = 'intercept';
call symput('randoms',parameter);
parameter1 = 'u0i';
call symput('linranlda',parameter1);
parameter2 = 'r10*u0i';
call symput('linransurv1',parameter2);
parameter3 = 'r20*u0i';
call symput('linransurv2',parameter3);
parameter4 = 'r30*u0i';
call symput('linransurv3',parameter4);
parameter5 = 'u0i ~ normal(0, tau0)';
call symput('ranspe',parameter5);
parameter6 = '';
call symput('covtype',parameter6);
end;
if &random = 2 then do;
parameter = "&ldatime intercept";
call symput('randoms',parameter);
parameter1 = 'u0i + &ldatime*u1i';
call symput('linranlda',parameter1);
parameter2 = "r10*u0i + r11*u1i";
call symput('linransurv1',parameter2);
parameter3 = "r20*u0i + r21*u1i";
call symput('linransurv2',parameter3);
parameter4 = "r30*u0i + r31*u1i";
call symput('linransurv3',parameter4);
parameter5 = 'u0i u1i ~ normal([0, 0],[tau0,tau2,tau1])';
call symput('ranspe',parameter5);
parameter6 = 'type = un';
call symput('covtype',parameter6);
end;
run;

%put &linransurv1;
%put &linranlda;
%put &ranspe;

****************************************************************
** Fit Separate model for comparison purpose and also generate**
** initial estimates for separate models                      **
****************************************************************;

*******************************************************************************************************
Step 2: Estimate a separate linear mixed model fit to the longitudinal data
Fit the LDA part with PROC MIXED (& ML estimation);
*******************************************************************************************************;

title1 'Separate Longitudinal Model - MIXED (MAR)'; 
proc mixed data=_dat noclprint method=ml cl covtest;
  class _id;
  model _Y = &LDApred / s cl; 
  random &randoms / subject=_id &covtype;
  ods output solutionf=beta;
  ods output covparms=cov;
  ods output ConvergenceStatus = con11 (keep = status);
run;

data con12;
  set con11;
  i = 1;
run;

** create dataset of initial LDA submodel parameters for NLMIXED;
** LDA mean model parameters;

data beta1; 
  length parameter $8; 
  set beta END=lastobs;
  format row 3.0 linpred $255.;
  row = _n_-1;
  parameter = "a" || left(row);
*create LDA linear predictor for NLMIXED;
  linpred =  "+" || cats(parameter) || "*" || effect;
  if _N_=1 then call symput('linpredLDA',parameter);
  else call symput('linpredLDA',trim(resolve('&linpredLDA'))||linpred);
run; 
%put &linpredLDA;

data cov1; 
  length parameter $8; set cov;
  format row 3.0;
  row = _n_-1;
  if covparm="Residual" then do;
     parameter="s2";
  end;
  else do;
     parameter = "tau" || left(1-row);
     if covparm="UN(1,1)" then parameter = "tau1";
     if covparm="UN(2,1)" then parameter = "tau2";
     if covparm="UN(2,2)" then parameter = "tau0";
  end;
run; 

*combine LDA parameters;
* Of note, for variance estimate, p value calculated from wald z test, so it is z value and not t value;

data LDAests; 
  length model $6. effect $15.; 
  set beta1 cov1(in = a rename=(covparm=Effect));
  model="LDA";
  row = _n_-1;
  if a then do;
    probt = probz;
    tvalue = zvalue;
  end;
  rename Probt=pvalue;
  i = 1;
  drop probz zvalue;
run; 

data LDAests;
  merge LDAests con12;
  by i;
  drop i;
run;

title;

*******************************************************************************************************
Step 3: Estimate separate parametric survival models to fit the competing risk survival data
Fit the survival model with PROC LIFEREG first then reparameterize with PROC NLMIXED;
*******************************************************************************************************;

title1 "Separate Survival Model 1 - &riskind = 1";

proc lifereg data=_dat(where=(last=1));
  model _time*_event(0)= &SURVpred / dist=&SURVdist1;
  ods output ParameterEstimates=SURVests10(rename=(Parameter=Effect));
  ods output ConvergenceStatus = con21(keep = status);
run;

data con22;
  set con21;
  i = 1;
run;

*Survival model parameters;

data SURVests11; 
  length model $6 parameter $8; 
  set SURVests10;
  format row 3.0;
  row = _n_-1;
  parameter = "b" || left(row);
  if effect="Scale" then parameter="scale1";
  if effect="Weibull Shape" then parameter="shape";
  if DF=0 then delete; *remove unestimated parameters;
  if effect="Weibull Shape" then delete; *for weib, only need scale: shape=1/scale;
  model="Surv1";
  rename LowerCL=Lower UpperCL=Upper ProbChisq=pvalue;
run;
 
*create Surv1 linear predictor for NLMIXED;

data SURVests12; 
  set SURVests11; format linpred $255.;
  if parameter NOT in("scale1", "Scale", "Shape") then do;
    linpred =  "+" || cats(parameter) || "*" || effect;
    if _N_=1 then call symput('linpredSURV1',parameter);
    else call symput('linpredSURV1',trim(resolve('&linpredSURV1'))||linpred);
  end;
  i = 1;
run;
%put &linpredSURV1;
title;

data SURVests12;
  merge SURVests12 con12;
  by i;
  drop i;
run;

*Add intial value for loading factors as 0;

data Rcorr10; 
  model="Surv1"; 
  Effect="Corr";
  parameter="r10"; 
  estimate=0; 
run; 

data Rcorr11; 
  set Rcorr10;
  parameter="r11"; 
  estimate=0; 
run; 

data Rcorr20; 
  set Rcorr10;
  model="Surv2";
  parameter="r20"; 
  estimate=0; 
run; 

data Rcorr21; 
  set Rcorr20;
  parameter="r21"; 
  estimate=0; 
run;
 
data Rcorr30; 
  set Rcorr10;
  model="Surv3";
  parameter="r30"; 
  estimate=0; 
run;
 
data Rcorr31; 
  set Rcorr30;
  parameter="r31"; 
  estimate=0; 
run; 

* Get initial value for re-parameterization to fit weibull/exponential distribution with proc nlmixed; 

%if &SurvDist1 = exponential %then %do;

  data temp13;
    set SURVests12;
    estimate = estimate * (-1);
    keep model effect parameter estimate;
  run;

%end;

%if &SurvDist1 = weibull %then %do;

  data temp11;
    set SURVests12;
    if parameter = 'scale1';
    keep estimate i;
    i = 1;
    rename estimate = estimate1;
  run;

  data temp12;
    set SURVests12;
    i = 1;
  run;

  data temp13;
    merge temp11 temp12;
    by i;
    if substr(parameter, 1, 1) = 'b' then estimate = -estimate / estimate1;
    keep model effect parameter estimate;
  run;

%end;

data _lastdat;
  set _dat;
  if last = 1;
run;

proc nlmixed data=_lastdat &nlmopts;

  parms / data=temp13;
  linpredSURV1  = &linpredSURV1;
  *exponential survival;
  %if &SurvDist1 = exponential %then %do;
    alph  = exp(linpredsurv1);
    loglikeSURV1 = (_event=1)*(linpredsurv1-alph*_time)+ (_event=0)*(-alph*_time);
  %end;

  *Weibull survival;
  %if &SurvDist1 = weibull %then %do;
    gam   = 1/scale1;
    alph  = exp(linpredsurv1);
    loglikeSURV1 = (_event=1)*(-alph*(_time**gam) + log(gam) + linpredsurv1 + (gam-1) * log(_time)) + (_event=0)*(-alph*(_time**gam));
  %end;

  model _time ~ general(loglikeSURV1);

  ods output ParameterEstimates=temp14 (drop = Gradient rename =(StandardError = StdErr Probt = pvalue));
  ods output ConvergenceStatus = temp15 (keep = status rename = (status = status));
run;

data temp16;
  set temp15;
  i = 1;
run;

data temp17;
  set temp14;
  i = 1;
run;

data temp18;
  merge temp17 temp16;
  by i;
  drop i;
run;

data temp19;
  format parameter $8.;
  merge temp18 temp13 (keep = parameter effect);
  by parameter;
run;

data _SPMpar; 
  set LDAests temp19 (in = a) Rcorr10 Rcorr11;
  row = _n_-1;
  if a then model = 'Surv1';
  drop df tvalue alpha subject;
run;

data _SPMinit;
  set _SPMpar;
  keep model effect parameter estimate;
run;

proc datasets lib = work;
  delete SURVests10 SURVests11 SURVests12 con11 con12 con21 con22 temp11 temp12 temp13 temp14 temp15 temp16 temp17 temp18 temp19;
run;

%if  &nrisk  = 2 %then %do;

  data _dat;
    set _dat;
    drop _time2 _event2;
  run;

  title1 "Separate Survival Model 2 - &riskind = 2";
  proc lifereg data=_dat(where=(last=1));
    model _time1*_event1(0)= &SURVpred / dist=&SURVdist2;
    ods output ParameterEstimates=SURVests20(rename=(Parameter=Effect));
    ods output ConvergenceStatus = con31(keep = status);
  run;

  data con32;
    set con31;
    i = 1;
  run;

  * create dataset of initial Survival submodel parameters for NLMIXED;
  *Survival model parameters;
  data SURVests21; 
    length model $6 parameter $8; set SURVests20;
    format row 3.0;
    row = _n_-1;
    parameter = "c" || left(row);
    if effect="Scale" then parameter="scale2";
    if effect="Weibull Shape" then parameter="shape";
    if DF=0 then delete; *remove unestimated parameters;
    if effect="Weibull Shape" then delete; *for weib, only need scale: shape=1/scale;
    model="Surv2";
    rename LowerCL=Lower UpperCL=Upper ProbChisq=pvalue;
  run;
 
  *create Surv2 linear predictor for NLMIXED;
  data SURVests22; 
    set SURVests21; 
    if parameter NOT in("scale2", "Scale", "Shape") then do;
      linpred =  "+" || cats(parameter) || "*" || effect;
      if _N_=1 then call symput('linpredSURV2',parameter);
      else call symput('linpredSURV2',trim(resolve('&linpredSURV2'))||linpred);
    end;
    i = 1;
  run;
  %put &linpredSURV2;
  title;

  data SURVests22;
    merge SURVests22 con32;
    by i;
    drop i;
  run;

  * Get initial value for re-parameterization to fit weibull/exponential distribution with proc nlmixed; 

  %if &SurvDist2 = exponential %then %do;

    data temp23;
      set SURVests22;
      estimate = estimate * (-1);
      keep model effect parameter estimate;
    run;

  %end;

  %if &SurvDist2 = weibull %then %do;

    data temp21;
      set SURVests22;
      if parameter = 'scale2';
      keep estimate i;
      i = 1;
      rename estimate = estimate1;
    run;

    data temp22;
      set SURVests22;
      i = 1;
    run;

    data temp23;
      merge temp21 temp22;
      by i;
      if substr(parameter, 1, 1) = 'c' then estimate = -estimate / estimate1;
      keep model effect parameter estimate;
    run;

  %end;

  proc nlmixed data=_lastdat &nlmopts;
    parms / data=temp23;
    linpredSURV2  = &linpredSURV2;
    *exponential survival;
    %if &SurvDist2 = exponential %then %do;
      alph  = exp(linpredSURV2);
      loglikeSURV2 = (_event1=1)*(linpredSURV2-alph*_time1)+ (_event1=0)*(-alph*_time1);
    %end;

    *Weibull survival;
    %if &SurvDist2 = weibull %then %do;
      gam   = 1/scale2;
      alph  = exp(linpredSURV2);
      loglikeSURV2 = (_event1=1)*(-alph*(_time1**gam) + log(gam) + linpredSURV2 + (gam-1) * log(_time1)) + (_event1=0)*(-alph*(_time1**gam));
    %end;

    model _time1 ~ general(loglikeSURV2);

    ods output ParameterEstimates=temp24 (drop = Gradient rename =(StandardError = StdErr Probt = pvalue));
    ods output ConvergenceStatus = temp25 (keep = status rename = (status = status));
  run;

  data temp26;
    set temp25;
    i = 1;
  run;

  data temp27;
	set temp24;
	i = 1;
  run;

  data temp28;
	merge temp27 temp26;
	by i;
	drop i;
  run;

  data temp29;
    format parameter $8.;
    merge temp28 temp23 (keep = parameter effect);
    by parameter;
  run;

  data _SPMpar; 
    set _SPMpar temp29 (in = a) Rcorr20 Rcorr21;
    row = _n_-1;
	if a then model = 'Surv2';
    drop df tvalue alpha;
  run;

  data _SPMinit;
    set _SPMpar;
    keep model effect parameter estimate;
  run;

  proc datasets lib = work;
    delete SURVests20 SURVests21 SURVests22 con31 con32 temp21 temp22 temp23 temp24 temp25 temp26 temp27 temp28 temp29;
  run;

%end;

%if  &nrisk  = 3 %then %do;

  title1 "Separate Survival Model 2 - &riskind = 2";
  proc lifereg data=_dat(where=(last=1));
    model _time1*_event1(0)= &SURVpred / dist=&SURVdist2;
    ods output ParameterEstimates=SURVests20(rename=(Parameter=Effect));
    ods output ConvergenceStatus = con31 (keep = status);
  run;

  data con32;
    set con31;
    i = 1;
  run;

  * create dataset of initial Survival submodel parameters for NLMIXED;
  *Survival model parameters;
  data SURVests21; 
    length model $6 parameter $8; set SURVests20;
    format row 3.0;
    row = _n_-1;
    parameter = "c" || left(row);
    if effect="Scale" then parameter="scale2";
    if effect="Weibull Shape" then parameter="shape";
    if DF=0 then delete; *remove unestimated parameters;
    if effect="Weibull Shape" then delete; *for weib, only need scale: shape=1/scale;
    model="Surv2";
    rename LowerCL=Lower UpperCL=Upper ProbChisq=pvalue;
  run;
 
  *create Surv2 linear predictor for NLMIXED;
  data SURVests22; 
    set SURVests21; 
    if parameter NOT in("scale2", "Scale", "Shape") then do;
      linpred =  "+" || cats(parameter) || "*" || effect;
      if _N_=1 then call symput('linpredSURV2',parameter);
      else call symput('linpredSURV2',trim(resolve('&linpredSURV2'))||linpred);
    end;
    i = 1;
  run;
  %put &linpredSURV2;
  title;

  data SURVests22;
    merge SURVests22 con32;
    by i;
    drop i;
  run;

  * Get initial value for re-parameterization to fit weibull/exponential distribution with proc nlmixed; 

  %if &SurvDist2 = exponential %then %do;

    data temp23;
      set SURVests22;
      estimate = estimate * (-1);
      keep model effect parameter estimate;
    run;

  %end;

  %if &SurvDist2 = weibull %then %do;

    data temp21;
      set SURVests22;
      if parameter = 'scale2';
      keep estimate i;
      i = 1;
      rename estimate = estimate1;
    run;

    data temp22;
      set SURVests22;
      i = 1;
    run;

    data temp23;
      merge temp21 temp22;
      by i;
      if substr(parameter, 1, 1) = 'c' then estimate = -estimate / estimate1;
      keep model effect parameter estimate;
    run;

  %end;

  proc nlmixed data=_lastdat &nlmopts;
    parms / data=temp23;
    linpredSURV2  = &linpredSURV2;
    /* exponential survival;*/
    %if &SurvDist2 = exponential %then %do;
      alph  = exp(linpredSURV2);
      loglikeSURV2 = (_event1=1)*(linpredSURV2-alph*_time1)+ (_event1=0)*(-alph*_time1);
    %end;

    /* Weibull survival;*/
    %if &SurvDist2 = weibull %then %do;
      gam   = 1/scale2;
      alph  = exp(linpredSURV2);
      loglikeSURV2 = (_event1=1)*(-alph*(_time1**gam) + log(gam) + linpredSURV2 + (gam-1) * log(_time1)) + (_event1=0)*(-alph*(_time1**gam));
    %end;

    model _time1 ~ general(loglikeSURV2);

    ods output ParameterEstimates=temp24 (drop = Gradient rename =(StandardError = StdErr Probt = pvalue));
    ods output ConvergenceStatus = temp25 (keep = status rename = (status = status));
  run;

  data temp26;
    set temp25;
    i = 1;
  run;

  data temp27;
	set temp24;
	i = 1;
  run;

  data temp28;
    merge temp27 temp26;
	by i;
	drop i;
  run;

  data temp29;
    format parameter $8.;
    merge temp28 temp23 (keep = parameter effect);
    by parameter;
  run;

  title1 "Separate Survival Model 3 - &riskind = 3";
  proc lifereg data=_dat(where=(last=1));
    model _time2*_event2(0)= &SURVpred / dist=&SURVdist3;
    ods output ParameterEstimates=SURVests30(rename=(Parameter=Effect));
    ods output ConvergenceStatus = con41 (keep = status);
  run;

  data con42;
    set con41;
    i = 1;
  run;

  * create dataset of initial Survival submodel parameters for NLMIXED;
  *Survival model parameters;
  data SURVests31; 
    length model $6 parameter $8; set SURVests30;
    format row 3.0;
    row = _n_-1;
    parameter = "d" || left(row);
    if effect="Scale" then parameter="scale3";
    if effect="Weibull Shape" then parameter="shape";
    if DF=0 then delete; *remove unestimated parameters;
    if effect="Weibull Shape" then delete; *for weib, only need scale: shape=1/scale;
    model="Surv3";
    rename LowerCL=Lower UpperCL=Upper ProbChisq=pvalue;
  run; 

  *create Surv3 linear predictor for NLMIXED;
  data SURVests32; 
    set SURVests31; 
    if parameter NOT in("scale3", "Scale", "Shape") then do;
      linpred =  "+" || cats(parameter) || "*" || effect;
      if _N_=1 then call symput('linpredSURV3',parameter);
      else call symput('linpredSURV3',trim(resolve('&linpredSURV3'))||linpred);
    end;
    i = 1;
  run;
  %put &linpredSURV3;
  title;

  data SURVests32;
    merge SURVests32 con42;
    by i;
    drop i;
  run;

  * Get initial value for re-parameterization to fit weibull/exponential distribution with proc nlmixed; 

  %if &SurvDist3 = exponential %then %do;
    data temp33;
      set SURVests32;
      estimate = estimate * (-1);
      keep model effect parameter estimate;
    run;
  %end;

  %if &SurvDist3 = weibull %then %do;
    data temp31;
      set SURVests32;
      if parameter = 'scale3';
      keep estimate i;
      i = 1;
      rename estimate = estimate1;
    run;

    data temp32;
      set SURVests32;
      i = 1;
    run;

    data temp33;
      merge temp31 temp32;
      by i;
      if substr(parameter, 1, 1) = 'd' then estimate = -estimate / estimate1;
      keep model effect parameter estimate;
    run;
  %end;

  proc nlmixed data=_lastdat &nlmopts;
    parms / data=temp33;
    linpredSURV3  = &linpredSURV3;
    /*exponential survival;*/
    %if &SurvDist3 = exponential %then %do;
      alph  = exp(linpredSURV3);
      loglikeSURV3 = (_event2=1)*(linpredSURV3-alph*_time2)+ (_event2=0)*(-alph*_time2);
    %end;
    /*Weibull survival;*/
    %if &SurvDist3 = weibull %then %do;
      gam   = 1/scale3;
      alph  = exp(linpredSURV3);
      loglikeSURV3 = (_event2=1)*(-alph*(_time2**gam) + log(gam) + linpredSURV3 + (gam-1) * log(_time2)) + (_event2=0)*(-alph*(_time2**gam));
    %end;

    model _time2 ~ general(loglikeSURV3);

    ods output ParameterEstimates=temp34 (drop = Gradient rename =(StandardError = StdErr Probt = pvalue));
    ods output ConvergenceStatus = temp35 (keep = status rename = (status = status));
  run;

  data temp36;
    set temp35;
    i = 1;
  run;

  data temp37;
	set temp34;
	i = 1;
  run;

  data temp38;
	merge temp37 temp36;
	by i;
	drop i;
  run;

  data temp39;
    format parameter $8.;
    merge temp38 temp33 (keep = parameter effect);
    by parameter;
  run;

  data _SPMpar;  
    set _SPMpar temp29 (in = a) Rcorr20 Rcorr21 temp39 (in = b) Rcorr30 Rcorr31 ;
    row = _n_-1;
	if a then model = 'Surv2';
 	if b then model = 'Surv3';
    drop df tvalue alpha;
  run;

  data _SPMinit;
    set _SPMpar;
    keep model effect parameter estimate;
  run;

  proc datasets lib = work;
    delete SURVests20 SURVests21 SURVests22 SURVests30 SURVests31 SURVests32 con31 con32 con41 con42 temp21 temp22 temp23 temp24 temp25 temp26 temp27 temp28 temp29
    temp31 temp32 temp33 temp34 temp35 temp36 temp37 temp38 temp39;
  run;

%end;

*end initialization if-block;

/* When random slope is not requested --- random = 1, remove factor loading and random intercept estimate*/

%if  &random  = 1 %then %do;

  data _SPMinit;
    set _SPMinit;
    if effect = 'Intercept' and parameter = 'tau1' then parameter = 'tau0';
    if parameter in ('r11', 'r21', 'r31', 'tau1', 'tau2') then delete;
  run;

  data _SPMpar;
    set _SPMpar;
    if effect = 'Intercept' and parameter = 'tau1' then parameter = 'tau0';
    if parameter in ('r11', 'r21', 'r31', 'tau1', 'tau2') then delete;
  run;

%end;

*******************************************************************************************************
Step 4: Estimate the full shared parameter model for longitudinal and survival outcomes jointly;
*******************************************************************************************************;

%if  &nrisk  = 1 %then %do;

  title1 'Shared Parameter Model - NLMIXED 1 - &riskind = 1';
  proc nlmixed data=_dat &nlmopts;

    **************************************************************************************************
    *0) initialize parameters, all parameters not assigned starting values explicitly are assigned
    *   the default value 1
    **************************************************************************************************;
    parms / data=_SPMinit;

    ***************************************************************************************************
    *1) Compute log likelihood contribution of the LDA part 
    ***************************************************************************************************;
    /*LDA: Gaussian:*/
    linpredLDA = &linpredLDA + &linranlda; 
    resid = (_Y-linpredLDA);
    if (abs(resid) > 1.3E100) or (s2 < 1e-60) then do;
      loglikeLDA = -1e20;
    end; 
    else do;
      loglikeLDA = -0.5*(1.837876 + resid**2 / s2  + log(s2));
    end;

    ****************************************************************************************************
    *2) Compute log likelihood contribution of the survival data part when the last observation of 
        a subject is reached
    ****************************************************************************************************;
    if (last) then do; 
      /*Note: the "last" var is created earlier datastep;*/
      /* Drop out Likelihood Specification */
      linpredSURV1  = &linpredSURV1 + &linransurv1; 

      /*exponential survival;*/
      %if &SurvDist1 = exponential %then %do;
        alph  = exp(linpredsurv1);
        loglikeSURV1 = (_event=1)*(linpredsurv1-alph*_time)+ (_event=0)*(-alph*_time);
      %end;

      /*Weibull survival;*/
      %if &SurvDist1 = weibull %then %do;
        gam   = 1/scale1;
        alph  = exp(linpredsurv1);
        loglikeSURV1 = (_event=1)*(-alph*(_time**gam) + log(gam) + linpredsurv1 + (gam-1) * log(_time)) + (_event=0)*(-alph*(_time**gam));
      %end;

    end; 
    else do;
      loglikeSURV1=0;
    end;

    ****************************************************************************************************
    *3) Compute log likelihood                          
    ****************************************************************************************************;
    model last ~ general(loglikeLDA + loglikeSURV1);

    random &ranspe subject=_id;
 
	/* Compute the reparameterized loading factor estimates using ESTIMATE statement of NLMIXED*/
	%if  &random  = 1 %then %do;
      estimate 'r10' sqrt(tau0) * r10;
    %end;

    %if  &random  = 2 %then %do;
      estimate 'r10' sqrt(tau0) * r10;
      estimate 'r11' sqrt(tau1) * r11;
    %end;

    ods output ParameterEstimates=SPMests;
    ods output FitStatistics = fit1; 
    ods output ConvergenceStatus = con51 (keep = status rename = (status = joint_status));
    ods output AdditionalEstimates = est1 (rename = (label = parameter));
  run;

%end;

%if  &nrisk  = 2 %then %do;

  title1 'Shared Parameter Model - NLMIXED 2 - &riskind = 2';
  proc nlmixed data=_dat &nlmopts;

    **************************************************************************************************
    *0) initialize parameters, all parameters not assigned starting values explicitly are assigned
    *   the default value 1
    **************************************************************************************************;
    parms / data=_SPMinit;

    ***************************************************************************************************
    *1) Compute log likelihood contribution of the LDA part 
    ***************************************************************************************************;
    /*LDA: Gaussian:*/
    linpredLDA = &linpredLDA + &linranlda; 
    resid = (_Y-linpredLDA);
    if (abs(resid) > 1.3E100) or (s2 < 1e-60) then do;
      loglikeLDA = -1e20;
    end; 
    else do;
      loglikeLDA = -0.5*(1.837876 + resid**2 / s2  + log(s2));
    end;

    ****************************************************************************************************
    *2) Compute log likelihood contribution of the survival data part when the last observation of 
        a subject is reached
    ****************************************************************************************************;
    if (last) then do; 
      /*Note: the "last" var is created earlier datastep;*/

      /*First competing Risk Likelihood Specification */
      linpredSURV1  = &linpredSURV1 + &linransurv1;

      /*exponential survival;*/
      %if &SurvDist1 = exponential %then %do;
        alph  = exp(linpredsurv1);
        loglikeSURV1 = (_event=1)*(linpredsurv1-alph*_time)+ (_event=0)*(-alph*_time);
      %end;

      /*Weibull survival;*/
      %if &SurvDist1 = weibull %then %do;
        gam   = 1/scale1;
        alph  = exp(linpredsurv1);
        loglikeSURV1 = (_event=1)*(-alph*(_time**gam) + log(gam) + linpredsurv1 + (gam-1) * log(_time)) + (_event=0)*(-alph*(_time**gam));
      %end;

      /*Second competing Risk Likelihood Specification */
      linpredsurv2  = &linpredsurv2 + &linransurv2; *need to macro-ize;

      /*exponential survival;*/
      %if &SurvDist2 = exponential %then %do;
        alph  = exp(linpredSURV2);
        loglikeSURV2 = (_event1=1)*(linpredSURV2-alph*_time1)+ (_event1=0)*(-alph*_time1);
      %end;

      /*Weibull survival;*/
      %if &SurvDist2 = weibull %then %do;
        gam   = 1/scale2;
        alph  = exp(linpredSURV2);
        loglikeSURV2 = (_event1=1)*(-alph*(_time1**gam) + log(gam) + linpredSURV2 + (gam-1) * log(_time1)) + (_event1=0)*(-alph*(_time1**gam));
      %end;

    end; 
    else do;
      loglikeSURV1=0;
      loglikeSURV2=0;
    end;

    ****************************************************************************************************
    *3) Compute log likelihood                          
    ****************************************************************************************************;
    model last ~ general(loglikeLDA + loglikeSURV1 + loglikeSURV2);

    random &ranspe subject=_id;

	/* Compute the reparameterized loading factor estimates using ESTIMATE statement of NLMIXED*/

	%if  &random  = 1 %then %do;
      estimate 'r10' sqrt(tau0) * r10;
      estimate 'r20' sqrt(tau0) * r20;
    %end;

	%if  &random  = 2 %then %do;
      estimate 'r10' sqrt(tau0) * r10;
      estimate 'r20' sqrt(tau0) * r20;
      estimate 'r11' sqrt(tau1) * r11;
      estimate 'r21' sqrt(tau1) * r21;
    %end;

    ods output ParameterEstimates=SPMests;
    ods output FitStatistics = fit1; 
    ods output ConvergenceStatus = con51 (keep = status rename = (status = joint_status));
    ods output AdditionalEstimates = est1 (rename = (label = parameter));
  run;

  %end;

%if  &nrisk  = 3 %then %do;

  title1 'Shared Parameter Model - NLMIXED 1 - &riskind = 3';
  proc nlmixed data=_dat &nlmopts;

    **************************************************************************************************
    *0) initialize parameters, all parameters not assigned starting values explicitly are assigned
    *   the default value 1
    **************************************************************************************************;
    parms / data=_SPMinit;

    ***************************************************************************************************
    *1) Compute log likelihood contribution of the LDA part 
    ***************************************************************************************************;
    /*LDA: Gaussian:*/
    linpredLDA = &linpredLDA + &linranlda; 
    resid = (_Y-linpredLDA);
    if (abs(resid) > 1.3E100) or (s2 < 1e-60) then do;
      loglikeLDA = -1e20;
    end; 
    else do;
      loglikeLDA = -0.5*(1.837876 + resid**2 / s2  + log(s2));
    end;

    ****************************************************************************************************
    *2) Compute log likelihood contribution of the survival data part when the last observation of 
        a subject is reached
    ****************************************************************************************************;
      if (last) then do; 
      /*Note: the "last" var is created earlier datastep;*/

      /* First competing Risk Likelihood Specification */
      linpredSURV1  = &linpredSURV1 + &linransurv1; *need to macro-ize;

      /*exponential survival;*/
      %if &SurvDist1 = exponential %then %do;
        alph  = exp(linpredsurv1);
        loglikeSURV1 = (_event=1)*(linpredsurv1-alph*_time)+ (_event=0)*(-alph*_time);
      %end;

      /*Weibull survival;*/
      %if &SurvDist1 = weibull %then %do;
        gam   = 1/scale1;
        alph  = exp(linpredsurv1);
        loglikeSURV1 = (_event=1)*(-alph*(_time**gam) + log(gam) + linpredsurv1 + (gam-1) * log(_time)) + (_event=0)*(-alph*(_time**gam));
      %end;

      /* Second competing Risk Likelihood Specification */
      linpredsurv2  = &linpredsurv2 + &linransurv2; *need to macro-ize;

      /*exponential survival;*/
      %if &SurvDist2 = exponential %then %do;
        alph  = exp(linpredSURV2);
        loglikeSURV2 = (_event1=1)*(linpredSURV2-alph*_time1)+ (_event1=0)*(-alph*_time1);
      %end;

      /*Weibull survival;*/
      %if &SurvDist2 = weibull %then %do;
        gam   = 1/scale2;
        alph  = exp(linpredSURV2);
        loglikeSURV2 = (_event1=1)*(-alph*(_time1**gam) + log(gam) + linpredSURV2 + (gam-1) * log(_time1)) + (_event1=0)*(-alph*(_time1**gam));
      %end;

      /* Third competing Risk Likelihood Specification */
      linpredsurv3  = &linpredsurv3 + &linransurv3; *need to macro-ize;

      /*exponential survival;*/
      %if &SurvDist3 = exponential %then %do;
        alph  = exp(linpredSURV3);
        loglikeSURV3 = (_event2=1)*(linpredSURV3-alph*_time2)+ (_event2=0)*(-alph*_time2);
      %end;

      /*Weibull survival;*/
      %if &SurvDist3 = weibull %then %do;
        gam   = 1/scale3;
        alph  = exp(linpredSURV3);
        loglikeSURV3 = (_event2=1)*(-alph*(_time2**gam) + log(gam) + linpredSURV3 + (gam-1) * log(_time2)) + (_event2=0)*(-alph*(_time2**gam));
      %end;

    end; 
    else do;
      loglikeSURV1=0;
      loglikeSURV2=0;
      loglikeSURV3=0;
    end;

    ****************************************************************************************************
    *3) Compute log likelihood                          
    ****************************************************************************************************;
    model last ~ general(loglikeLDA + loglikeSURV1 + loglikeSURV2 + loglikeSURV3);

    random &ranspe subject=_id;

	/* Compute the reparameterized loading factor estimates using ESTIMATE statement of NLMIXED*/

	%if  &random  = 1 %then %do;
      estimate 'r10' sqrt(tau0) * r10;
      estimate 'r20' sqrt(tau0) * r20;
      estimate 'r30' sqrt(tau0) * r30;
    %end;

	%if  &random  = 2 %then %do;
      estimate 'r10' sqrt(tau0) * r10;
      estimate 'r20' sqrt(tau0) * r20;
      estimate 'r30' sqrt(tau0) * r30;
      estimate 'r11' sqrt(tau1) * r11;
      estimate 'r21' sqrt(tau1) * r21;
      estimate 'r31' sqrt(tau1) * r31;
    %end;


    ods output ParameterEstimates=SPMests;
    ods output FitStatistics = fit1; 
    ods output ConvergenceStatus = con51 (keep = status rename = (status = joint_status));
    ods output AdditionalEstimates = est1 (rename = (label = parameter));
  run;

%end;

proc sort data = est1;
  by parameter;
run;

proc sort data = SPMests;
  by parameter;
run;

data SPMests;
  merge SPMests est1;
  by parameter;
run;

data con52;
  set con51;
  i = 1;
run;

data SPMests;
  set SPMests;
  i = 1;
run;

data SPMests;
  merge SPMests con52;
  by i;
  drop i;
run;

*******************************************************************************************************
Step 5: Organize the resulting estimates and output data sets;
*******************************************************************************************************;

proc transpose data = fit1 out = fit2 (rename = (AIC__smaller_is_better_ = aic BIC__smaller_is_better_ = bic) drop = AICC__smaller_is_better_  _name_);
  id descr;
  var value;
run;

data fit2;
  merge fit2 con52;
  drop i;
run;

proc sort data=_SPMpar; 
  by Parameter; 
run;
 
proc sort data=SPMests; 
  by Parameter; 
run;
 
data SPMests; 
  merge _SPMpar(keep=model Parameter effect row) SPMests; 
  by Parameter; 
run;

data Estimates1; 
  merge _SPMpar (rename=(estimate=SeparateEstimate
                        stderr=SeparateStdErr
                        Lower=SeparateLowerCL
                        Upper=SeparateUpperCL
                        pvalue=SeparatePvalue)) 
   SPMests (rename=(estimate=JointEstimate
                    StandardError=JointStdErr
                    Lower=JointLowerCL
                    Upper=JointUpperCL
                    Probt=JointPvalue)); 
   by Parameter; 
   format SeparateStatus JointStatus $15.;
   if status = 0 then SeparateStatus = 'Convergence';
   else if status = 1 then SeparateStatus = 'Non_Converg';
   if joint_status = 0 then JointStatus = 'Convergence';
   else if joint_status = 1 then JointStatus = 'Non_Converg';
   drop df tvalue alpha gradient linpred;
run;

proc sort data=Estimates1 out=Estimates(drop=row); 
  by row; 
run;

proc sort data=_SPMpar; 
  by row; 
run; 

proc sort data=SPMests; 
  by row; 
run; 

data &outaic;
  set fit2;
run;

data &outpar;
  set estimates;
run;

title;

proc datasets lib = work;
  delete beta beta1 cov cov1 estimates1 estimates fit1 fit2 ldaests rcorr rcorr10 rcorr11 rcorr20 rcorr21 rcorr30 rcorr31 spmests
  _dat _spminit _spmpar con51 con52 _lastdat est1;
run;

%mend;

** SPM is complied Macro to call macro SPM1 using all possible combinations of baseline parametric hazard functions, SURVdist1, SURVdist2 or SURVdist3. e.g. when nrisk = 2, 
SPM1 will be called four times with the following combinations, 

1. SURVdist1   = exponential and SURVdist2   = exponential,
2. SURVdist1   = exponential and SURVdist2   = weibull,
3. SURVdist1   = weibull and SURVdist2   = exponential,
4. SURVdist1   = weibull and SURVdist2   = weibull;

%macro SPM(data        = , 
           id          = ,
           nrisk       = ,
           LDAoutcome  = ,
           LDApred     = , 
           LDAtime     = , 
           riskind     = ,
           SURVtime    = ,
           SURVpred    = ,
           Random      = ,
           nlmopts     = ,
           printall    = T,
           allfit      = ,
           allpar      = ,
           leastaicpar =
           );

* nrisk = 1, call SPM1 twice, SURVdist1 = exponential or weibull;

%if  &nrisk  = 1 %then %do;

  %let SURVdist1   = exponential;

  %SPM1(data       = &data, 
        id         = &id,
        nrisk      = &nrisk,
        LDAoutcome = &LDAoutcome,
        LDApred    = &LDApred, 
        LDAtime    = &LDAtime, 
        riskind    = &riskind,
        SURVtime   = &SURVtime,
        SURVpred   = &SURVpred,
        Random     = &Random,
        nlmopts    = &nlmopts,
        outaic     = yy01,
        outpar     = zz01
        );

  %let SURVdist1   = weibull;

  %SPM1(data       = &data, 
        id         = &id,
        nrisk      = &nrisk,
        LDAoutcome = &LDAoutcome,
        LDApred    = &LDApred, 
        LDAtime    = &LDAtime, 
        riskind    = &riskind,
        SURVtime   = &SURVtime,
        SURVpred   = &SURVpred,
        Random     = &Random,
        nlmopts    = &nlmopts,
        outaic     = yy02,
        outpar     = zz02
        );

  data final1;
    set yy01 (in = a) yy02 (in = b);
    if a then do;
      surv = 'Exponential';
    end;
    if b then do;
      surv = 'Weibull';
    end;
    ind = _n_;
    n = 1;
    if joint_status = 0 then JointStatus = 'Convergence';
    else if joint_status = 1 then JointStatus = 'Non_Converg';
    drop joint_status;
  run;

  data &allfit;
    set final1;
    drop ind n;
  run;

  data final1;
    set final1;
    if JointStatus = 'Convergence';
  run;

  proc sort data = final1;
    by aic;
  run;

  data temp1;
    set final1;
    by n;
    if first.n;
    keep ind;
  run;

  data final2;
    set zz01 (in = a) zz02 (in = b);
    if a then do;
      ind = 1;
      surv = 'Exponential';
    end;
    if b then do;
      ind = 2;
      surv = 'Weibull';
    end;
  run;

  data &allpar;
    set final2;
    drop ind status joint_status;
  run;

  data &leastaicpar;
    merge final2 temp1 (in = a);
    by ind;
    if a;
    drop status joint_status;
  run;

  proc datasets lib = work;
    delete final1 final2 temp1 yy01 yy02 zz01 zz02;
  run;

%end;

* nrisk = 2, call SPM1 four times, SURVdist1 = exponential or Weibull, SURVdist2 = exponential or weibull;

%if  &nrisk  = 2 %then %do;

  %let SURVdist1   = exponential;
  %let SURVdist2   = exponential;

  %SPM1(data       = &data, 
        id         = &id,
        nrisk      = &nrisk,
        LDAoutcome = &LDAoutcome,
        LDApred    = &LDApred, 
        LDAtime    = &LDAtime, 
        riskind    = &riskind,
        SURVtime   = &SURVtime,
        SURVpred   = &SURVpred,
        Random     = &Random,
        nlmopts    = &nlmopts,
        outaic     = yy01,
        outpar     = zz01
        );

  %let SURVdist1   = exponential;
  %let SURVdist2   = weibull;

  %SPM1(data       = &data, 
        id         = &id,
        nrisk      = &nrisk,
        LDAoutcome = &LDAoutcome,
        LDApred    = &LDApred, 
        LDAtime    = &LDAtime, 
        riskind    = &riskind,
        SURVtime   = &SURVtime,
        SURVpred   = &SURVpred,
        Random     = &Random,
        nlmopts    = &nlmopts,
        outaic     = yy02,
        outpar     = zz02
        );

  %let SURVdist1   = weibull;
  %let SURVdist2   = exponential;

  %SPM1(data       = &data, 
        id         = &id,
        nrisk      = &nrisk,
        LDAoutcome = &LDAoutcome,
        LDApred    = &LDApred, 
        LDAtime    = &LDAtime, 
        riskind    = &riskind,
        SURVtime   = &SURVtime,
        SURVpred   = &SURVpred,
        Random     = &Random,
        nlmopts    = &nlmopts,
        outaic     = yy03,
        outpar     = zz03
        );

  %let SURVdist1   = weibull;
  %let SURVdist2   = weibull;

  %SPM1(data       = &data, 
        id         = &id,
        nrisk      = &nrisk,
        LDAoutcome = &LDAoutcome,
        LDApred    = &LDApred, 
        LDAtime    = &LDAtime, 
        riskind    = &riskind,
        SURVtime   = &SURVtime,
        SURVpred   = &SURVpred,
        Random     = &Random,
        nlmopts    = &nlmopts,
        outaic     = yy04,
        outpar     = zz04
        );

  data final1;
    set yy01 (in = a) yy02 (in = b) yy03 (in = c) yy04 (in = d);
    if a then do;
      surv1 = 'Exponential';
      surv2 = 'Exponential';
    end;
    if b then do;
      surv1 = 'Exponential';
      surv2 = 'Weibull';
    end;
    if c then do;
      surv1 = 'Weibull';
      surv2 = 'Exponential';
    end;
    if d then do;
      surv1 = 'Weibull';
      surv2 = 'Weibull';
    end;
    ind = _n_;
    n = 1;
    if joint_status = 0 then JointStatus = 'Convergence';
    else if joint_status = 1 then JointStatus = 'Non_Converg';
    drop joint_status;
  run;

  data &allfit;
    set final1;
    drop ind n;
  run;

  data final1;
    set final1;
    if JointStatus = 'Convergence';
  run;

  proc sort data = final1;
    by aic;
  run;

  data temp1;
    set final1;
    by n;
    if first.n;
    keep ind;
  run;

  data final2;
    set zz01 (in = a) zz02 (in = b) zz03 (in = c) zz04 (in = d);
    if a then do;
      ind = 1;
      surv1 = 'Exponential';
      surv2 = 'Exponential';
    end;
    if b then do;
      ind = 2;
      surv1 = 'Exponential';
      surv2 = 'Weibull';
    end;
    if c then do;
      ind = 3;
      surv1 = 'Weibull';
      surv2 = 'Exponential';
    end;
    if d then do;
      ind = 4;
      surv1 = 'Weibull';
      surv2 = 'Weibull';
    end;
  run;

  data &allpar;
    set final2;
    drop ind status joint_status;
  run;

  data &leastaicpar;
    merge final2 temp1 (in = a);
    by ind;
    if a;
    drop status joint_status;
  run;

  proc datasets lib = work;
    delete final1 final2 temp1 yy01 yy02 yy03 yy04 zz01 zz02 zz03 zz04;
  run;

%end;

* nrisk = 3, call SPM1 eight times, SURVdist1 = exponential or Weibull, SURVdist2 = exponential or weibull, SURVdist3 = exponential or weibull;

%if  &nrisk  = 3 %then %do;

  %let SURVdist1   = exponential;
  %let SURVdist2   = exponential;
  %let SURVdist3   = exponential;

  %SPM1(data       = &data, 
        id         = &id,
        nrisk      = &nrisk,
        LDAoutcome = &LDAoutcome,
        LDApred    = &LDApred, 
        LDAtime    = &LDAtime, 
        riskind    = &riskind,
        SURVtime   = &SURVtime,
        SURVpred   = &SURVpred,
        Random     = &Random,
        nlmopts    = &nlmopts,
        outaic     = yy01,
        outpar     = zz01
        );

  %let SURVdist1   = exponential;
  %let SURVdist2   = weibull;
  %let SURVdist3   = exponential;

  %SPM1(data       = &data, 
        id         = &id,
        nrisk      = &nrisk,
        LDAoutcome = &LDAoutcome,
        LDApred    = &LDApred, 
        LDAtime    = &LDAtime, 
        riskind    = &riskind,
        SURVtime   = &SURVtime,
        SURVpred   = &SURVpred,
        Random     = &Random,
        nlmopts    = &nlmopts,
        outaic     = yy02,
        outpar     = zz02
        );

  %let SURVdist1   = weibull;
  %let SURVdist2   = exponential;
  %let SURVdist3   = exponential;

  %SPM1(data       = &data, 
        id         = &id,
        nrisk      = &nrisk,
        LDAoutcome = &LDAoutcome,
        LDApred    = &LDApred, 
        LDAtime    = &LDAtime, 
        riskind    = &riskind,
        SURVtime   = &SURVtime,
        SURVpred   = &SURVpred,
        Random     = &Random,
        nlmopts    = &nlmopts,
        outaic     = yy03,
        outpar     = zz03
        );

  %let SURVdist1   = weibull;
  %let SURVdist2   = weibull;
  %let SURVdist3   = exponential;

  %SPM1(data       = &data, 
        id         = &id,
        nrisk      = &nrisk,
        LDAoutcome = &LDAoutcome,
        LDApred    = &LDApred, 
        LDAtime    = &LDAtime, 
        riskind    = &riskind,
        SURVtime   = &SURVtime,
        SURVpred   = &SURVpred,
        Random     = &Random,
        nlmopts    = &nlmopts,
        outaic     = yy04,
        outpar     = zz04
        );

  %let SURVdist1   = exponential;
  %let SURVdist2   = exponential;
  %let SURVdist3   = weibull;

  %SPM1(data       = &data, 
        id         = &id,
        nrisk      = &nrisk,
        LDAoutcome = &LDAoutcome,
        LDApred    = &LDApred, 
        LDAtime    = &LDAtime, 
        riskind    = &riskind,
        SURVtime   = &SURVtime,
        SURVpred   = &SURVpred,
        Random     = &Random,
        nlmopts    = &nlmopts,
        outaic     = yy05,
        outpar     = zz05
        );

  %let SURVdist1   = exponential;
  %let SURVdist2   = weibull;
  %let SURVdist3   = weibull;

  %SPM1(data       = &data, 
        id         = &id,
        nrisk      = &nrisk,
        LDAoutcome = &LDAoutcome,
        LDApred    = &LDApred, 
        LDAtime    = &LDAtime, 
        riskind    = &riskind,
        SURVtime   = &SURVtime,
        SURVpred   = &SURVpred,
        Random     = &Random,
        nlmopts    = &nlmopts,
        outaic     = yy06,
        outpar     = zz06
        );

  %let SURVdist1   = weibull;
  %let SURVdist2   = exponential;
  %let SURVdist3   = weibull;

  %SPM1(data       = &data, 
        id         = &id,
        nrisk      = &nrisk,
        LDAoutcome = &LDAoutcome,
        LDApred    = &LDApred, 
        LDAtime    = &LDAtime, 
        riskind    = &riskind,
        SURVtime   = &SURVtime,
        SURVpred   = &SURVpred,
        Random     = &Random,
        nlmopts    = &nlmopts,
        outaic     = yy07,
        outpar     = zz07
        );

  %let SURVdist1   = weibull;
  %let SURVdist2   = weibull;
  %let SURVdist3   = weibull;

  %SPM1(data       = &data, 
        id         = &id,
        nrisk      = &nrisk,
        LDAoutcome = &LDAoutcome,
        LDApred    = &LDApred, 
        LDAtime    = &LDAtime, 
        riskind    = &riskind,
        SURVtime   = &SURVtime,
        SURVpred   = &SURVpred,
        Random     = &Random,
        nlmopts    = &nlmopts,
        outaic     = yy08,
        outpar     = zz08
        );

  data final1;
    set yy01 (in = a) yy02 (in = b) yy03 (in = c) yy04 (in = d) yy05 (in = e) yy06 (in = f) yy07 (in = g) yy08 (in = h);
    if a then do;
      surv1 = 'Exponential';
      surv2 = 'Exponential';
      surv3 = 'Exponential';
    end;
    if b then do;
      surv1 = 'Exponential';
      surv2 = 'Weibull';
      surv3 = 'Exponential';
    end;
    if c then do;
      surv1 = 'Weibull';
      surv2 = 'Exponential';
      surv3 = 'Exponential';
    end;
    if d then do;
      surv1 = 'Weibull';
      surv2 = 'Weibull';
      surv3 = 'Exponential';
    end;
    if e then do;
      surv1 = 'Exponential';
      surv2 = 'Exponential';
      surv3 = 'Weibull';
    end;
    if f then do;
      surv1 = 'Exponential';
      surv2 = 'Weibull';
      surv3 = 'Weibull';
    end;
    if g then do;
      surv1 = 'Weibull';
      surv2 = 'Exponential';
      surv3 = 'Weibull';
    end;
    if h then do;
      surv1 = 'Weibull';
      surv2 = 'Weibull';
      surv3 = 'Weibull';
    end;
    ind = _n_;
    n = 1;
    if joint_status = 0 then JointStatus = 'Convergence';
    else if joint_status = 1 then JointStatus = 'Non_Converg';
    drop joint_status;
  run;

  data &allfit;
    set final1;
    drop ind n;
  run;

  data final1;
    set final1;
    if JointStatus = 'Convergence';
  run;

  proc sort data = final1;
    by aic;
  run;

  data temp1;
    set final1;
    by n;
    if first.n;
    keep ind;
  run;

  data final2;
    set zz01 (in = a) zz02 (in = b) zz03 (in = c) zz04 (in = d) zz05 (in = e) zz06 (in = f) zz07 (in = g) zz08 (in = h);
    if a then do;
      ind = 1;
      surv1 = 'Exponential';
      surv2 = 'Exponential';
      surv3 = 'Exponential';
    end;
    if b then do;
      ind = 2;
      surv1 = 'Exponential';
      surv2 = 'Weibull';
      surv3 = 'Exponential';
    end;
    if c then do;
      ind = 3;
      surv1 = 'Weibull';
      surv2 = 'Exponential';
      surv3 = 'Exponential';
    end;
    if d then do;
      ind = 4;
      surv1 = 'Weibull';
      surv2 = 'Weibull';
      surv3 = 'Exponential';
    end;
    if e then do;
      ind = 5;
      surv1 = 'Exponential';
      surv2 = 'Exponential';
      surv3 = 'Weibull';
    end;
    if f then do;
      ind = 6;
      surv1 = 'Exponential';
      surv2 = 'Weibull';
      surv3 = 'Weibull';
    end;
    if g then do;
      ind = 7;
      surv1 = 'Weibull';
      surv2 = 'Exponential';
      surv3 = 'Weibull';
    end;
    if h then do;
      ind = 8;
      surv1 = 'Weibull';
      surv2 = 'Weibull';
      surv3 = 'Weibull';
    end;
  run;

  data &allpar;
    set final2;
    drop ind status joint_status;
  run;

  data &leastaicpar;
    merge final2 temp1 (in = a);
    by ind;
    if a;
    drop status joint_status;
  run;

  proc datasets lib = work;
    delete final1 final2 temp1 yy01 yy02 yy03 yy04 yy05 yy06 yy07 yy08 zz01 zz02 zz03 zz04 zz05 zz06 zz07 zz08;
  run;

%end;

%if &printall=T %then %do;

  title1 'Estimates and Standard Errors from Separate (MAR) and SPM (MNAR) Analysis';
  proc print data=&leastaicpar label noobs;
    label Effect = "Effect"
    SeparateEstimate = 'MLE in separate models'
    SeparateStdErr   = 'Asymptotic Std. Error of MLE in separate models'
    SeparatePvalue   = 'Pvalue in separate models'
    SeparateStatus = 'Converge status in separate models'
    JointEstimate    = 'MLE in joint model'
    JointStdErr      = 'Asymptotic Std. Error of MLE in joint model'
    JointPvalue      = 'Pvalue in joint model'
    JointStatus = 'Converge status in joint model';
    var model effect Parameter
        SeparateEstimate SeparateStdErr SeparatePvalue SeparateStatus
        JointEstimate    JointStdErr    JointPvalue JointStatus;
  run;

  title1 'Separate (MAR) and SPM (MNAR) Analysis Results';
  proc print data=&leastaicpar noobs; 
    var  model effect Parameter 
         SeparateEstimate SeparateStdErr SeparateLowerCL SeparateUpperCL SeparatePvalue SeparateStatus; 
  run;

  proc print data=&leastaicpar noobs; 
    var  model effect Parameter 
         JointEstimate JointStdErr JointLowerCL JointUpperCL JointPvalue JointStatus; 
  run;

%end;  *-end printall=T option;

title;

%mend;

proc import datafile = 'C:\Users\wwang\Box Sync\manuscripts\meth001griswoldspm\4-docs\paper\SPM Macro Paper\Github Upload SPM Macro\example.xlsx' out = sim replace;
run;

%SPM( data         = sim, 
      id           = idnum,
      nrisk        = 2,
      LDAoutcome   = total_z,
      LDApred      = time agev2c timeagev2c hypert25 hypertime educ1 educ2, 
      LDAtime      = time, 
      riskind      = riskind,
      SURVtime     = _t,
      SURVpred     = agev2c educ1 educ2 hypert25,
      Random       = 2,
      nlmopts      = maxfunc = 50000 maxiter = 5000,
      printall     = T,
      allfit       = gg51,
      allpar       = gg52,
	  leastaicpar  = gg53
  );
