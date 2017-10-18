function [fvalin,fval] = EulerEquation(x)
% EulerEquation
% Calculates the first order conditions

global exstate endstate1 xgrid
global FL FH
global gamma1 gamma2 beta
global div1 end1
global passet11 passet12 passet13 passet14 passet15 passet16 passet17 passet18 
global hh1house1 hh1house2 hh1house3 hh1house4 hh1house5 hh1house6 hh1house7 hh1house8
global phouse1 phouse2 phouse3 phouse4 phouse5 phouse6 phouse7 phouse8
global Value101 Value102 Value103 Value104 Value105 Value106 Value107 Value108
global Value201 Value202 Value203 Value204 Value205 Value206 Value207 Value208
global ch101 ch102 ch103 ch104 ch105 ch106 ch107 ch108
global ch201 ch202 ch203 ch204 ch205 ch206 ch207 ch208
global ev1 ev2 psi1 psi2 growth
global margin1 margin2 alpha
global pix1 pix2

fval=zeros(22,1);
fvalin=[];
%---------------------------------------------------------
% first set of variables: guesses of future wealth share
%---------------------------------------------------------
hh1wshareguess=zeros(8,1);
for ii=1:8
    hh1wshareguess(ii)=x(ii);
end
%-------------------------------------------------------
% guesses for portfolio choices
%-------------------------------------------------------
hh1asset1guess      =x(9);
hh2asset1guess      =1-x(9);
hh1bondguess        =x(10);
hh2bondguess        =-x(10);
hh1houseguess       =x(11);
hh2houseguess       =1-x(11);



%------------------------------------------------------
% guesses for consumption
%------------------------------------------------------
ch1=exp(x(12));
ch2=exp(x(13));
%------------------------------------------------------
% guesses for lagrange multiplier
%------------------------------------------------------
lag1collpos         =max(0,x(14))^2;
lag1collneg         =max(0,-x(14))^2;
lag2collpos         =max(0,x(15))^2;
lag2collneg         =max(0,-x(15))^2;

lag1assetshortpos   =max(0,x(16))^2;
lag1assetshortneg  =max(0,-x(16))^2;
lag2assetshortpos   =max(0,x(17))^2;
lag2assetshortneg  =max(0,-x(17))^2;
lag1houseshortpos  =max(0,x(18))^2;
lag1houseshortneg  =max(0,-x(18))^2;
lag2houseshortpos  =max(0,x(19))^2;
lag2houseshortneg  =max(0,-x(19))^2;
%-------------------------------------------------------
% guesses for prices
%------------------------------------------------------
q1g=exp(x(20));
qbg=exp(x(21));
qhg=exp(x(22));
%---------------------------------------------------------
% Interpolate future prices and consumption
%---------------------------------------------------------
ch1t=zeros(8,1);
ch2t=zeros(8,1);
ch1tilde=zeros(8,1);
ch2tilde=zeros(8,1);
hh1houset=zeros(8,1);
phouset=zeros(8,1);
passet1t=zeros(8,1);
Value1t=zeros(8,1);
Value2t=zeros(8,1);

ii=1;
passet1t(ii)        =exp(ppval(passet11,hh1wshareguess(ii)));
hh1houset(ii)       =ppval(hh1house1,hh1wshareguess(ii));
phouset(ii)         =ppval(phouse1,hh1wshareguess(ii));
%pbondt(ii)        =ppval(pbond1,hh1wshareguess(ii));
Value1t(ii)     =ppval(Value101,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value201,hh1wshareguess(ii));
ch1t(ii)        =ppval(ch101,hh1wshareguess(ii));
ch2t(ii)        =ppval(ch201,hh1wshareguess(ii));

ii=2;
passet1t(ii)    =exp(ppval(passet12,hh1wshareguess(ii)));
hh1houset(ii)       =ppval(hh1house2,hh1wshareguess(ii));
phouset(ii)         =ppval(phouse2,hh1wshareguess(ii));
%pbondt(ii)        =ppval(pbond2,hh1wshareguess(ii));
Value1t(ii)     =ppval(Value102,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value202,hh1wshareguess(ii));
ch1t(ii)        =ppval(ch102,hh1wshareguess(ii));
ch2t(ii)        =ppval(ch202,hh1wshareguess(ii));

ii=3;
passet1t(ii)    =exp(ppval(passet13,hh1wshareguess(ii)));
hh1houset(ii)       =ppval(hh1house3,hh1wshareguess(ii));
phouset(ii)         =ppval(phouse3,hh1wshareguess(ii));
%pbondt(ii)        =ppval(pbond3,hh1wshareguess(ii));
Value1t(ii)     =ppval(Value103,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value203,hh1wshareguess(ii));
ch1t(ii)        =ppval(ch103,hh1wshareguess(ii));
ch2t(ii)        =ppval(ch203,hh1wshareguess(ii));

ii=4;
passet1t(ii)    =exp(ppval(passet14,hh1wshareguess(ii)));
hh1houset(ii)       =ppval(hh1house4,hh1wshareguess(ii));
phouset(ii)         =ppval(phouse4,hh1wshareguess(ii));
%pbondt(ii)        =ppval(pbond4,hh1wshareguess(ii));
Value1t(ii)     =ppval(Value104,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value204,hh1wshareguess(ii));
ch1t(ii)        =ppval(ch104,hh1wshareguess(ii));
ch2t(ii)        =ppval(ch204,hh1wshareguess(ii));

ii=5;
passet1t(ii)    =exp(ppval(passet15,hh1wshareguess(ii)));
hh1houset(ii)       =ppval(hh1house5,hh1wshareguess(ii));
phouset(ii)         =ppval(phouse5,hh1wshareguess(ii));
%pbondt(ii)        =ppval(pbond5,hh1wshareguess(ii));
Value1t(ii)     =ppval(Value105,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value205,hh1wshareguess(ii));
ch1t(ii)        =ppval(ch105,hh1wshareguess(ii));
ch2t(ii)        =ppval(ch205,hh1wshareguess(ii));

ii=6;
passet1t(ii)    =exp(ppval(passet16,hh1wshareguess(ii)));
hh1houset(ii)       =ppval(hh1house6,hh1wshareguess(ii));
phouset(ii)         =ppval(phouse6,hh1wshareguess(ii));
%pbondt(ii)        =ppval(pbond6,hh1wshareguess(ii));
Value1t(ii)     =ppval(Value106,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value206,hh1wshareguess(ii));
ch1t(ii)        =ppval(ch106,hh1wshareguess(ii));
ch2t(ii)        =ppval(ch206,hh1wshareguess(ii));

ii=7;
passet1t(ii)    =exp(ppval(passet17,hh1wshareguess(ii)));
hh1houset(ii)       =ppval(hh1house7,hh1wshareguess(ii));
phouset(ii)         =ppval(phouse7,hh1wshareguess(ii));
%pbondt(ii)        =ppval(pbond7,hh1wshareguess(ii));
Value1t(ii)     =ppval(Value107,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value207,hh1wshareguess(ii));
ch1t(ii)        =ppval(ch107,hh1wshareguess(ii));
ch2t(ii)        =ppval(ch207,hh1wshareguess(ii));

ii=8;
passet1t(ii)    =exp(ppval(passet18,hh1wshareguess(ii)));
hh1houset(ii)       =ppval(hh1house8,hh1wshareguess(ii));
phouset(ii)         =ppval(phouse8,hh1wshareguess(ii));
%pbondt(ii)        =ppval(pbond8,hh1wshareguess(ii));
Value1t(ii)     =ppval(Value108,hh1wshareguess(ii));
Value2t(ii)     =ppval(Value208,hh1wshareguess(ii));
ch1t(ii)        =ppval(ch108,hh1wshareguess(ii));
ch2t(ii)        =ppval(ch208,hh1wshareguess(ii));

ch1t=exp(ch1t);
ch2t=exp(ch2t);
hh2houset=1-hh1houset;
phouset=exp(phouset);
pmin=min(phouset);
pamin=min(passet1t+0.15);

for ii=1:8
%    ch1t(ii)    =end1(ii)+hh1asset1guess*(passet1t(ii)+div1(ii))+hh1bondguess/growth(ii);
%    ch1t(ii)    =ch1t(ii)-hh1asset1t(ii)*passet1t(ii)-hh1bondt(ii)*pbondt(ii);
%    ch2t(ii)    =end2(ii)+hh2asset1guess*(passet1t(ii)+div1(ii))+hh2bondguess/growth(ii);
%    ch2t(ii)    =ch2t(ii)-(1-hh1asset1t(ii))*passet1t(ii)+hh1bondt(ii)*pbondt(ii);
    ch1tilde(ii) = (alpha*ch1t(ii)^((pix1-1)/pix1)+(1-alpha)*hh1houset(ii)^((pix1-1)/pix1))^(pix1/(pix1-1));
    ch2tilde(ii) = (alpha*ch2t(ii)^((pix2-1)/pix2)+(1-alpha)*hh2houset(ii)^((pix2-1)/pix2))^(pix2/(pix2-1));
end

%--------------------------------------------------------
% First order conditions
%--------------------------------------------------------
Value1t=exp(Value1t);
Value2t=exp(Value2t);

q11=0;
q21=0;
qh1=0;
qh2=0;
qb1=0;
qb2=0;
switch exstate
    case 1
        transMatrix1=FH;
        transMatrix2=FH;
    case 2
        transMatrix1=FH;
        transMatrix2=FL;
    case 3
        transMatrix1=FL;
        transMatrix2=FH;
    case 4
        transMatrix1=FL;
        transMatrix2=FL;
    case 5
        transMatrix1=FH;
        transMatrix2=FH;
    case 6
        transMatrix1=FH;
        transMatrix2=FL;
    case 7
        transMatrix1=FL;
        transMatrix2=FH;
    case 8
        transMatrix1=FL;
        transMatrix2=FL;
end
EValue1=0;
EValue2=0;

for ii=1:8
    q11=q11+transMatrix1(exstate,ii)*Value1t(ii)^(1/psi1-gamma1)*growth(ii)^(1-gamma1)*(1-beta)*ch1tilde(ii)^(-1/psi1)*ch1tilde(ii)^(1/pix1)*alpha*ch1t(ii)^(1/(1-pix1))*(passet1t(ii)+div1(ii));
    q21=q21+transMatrix2(exstate,ii)*Value2t(ii)^(1/psi1-gamma1)*growth(ii)^(1-gamma2)*(1-beta)*ch2tilde(ii)^(-1/psi2)*ch2tilde(ii)^(1/pix2)*alpha*ch2t(ii)^(1/(1-pix2))*(passet1t(ii)+div1(ii));
    
    qb1 = qb1+transMatrix1(exstate,ii)*Value1t(ii)^(1/psi1-gamma1)*growth(ii)^(1-gamma1)*(1-beta)*ch1tilde(ii)^(-1/psi1)*ch1tilde(ii)^(1/pix1)*alpha*ch1t(ii)^(1/(1-pix1))/growth(ii);
    qb2 = qb2+transMatrix2(exstate,ii)*Value2t(ii)^(1/psi2-gamma1)*growth(ii)^(1-gamma2)*(1-beta)*ch2tilde(ii)^(-1/psi2)*ch2tilde(ii)^(1/pix2)*alpha*ch2t(ii)^(1/(1-pix2))/growth(ii);
    
    qh1=qh1+transMatrix1(exstate,ii)*Value1t(ii)^(1/psi1-gamma1)*growth(ii)^(1-gamma1)*(1-beta)*ch1tilde(ii)^(-1/psi2)*ch1tilde(ii)^(1/pix1)*alpha*ch1t(ii)^(1/(1-pix1))*(phouset(ii));
    qh2=qh2+transMatrix2(exstate,ii)*Value2t(ii)^(1/psi1-gamma1)*growth(ii)^(1-gamma2)*(1-beta)*ch2tilde(ii)^(-1/psi2)*ch2tilde(ii)^(1/pix2)*alpha*ch2t(ii)^(1/(1-pix2))*(phouset(ii));
    
    EValue1=EValue1+transMatrix1(exstate,ii)*(Value1t(ii)*growth(ii))^(1-gamma1);
    EValue2=EValue2+transMatrix2(exstate,ii)*(Value2t(ii)*growth(ii))^(1-gamma2);
end
ch1til=(alpha*ch1^((pix1-1)/pix1)+(1-alpha)*hh1houseguess^((pix1-1)/pix1))^(pix1/(pix1-1));
ch2til=(alpha*ch2^((pix2-1)/pix2)+(1-alpha)*hh2houseguess^((pix2-1)/pix2))^(pix2/(pix2-1));
phi1=(1-gamma1)/(1-1/psi1);
phi2=(1-gamma2)/(1-1/psi2);
Value1today=((1-beta)*ch1til^(1-1/psi1)+beta*EValue1^(1/phi1))^(phi1/(1-gamma1));
Value2today=((1-beta)*ch2til^(1-1/psi2)+beta*EValue2^(1/phi2))^(phi2/(1-gamma2));
MValue1today=(1-beta)*Value1today^(1/psi1)*ch1til^(-1/psi1)*alpha*ch1til^(1/pix1)*ch1^(1/(1-pix1));
MValue2today=(1-beta)*Value2today^(1/psi2)*ch2til^(-1/psi2)*alpha*ch2til^(1/pix2)*ch2^(1/(1-pix2));
MValueh1today=(1-beta)*Value1today^(1/psi1)*ch1til^(-1/psi1)*(1-alpha)*ch1til^(1/pix1)*hh1houseguess^(1/(1-pix1));
MValueh2today=(1-beta)*Value2today^(1/psi2)*ch2til^(-1/psi2)*(1-alpha)*ch2til^(1/pix1)*hh2houseguess^(1/(1-pix2));

ev1=Value1today;
ev2=Value2today;
%--------------------------------------------------------
% Calculate financial wealth share
%-------------------------------------------------------
hh1wsharet=zeros(16,1);
for ii=1:8
    hh1wsharet(ii)=(hh1asset1guess*(passet1t(ii)+div1(ii))+hh1houseguess*phouset(ii)+hh1bondguess/growth(ii))/(passet1t(ii)+div1(ii)+phouset(ii));
    fval(ii)=hh1wsharet(ii)-hh1wshareguess(ii);
    %fval(16+ii)=1-hh1wshareguess(ii)-hh2wshareguess(ii);
end
%---------------------------------------------------------------------
% First order conditions
%---------------------------------------------------------------------
fval(9)=q1g*MValue1today-(lag1assetshortpos+beta*Value1today^(1/psi1)*EValue1^((1-phi1)/phi1)*q11);
fval(10)=q1g*MValue2today-(lag2assetshortpos+beta*Value2today^(1/psi2)*EValue2^((1-phi2)/phi2)*q21);

fval(11)=qhg*MValue1today-(lag1houseshortpos+lag1collpos*margin1*pmin+MValueh1today+beta*Value1today^(1/psi1)*EValue1^((1-phi1)/phi1)*qh1);
fval(12)=qhg*MValue2today-(lag2houseshortpos+lag2collpos*margin1*pmin+MValueh2today+beta*Value2today^(1/psi2)*EValue2^((1-phi2)/phi2)*qh2);

fval(13)=qbg*MValue1today-(lag1collpos+beta*Value1today^(1/psi1)*EValue1^((1-phi1)/phi1)*qb1);
fval(14)=qbg*MValue2today-(lag2collpos+beta*Value2today^(1/psi2)*EValue2^((1-phi2)/phi2)*qb2);

%----------------------------------------------------------------------
% Budget constraint
%----------------------------------------------------------------------
income1=end1(exstate)+xgrid(endstate1)*(q1g+div1(exstate)+qhg)-hh1asset1guess*q1g-hh1bondguess*qbg-hh1houseguess*qhg;
fval(15)=ch1-income1;
fval(16)=1-ch1-ch2;
%-------------------------------------------------------------------
% Complementary slackness conditions
%-------------------------------------------------------------------
fval(17)=lag1collneg-(hh1houseguess*margin1*pmin+hh1asset1guess*margin2*pamin+hh1bondguess);
fval(18)=lag2collneg-(hh2houseguess*margin1*pmin+hh2asset1guess*margin2*pamin+hh2bondguess);
fval(19)=lag1assetshortneg-hh1asset1guess;
fval(20)=lag2assetshortneg-hh2asset1guess;
fval(21)=lag1houseshortneg-hh1houseguess;
fval(22)=lag2houseshortneg-hh2houseguess;
end

