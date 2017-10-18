clear

global exstate endstate1 xgrid
global FL FH
global gamma1 gamma2 beta
global div1 end1 end2
global passet11 passet12 passet13 passet14 passet15 passet16 passet17 passet18 
global hh1asset11 hh1asset12 hh1asset13 hh1asset14 hh1asset15 hh1asset16 hh1asset17 hh1asset18
global hh1bond1 hh1bond2 hh1bond3 hh1bond4 hh1bond5 hh1bond6 hh1bond7 hh1bond8
global pbond1 pbond2 pbond3 pbond4 pbond5 pbond6 pbond7 pbond8
global hh1house1 hh1house2 hh1house3 hh1house4 hh1house5 hh1house6 hh1house7 hh1house8
global phouse1 phouse2 phouse3 phouse4 phouse5 phouse6 phouse7 phouse8
global Value101 Value102 Value103 Value104 Value105 Value106 Value107 Value108
global Value201 Value202 Value203 Value204 Value205 Value206 Value207 Value208
global ch101 ch102 ch103 ch104 ch105 ch106 ch107 ch108
global ch201 ch202 ch203 ch204 ch205 ch206 ch207 ch208
global ev1 ev2 psi1 psi2 growth
global margin1 margin2 alpha
global pix1 pix2
global phousem
nendog=101;
%---------------------------------------------------------------------------
% Variables for the time iteration
%---------------------------------------------------------------------------
passet1n    =ones(nendog,8)*(-200);
pbondn      =zeros(nendog,8);
hh1asset1n  =zeros(nendog,8);
hh1bondn    =zeros(nendog,8);
hh1housen   =ones(nendog,8)*0.5;
phousen     =zeros(nendog,8);
ch1n        =zeros(nendog,8);
ch2n        =zeros(nendog,8);
Value1n     =ones(nendog,8);
Value2n     =ones(nendog,8);

passet1     =ones(nendog,8)*(-1);
pbond       =zeros(nendog,8);
hh1asset1   =zeros(nendog,8);
hh1bond     =zeros(nendog,8);
hh1house    =ones(nendog,8)*0.05;
phouse      =ones(nendog,8)*log(0.01);
ch1         =ones(nendog,8)*log(0.5);
ch2         =ones(nendog,8)*log(0.5);
Value1      =ones(nendog,8)*1;
Value2      =ones(nendog,8)*1;
optimSolutions=ones(nendog,8,19)*(0.01);
%------------------------------
% Parameterize the Economy
%-----------------------------

%-----------------------------
% Beliefs
%-----------------------------

a1=0.50;
a2=0.14;
a3=0.14;
a4=0.14;

phix=0.43;
alpha1=0.57;
alpha2=0.57;
etah=1.60;

A=[a1, alpha1-a1, alpha2-a1, 1+a1-alpha1-alpha2;
   a2, alpha1-a2, alpha2-a2, 1+a2-alpha1-alpha2;
   a3, alpha1-a3, alpha2-a3, 1+a3-alpha1-alpha2;
   a4, alpha1-a4, alpha2-a4, 1+a4-alpha1-alpha2];

transMat=[0.43 0.57;
          0.57 0.43];
      
Gamma=[transMat(1,1)*A transMat(1,2)*A;
       transMat(2,1)*A transMat(2,2)*A];
   

FH=[   etah*transMat(1,1)*A transMat(1,2)*(1-etah*(transMat(1,1)))/(transMat(1,2))*A;
       etah*transMat(2,1)*A transMat(2,2)*(1-etah*(transMat(2,1)))/(transMat(2,2))*A];

FL=1/(1-alpha1)*(Gamma-alpha1*FH);
%-----------------------------
% Preferences
%-----------------------------
alpha=0.90;
gamma1=2.00;
gamma2=2.00;
psi1=1.5;
psi2=1.5;
beta=0.96;
pix1=1.15; %1.25 standard case
pix2=1.15; %1.25 standard case

Value1=log((1-beta)*0.5)*Value1;
Value2=log((1-beta)*0.5)*Value2;

%---------------------------------
% Economic Fundamentals
%---------------------------------
aggGrowthHigh=1.054;
aggGrowthLow=0.982;
growth=[aggGrowthHigh aggGrowthHigh aggGrowthHigh aggGrowthHigh aggGrowthLow aggGrowthLow aggGrowthLow aggGrowthLow];
div1=[0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15];
eh_growth=0.5;
el_growth=1-eh_growth;
eh_recess=0.5;
el_recess=1-eh_recess;
margin1=0.50; % Collateral of housing
margin2=0.00; % Collateral of stock
end1=[eh_growth eh_growth eh_growth eh_growth eh_growth eh_growth eh_growth eh_growth];
end2=1-end1;
end1=end1.*(1-div1);
end2=end2.*(1-div1);
xgrid=linspace(0.00,1.000,nendog);

%--------------------------------------
% Set parameters for solvers
%--------------------------------------
opts=optiset('solver','ipopt','maxiter',10000,'display','iter','tolafun',1e-7,'tolrfun',1e-7);
options = optimoptions(@fmincon,'MaxIter',5000000,'MaxFunEvals',500000000,'Display','iter','Algorithm','interior-point','SubproblemAlgorithm','cg','HessianApproximation','finite-difference','FinDiffType','central','TolCon',1e-7);
%-------------------------------------
%Load old solutions
%-------------------------------------
load('Solutions_tempa.mat');
%load('Solutions_2_g200_m050_p150_e160.mat');
max_error=100;
time_iter=0;
%-------------------------------
% Start time iteration
%-------------------------------
while max_error>1e-4
    time_iter=time_iter+1;
    disp(time_iter);
    passet11=pchip(xgrid,squeeze(passet1(:,1)));
	passet12=pchip(xgrid,squeeze(passet1(:,2)));
	passet13=pchip(xgrid,squeeze(passet1(:,3)));
	passet14=pchip(xgrid,squeeze(passet1(:,4)));
	passet15=pchip(xgrid,squeeze(passet1(:,5)));
	passet16=pchip(xgrid,squeeze(passet1(:,6)));
	passet17=pchip(xgrid,squeeze(passet1(:,7)));
	passet18=pchip(xgrid,squeeze(passet1(:,8)));
   
    Value101=pchip(xgrid,squeeze(Value1(:,1)));
    Value102=pchip(xgrid,squeeze(Value1(:,2)));
    Value103=pchip(xgrid,squeeze(Value1(:,3)));
    Value104=pchip(xgrid,squeeze(Value1(:,4)));
    Value105=pchip(xgrid,squeeze(Value1(:,5)));
    Value106=pchip(xgrid,squeeze(Value1(:,6)));
    Value107=pchip(xgrid,squeeze(Value1(:,7)));
    Value108=pchip(xgrid,squeeze(Value1(:,8)));

    Value201=pchip(xgrid,squeeze(Value2(:,1)));
    Value202=pchip(xgrid,squeeze(Value2(:,2)));
    Value203=pchip(xgrid,squeeze(Value2(:,3)));
    Value204=pchip(xgrid,squeeze(Value2(:,4)));
    Value205=pchip(xgrid,squeeze(Value2(:,5)));
    Value206=pchip(xgrid,squeeze(Value2(:,6)));
    Value207=pchip(xgrid,squeeze(Value2(:,7)));
    Value208=pchip(xgrid,squeeze(Value2(:,8)));

    ch101=pchip(xgrid,squeeze(ch1(:,1)));
    ch102=pchip(xgrid,squeeze(ch1(:,2)));
    ch103=pchip(xgrid,squeeze(ch1(:,3)));
    ch104=pchip(xgrid,squeeze(ch1(:,4)));
    ch105=pchip(xgrid,squeeze(ch1(:,5)));
    ch106=pchip(xgrid,squeeze(ch1(:,6)));
    ch107=pchip(xgrid,squeeze(ch1(:,7)));
    ch108=pchip(xgrid,squeeze(ch1(:,8)));

    ch201=pchip(xgrid,squeeze(ch2(:,1)));
    ch202=pchip(xgrid,squeeze(ch2(:,2)));
    ch203=pchip(xgrid,squeeze(ch2(:,3)));
    ch204=pchip(xgrid,squeeze(ch2(:,4)));
    ch205=pchip(xgrid,squeeze(ch2(:,5)));
    ch206=pchip(xgrid,squeeze(ch2(:,6)));
    ch207=pchip(xgrid,squeeze(ch2(:,7)));
    ch208=pchip(xgrid,squeeze(ch2(:,8)));
    
    hh1asset11=pchip(xgrid,squeeze(hh1asset1(:,1)));
    hh1asset12=pchip(xgrid,squeeze(hh1asset1(:,2)));
    hh1asset13=pchip(xgrid,squeeze(hh1asset1(:,3)));
    hh1asset14=pchip(xgrid,squeeze(hh1asset1(:,4)));
    hh1asset15=pchip(xgrid,squeeze(hh1asset1(:,5)));
    hh1asset16=pchip(xgrid,squeeze(hh1asset1(:,6)));
    hh1asset17=pchip(xgrid,squeeze(hh1asset1(:,7)));
    hh1asset18=pchip(xgrid,squeeze(hh1asset1(:,8)));
   
    hh1bond1 = pchip(xgrid,squeeze(hh1bond(:,1)));
    hh1bond2 = pchip(xgrid,squeeze(hh1bond(:,2)));
    hh1bond3 = pchip(xgrid,squeeze(hh1bond(:,3)));
    hh1bond4 = pchip(xgrid,squeeze(hh1bond(:,4)));
    hh1bond5 = pchip(xgrid,squeeze(hh1bond(:,5)));
    hh1bond6 = pchip(xgrid,squeeze(hh1bond(:,6)));
    hh1bond7 = pchip(xgrid,squeeze(hh1bond(:,7)));
    hh1bond8 = pchip(xgrid,squeeze(hh1bond(:,8)));
    
    pbond1 = pchip(xgrid,squeeze(pbond(:,1)));
    pbond2 = pchip(xgrid,squeeze(pbond(:,2)));
    pbond3 = pchip(xgrid,squeeze(pbond(:,3)));
    pbond4 = pchip(xgrid,squeeze(pbond(:,4)));
    pbond5 = pchip(xgrid,squeeze(pbond(:,5)));
    pbond6 = pchip(xgrid,squeeze(pbond(:,6)));
    pbond7 = pchip(xgrid,squeeze(pbond(:,7)));
    pbond8 = pchip(xgrid,squeeze(pbond(:,8)));
    
    hh1house1=pchip(xgrid,squeeze(hh1house(:,1)));
    hh1house2=pchip(xgrid,squeeze(hh1house(:,2)));
    hh1house3=pchip(xgrid,squeeze(hh1house(:,3)));
    hh1house4=pchip(xgrid,squeeze(hh1house(:,4)));
    hh1house5=pchip(xgrid,squeeze(hh1house(:,5)));
    hh1house6=pchip(xgrid,squeeze(hh1house(:,6)));
    hh1house7=pchip(xgrid,squeeze(hh1house(:,7)));
    hh1house8=pchip(xgrid,squeeze(hh1house(:,8)));
    
    phouse1=pchip(xgrid,squeeze(phouse(:,1)));
    phouse2=pchip(xgrid,squeeze(phouse(:,2)));
    phouse3=pchip(xgrid,squeeze(phouse(:,3)));
    phouse4=pchip(xgrid,squeeze(phouse(:,4)));
    phouse5=pchip(xgrid,squeeze(phouse(:,5)));
    phouse6=pchip(xgrid,squeeze(phouse(:,6)));
    phouse7=pchip(xgrid,squeeze(phouse(:,7)));
    phouse8=pchip(xgrid,squeeze(phouse(:,8)));
   
	for exstate=1:8
        disp(exstate);
		xguess=squeeze(optimSolutions(1,exstate,:));
        if exstate>1
            xguess=squeeze(optimSolutions(1,exstate-1,:));
        end
        
        %xguess=ones(40,1)*(0.01);
		for endstate1=1:nendog
            switch exstate
                case 1
                    phousem=exp(ppval(phouse1,endstate1));
                case 2
                    phousem=exp(ppval(phouse2,endstate1));
                case 3
                    phousem=exp(ppval(phouse3,endstate1));
                case 4
                    phousem=exp(ppval(phouse4,endstate1));
                case 5
                    phousem=exp(ppval(phouse5,endstate1));
                case 6
                    phousem=exp(ppval(phouse6,endstate1));
                case 7
                    phousem=exp(ppval(phouse7,endstate1));
                case 8
                    phousem=exp(ppval(phouse8,endstate1));
            end
            %disp(endstate1);
            %f=endstate1;
            %endstate=20;
            %if time_iter>50
               xguess=squeeze(optimSolutions(endstate1,exstate,:));
               %xguess=ones(42,1)*(-0.01);
            %end
            if time_iter<3
                disp(endstate1);
            end
            %[x,z,exitflag]=opti_fmincon(@EulerEquation,xguess,[],[],[],[],[],[],[],opts);
            %xguess=x;
            %pause;
            [x,z,exitflag]=fmincon(@FDummy,xguess,[],[],[],[],[],[],@EulerEquation,options);
            %[x,z,exitflag]=fsolve(@EulerEquation,xguess,options);
            xguess=x;
            %endstate1=z;
			passet1n(endstate1,exstate)=x(17);
            pbondn(endstate1,exstate)=x(18);
            phousen(endstate1,exstate)=x(19);
			hh1asset1n(endstate1,exstate)=1;
            hh1bondn(endstate1,exstate)=x(9);
            hh1housen(endstate1,exstate)=x(10);
            ch1n(endstate1,exstate)=x(11);
            ch2n(endstate1,exstate)=x(12);
            Value1n(endstate1,exstate)=log(ev1);
            Value2n(endstate1,exstate)=log(ev2);
            optimSolutions(endstate1,exstate,:)=x;
            %[a b]=EulerEquation(x);
            %disp(b);
            
            %pause;
        end
        %pause;
    end
    diff_asset=max(max(abs(exp(passet1n)-exp(passet1))));
	diff_bond = max(max(abs(exp(pbondn)-exp(pbond))));
    diff_house = max(max(abs(exp(phousen)-exp(phouse))));
    diff_hh1a1= max(max(abs((hh1asset1n)-(hh1asset1))));
	diff_hh1b1= max(max(abs((hh1bondn)-(hh1bond))));
    diff_hh1h1 =max(max(abs(hh1housen-hh1house)));
    diff=[diff_asset diff_bond diff_house diff_hh1a1 diff_hh1b1 diff_hh1h1];
	max_error=max(diff);
	disp('Maximum Differences in Prices and Portfolios');
	disp(diff);
    %pause;
	passet1=passet1n;
	pbond=pbondn;
    phouse=phousen;
	ch1=ch1n;
	ch2=ch2n;
    Value1=Value1n;
    Value2=Value2n;
	hh1asset1=hh1asset1n;
	hh1bond=hh1bondn;
    hh1house=hh1housen;
    

    save('Solutions_tempa.mat','passet1','pbond','phouse','hh1asset1','hh1bond','hh1house','ch1','ch2','Value1','Value2','optimSolutions','max_error','diff');
    %pause;
    %options = optimset('Display','iter','TolCon',1e-10,'TolX',1e-30,'MaxIter',50000);
end
save('Solutions_2_g200_m050_p150_pix105_e160.mat','passet1','pbond','phouse','hh1asset1','hh1bond','hh1house','ch1','ch2','Value1','Value2','optimSolutions','max_error');

		