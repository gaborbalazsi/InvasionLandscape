clearvars; close all;

%*****************************************************************************80
%
%% OU_EXACT applies the exact method to the Ornstein-Uhlenbeck SDE.
%
%  Original code by John Burkardt was modified by Gábor Balázsi.   
%
%  Reference:
%
%    Daniel T. Gillespie,
%    Exact numerical simulation of the Ornstein-Uhlenbeck process 
%    and its integral,
%    Physical review. E, 54(2):2084-2091 · September 1996.
%
%  Parameters:
%
%    Input, real TMAX, NTSTEPS, NCELLS, THRESHOLD, SEED, the value of problem parameters.
%
%
%    Input, real TMAX, the final time.
%
%    Input, integer NTSTEPS, the number of time steps.
%
%    Input, integer NCELLS, the maximum number of cells.
%
%    Input, real THRESHOLD, the threshold for survival.
%
%    Input, integer SEED, a seed for the random number generator.
%
% for example: ou_second_NFPFstats(150, 10000, 1000, 0.24, 'shuffle');

tmax=50;            %simulation (end) time (hours)
ntsteps=tmax*10;    %number of time steps
ncells=5000;        %number of cells
% thetaNF = 10/(tmax/ntsteps);
thetaNF = 1000.0;   %relaxation (autocorrelation) time in units of tsteps
ProbAmp=0.25;       %invasion probability factor
seed='shuffle';     %for random number generator

% quadratic fit estimates experimental CV from experiemntal mean
p1 =    0.008908;
p2 =    -0.05025;
p3 =      0.1378;
quadfit=@(x) p1*x.^2 + p2*x + p3;

doxRange=[    0    0.1    0.3    0.35    0.5    0.6    1    2    10];

%binned ladscape estimates
BLdata=load('LandscIndBL.txt'); %Landscape data for each individual Dox
%BLcons=load('aveLandscBL.txt');
BLcons=load('smLandscBL.txt');
xCentExp=load('xCentersBL.txt');

%plot individual landscapes
figure;
plot(xCentExp,BLdata);
legCell=cellstr(num2str(doxRange', 'Dox=%.2f'));legend(legCell,'Location','EastOutside')
title('Individual landscapes');
xlabel('lg(BACH1)');ylabel('Invasiveness');
set(gca,'YLim',[0 1])

%inducer concentrations (padded at both ends)
dox=[...
    -1
    0
    0.1
    0.3
    0.35
    0.5
    0.6
    1
    2
    10
    20];

%average fluorescence of Seeded BL cells (padded at both ends)
GFPave=[...
    2.0
    2.979016133
    3.398626132
    3.648960343
    3.788605752
    3.809999739
    3.823816377
    4.062188809
    4.202669121
    4.322541832
    10];

%CV of Seeded BL cells (padded at both ends)
GFPcv=[...
    0.05    
    0.0670
    0.0690
    0.0757
    0.0637
    0.0751
    0.0764
    0.0797
    0.0802
    0.0902
    0.1];

GFPstd=GFPcv.*GFPave;

%invasiveness (Low mNF-BACH1)  (padded at both ends with bInv, sInv)
bInv=0;
sInv=100;
InvFL=[...
    bInv    bInv    bInv
    27.55 	25.66 	24.27 
    30.92 	30.29 	33.33 
    35.93 	31.15 	35.86 
    30.39 	26.63 	25.13 
    13.97 	11.42 	15.74 
    16.25 	17.58 	18.06 
    16.62 	21.42 	22.65 
    23.79 	29.70 	31.43 
    38.35 	36.83 	39.37 
    sInv    sInv    sInv];

avInv=mean(InvFL(2:end-1,:)');   % mean of 3 invasion replicates, unpadded
stInv=std(InvFL(2:end-1,:)');  % std of 3 invasion replicates, unpadded
minInv=min(avInv);  % overall landscape miminmum, unpadded
nInv=(avInv-minInv);  %normalized landscape with minimum reset to = 0
mInv=avInv;

%average GFP fluorescence of Invaded BL cells (not padded)
invGFPave=[...
    3.37356879
    3.533857401
    3.609787441
    3.735868215
    3.785768841
    3.859890899
    4.394765724
    4.352679483
    4.502049469];

%CV of Invaded BL cells (not padded)
invGFPcv=[...
    0.054
    0.0682
    0.0744
    0.0427
    0.1238
    0.1140
    0.0377
    0.0701
    0.0611];

GFPstd=GFPcv.*GFPave;

%these will be the simulated values for invading cells
percInv=zeros(9,10);
aveInv=zeros(9,10);
stdInv=zeros(9,10);
cvrInv=zeros(9,10);

xx=linspace(0,11,1000);

ctr=1;
for doxC=doxRange;
%current Dox concentration
    indD=find(dox==doxC);   %where on the landscape?
    currind=indD-1;
    gfpC=GFPave(indD);      %what is the current mean GFP expression level?
    sigmaNF=GFPstd(indD);   %what is the current std of GFP expression level?
    muNF = gfpC;
    varNF=sigmaNF*sigmaNF;  %GFP expression variance
    
    ExpLandsc=BLdata(currind,:);    %individual landscape
    %ExpLandsc=BLcons;              %consolidated landscape
    
    %read in actual data and create histograms
    xEdges=linspace(2.5,5,50);yEdges=xEdges;
    sheetname=sprintf('%d',currind);
    ExpData=xlsread('./MB231_1.1BLdataGB.xlsx',sheetname);
    ExpSeeded=ExpData(:,1:3);ExpSeeded=ExpSeeded(~isnan(ExpSeeded(:)));
    [sBLexp(currind,:),~]=histcounts(ExpSeeded,xEdges,'Normalization','pdf');
    ExpInvded=ExpData(:,4:6);ExpInvded=ExpInvded(~isnan(ExpInvded(:)));
    [iBLexp(currind,:),~]=histcounts(ExpInvded,xEdges,'Normalization','pdf');
    iBLexp(currind,:)=mInv(currind).*iBLexp(currind,:);

    %a quick estimate for cov_z_z2
    ExpSeeded=muNF+sigmaNF*randn(21000,1);
    cov0=cov(ExpSeeded,ExpSeeded.^2);
    cov_z_z2=cov0(1,2);
    [cov_z_z2 2*muNF*varNF]
    
    cov_z_z2=2*muNF*varNF;
   
    xRange=xx(xx >= muNF-sigmaNF/0.5 & xx < muNF+sigmaNF/0.5);
    LandSc=pchip(xCentExp,ExpLandsc,xRange);
    [c1,q1]=polyfit(xRange,LandSc,1);   %linear fit to current landsc. region
    [c2,q2]=polyfit(xRange,LandSc,2);   %quadratic fit to current region

    a1=c1(1);b1=c1(2);  %slope and intercept from linear fit
    a2=c2(1);b2=c2(2);o2=c2(3);  %quadratic fit coefficients
    
    cov11=cov(a1*ExpSeeded+b1,ExpSeeded);
    cov1wz1=cov11(1,2);
    cov12=cov(a1*ExpSeeded+b1,ExpSeeded.^2);
    cov1wz2=cov12(1,2);
    
    cov21=cov(a2*ExpSeeded.^2+b2*ExpSeeded+o2,ExpSeeded);
    cov2wz1=cov21(1,2);
    cov22=cov(a2*ExpSeeded.^2+b2*ExpSeeded+o2,ExpSeeded.^2);
    cov2wz2=cov22(1,2);
    
    Dmean1=a1*varNF/mean(a1*ExpSeeded+b1);    %1st expected change in mean (from theory: linear landscape)
    Dmean2=(2*a2*muNF+b2)*varNF/mean(a2*ExpSeeded.^2+b2*ExpSeeded+o2);   %2nd expect. change in mean (quadratic landscape)
%    Dmean2=(a2*cov_z_z2+b2*varNF)/mean(a2*ExpSeeded.^2+b2*ExpSeeded+o2);   %2nd expect. change in mean (quadratic landscape)

%either of the next 2 works:
    Dvar1=a1*(cov_z_z2-2*muNF*varNF)/mean(a1*ExpSeeded+b1)-(a1*varNF)^2/mean(a1*ExpSeeded+b1)^2;    %1st expected change in variance (from theory: linear landscape)
%    Dvar1=cov1wz2/(a1*mean(ExpSeeded)+b1) - cov1wz1*(mean((a1*ExpSeeded+b1).*ExpSeeded) + mean(a1*ExpSeeded+b1)*mean(ExpSeeded))/mean(a1*ExpSeeded+b1)^2;
    
%next one works (2nd expected change in variance: quadratic landscape)):
    Dvar2=cov2wz2/mean(a2*ExpSeeded.^2+b2*ExpSeeded+o2) - cov2wz1*(mean((a2*ExpSeeded.^2+b2*ExpSeeded+o2).*ExpSeeded) + mean(a2*ExpSeeded.^2+b2*ExpSeeded+o2)*mean(ExpSeeded))/mean(a2*ExpSeeded.^2+b2*ExpSeeded+o2)^2;
    
    Nmean1(ctr)=muNF+Dmean1;  %new mean (from linear landscape theory) for invaded cells
    Nvar1(ctr)=varNF+Dvar1;   %new variance (from linear landscape theory)
    
    Nmean2(ctr)=muNF+Dmean2;  %new mean (from quad theory) for invaded cells
    Nvar2(ctr)=varNF+Dvar2;   %new variance (from quad theory)
    
   figure;hold on;
   plot(xRange,LandSc,'LineWidth',2);
   plot(xRange,a1*xRange+b1,'b--');
   plot(xRange,a2*xRange.^2+b2*xRange+o2,'c--');
   plot(gfpC,mInv(currind)/100,'o')
   legend('orig','lin','quad')
   title(sprintf('Dox=%.2f',doxC));
   

    NMAX=ncells*2;  %maximum number of cells
    growthRate=0.0; %growth turned off if = 0

    % simulate NF
    for rep=1:10
    %
    %  Initialize the random number generator.
    %  The following way to initialize the random number generator 
    %  may not be available in older versions of MATLAB.
    %
        rng ( seed )
    %
    %  Set the discrete time stepsize.
    %
        dt = tmax / ntsteps;
    %
    %  Compute the Brownian increments.
    %
        dw = randn ( NMAX, ntsteps );
    %
    %  Carry out the exact simulation for NF.
    %
        t = linspace ( 0, tmax, ntsteps + 1 );
        xNF = NaN*zeros ( NMAX, ntsteps + 1 );
        NcurrentNF = zeros ( 1 , ntsteps + 1 );
        NinvNF = zeros ( 1 , ntsteps + 1 );

        xNF(1:ncells,1) = muNF + sigmaNF * randn(ncells,1); %prepopulate gene expression 
        yNF=[];
        NcurrentNF(1)=ncells;
        for j = 1 : ntsteps
            xNF(:,j+1) = muNF + (xNF(:,j)-muNF)*exp(-dt/thetaNF) + sqrt( (1-exp(-2*dt/thetaNF))*sigmaNF^2 ) * dw(:,j);
            yNF = muNF + (yNF-muNF)*exp(-dt/thetaNF) + sqrt( (1-exp(-2*dt/thetaNF))*sigmaNF^2 ) * randn ( size(yNF) );
            prbInv=ProbAmp * 0.01*pchip(xCentExp,ExpLandsc,xNF(:,j+1));
            %prbInv=-ProbAmp *log(1-pchip(xCentExp,ExpLandsc,xNF(:,j+1)));
            rndInv=rand(size(prbInv));
            indInv=find(prbInv>rndInv);
            yNF=[yNF; xNF(indInv,j+1)]; % add next set of invading cells to yNF
            xNF(indInv,j+1)=NaN;        % remove next set of invading cells from xNF
            NcurrentNF(j+1)=sum(isfinite(xNF(:,j+1)));
            NinvNF(j+1)=sum(isfinite(yNF));
        end
%         plot ( t, NcurrentNF, 'b-', 'LineWidth', 2 );
%         plot ( t, NinvNF, 'c-', 'LineWidth', 2 );
        percInv(ctr,rep)=NinvNF(end)/ncells;
        aveInv(ctr,rep)=mean(yNF);
        stdInv(ctr,rep)=std(yNF);
        cvrInv(ctr,rep)=std(yNF)/mean(yNF);
        xNFp=reshape(xNF,numel(xNF),1);xNFp=xNFp(~isnan(xNFp));
        [xNFh(rep,:),xEdges]=histcounts(xNFp,linspace(2.5,5,100),'Normalization','Probability');
        [yNFh(rep,:),yEdges]=histcounts(yNF,linspace(2.5,5,100),'Normalization','Probability');
    end
    
    %plot control and invasion histograms
    figure; hold on;
    xCenters=(xEdges(1:end-1)+xEdges(2:end))/2;yCenters=(yEdges(1:end-1)+yEdges(2:end))/2;
    plot(xCenters,mean(xNFh),'LineWidth',2,'Color',[0.5 0.5 0.5]);
    plot(xCenters,mean(yNFh),'LineWidth',2,'Color',[0.53 0 1]);
    
    legend('simSEED','simINV','gftSEED','gxpSEED','gxpINV','expSEED','expINV','Location','NorthWest')
    xlabel ( 'BACH1 level (arb. units)', 'FontSize', 24)
    ylabel ( 'Probability', 'FontSize', 24);
    tstr=sprintf("Dox=%.2f",doxC);
    title ( tstr, 'FontSize', 24 )
    grid on; hold on;
    set(gca,'XLim',[2.5 5],'YLim',[0 0.1], 'FontSize', 24)
    
    %experimental pre-invasion values
    aveCtr(ctr)=muNF;
    stdCtr(ctr)=sigmaNF;
    cvrCtr(ctr)=sigmaNF/muNF;
    doxC
    ctr=ctr+1;
end

%prepare a normal distribution to see where it is on landscape
xx=linspace(0,7,1000);

BHchip=load('BHchip.mat');
BHexp=load('BHexp.mat');

%plot invasion data and interpolated landscape
figure;hold on;
errorbar(GFPave(2:end-1),mInv,stInv,'bo','LineWidth',2);
plot(GFPave(2:end-1),100*mean(percInv'),'c*','LineWidth',2);
plot(xCentExp,100*BLcons,'b--','LineWidth',2);
xlabel('log_{10}(GFP MFI) (arb. units)');ylabel('Invasiveness(%)');
legend('expt.','sim.','inf.','Location','SouthEast');
box on;yMax=50;
set(gca,'FontSize',24,'XLim',[2.5 5.0],'YLim',[0 yMax]);
text(2.5,yMax+2,'Dox=','FontSize',16,'Color','b');
%text(GFPave(2)-0.1,52,'0.0','FontSize',16,'Color','b');
text(GFPave(2)-0.08,yMax+2,'0.0','FontSize',16,'Color','b');
% text(GFPave(4)-0.05,105,'0.3','FontSize',20);
text(GFPave(6)-0.08,yMax+2,'0.5','FontSize',16,'Color','b');
% text(GFPave(7)-0.05,105,'1.0','FontSize',20);
text(GFPave(end-1)-0.06,yMax+2,'10','FontSize',16,'Color','b');

%inducer post-invasion
doxPL=[...
    0
    0.1
    0.3
    0.35
    0.5
    0.6
    1
    2
    10];

%average fluorescence post-invasion
GFPavePL=[...
    3.37356879
    3.533857401
    3.609787441
    3.735868215
    3.785768841
    3.859890899
    4.394765724
    4.352679483
    4.502049469];

%plot observed & predicted values for the MEAN before & after invasion,
%linear landscape theory
figure;hold on;
plot([doxRange(1:end-1) 10],aveCtr,'ko','LineWidth',2);
errorbar([doxRange(1:end-1) 10],mean(aveInv'),std(aveInv'),'ms','LineWidth',2);
plot([doxRange(1:end-1) 10],Nmean1,'rs','MarkerSize',12); %theory
plot([doxPL(1:end-1)' 10],GFPavePL,'o','MarkerSize',12,'Color',[0.53 0 1]);
xlabel ( 'Dox (ng/ml)', 'FontSize', 24)
ylabel ( 'mean of lg(GFP)', 'FontSize', 24);
grid on;
%set(gca,'FontSize', 24,'XLim',[0 2.5],'YLim',[3.3 4.65],'XTick',[0 1 2],'XTickLabel',[0 1 10]);
set(gca,'FontSize', 24);
L1=legend('ctr.','sim.','th1.','expt.','Location','EastOutside');
% line([1.2 1.4],[3.3 3.5],'LineWidth',2,'Color','k')
% line([1.3 1.5],[3.3 3.5],'LineWidth',2,'Color','k')
% line([1.2 1.4],[4.45 4.65],'LineWidth',2,'Color','k')
% line([1.3 1.5],[4.45 4.65],'LineWidth',2,'Color','k')
set(L1,'String',{'ctr.' 'sim.' 'th1.' 'expt.' })
title('mNF-BL, linear mean');

%plot observed & predicted values for the MEAN before & after invasion,
%quad landscape theory
figure;hold on;
plot([doxRange(1:end-1) 10],aveCtr,'ko','LineWidth',2);
errorbar([doxRange(1:end-1) 10],mean(aveInv'),std(aveInv'),'ms','LineWidth',2);
plot([doxRange(1:end-1) 10],Nmean2,'rs','MarkerSize',12); %theory
plot([doxPL(1:end-1)' 10],GFPavePL,'o','MarkerSize',12,'Color',[0.53 0 1]);
xlabel ( 'Dox (ng/ml)', 'FontSize', 24)
ylabel ( 'mean of lg(GFP)', 'FontSize', 24);
grid on;
set(gca,'FontSize', 24);
L1=legend('ctr.','sim.','th2.','expt.','Location','EastOutside');
set(L1,'String',{'ctr.' 'sim.' 'th2.' 'expt.' })
title('mNF-BL, quad mean');

%fluorescence CV post-invasion
GFPcvPL=[...
    0.0540
    0.0682
    0.0744
    0.0427
    0.1238
    0.1140
    0.0377
    0.0701
    0.0611];

GFPcvPL=GFPcvPL;

%plot observed & predicted values for the CV before & after invasion,
%linear landscape theory
figure;hold on;
plot([doxRange(1:end-1) 10],cvrCtr,'ko','LineWidth',2);
errorbar([doxRange(1:end-1) 10],mean(cvrInv'),std(cvrInv'),'ms','LineWidth',2);
plot([doxRange(1:end-1) 10],sqrt(Nvar1)./Nmean1,'rs','MarkerSize',12);
plot([doxPL(1:end-1)' 10],GFPcvPL,'o','MarkerSize',12,'Color',[0.53 0 1]);
xlabel ( 'Dox (ng/ml)', 'FontSize', 24)
ylabel ( 'CV of lg(GFP)', 'FontSize', 24);
grid on;
%set(gca,'FontSize', 24,'XLim',[0 2.5],'XTick',[0 1 2],'XTickLabel',[0 1 10])
set(gca,'FontSize', 24)
L2=legend('ctr.','sim.','th1.','expt.','Location','EastOutside');
% line([1.2 1.4],[0.0 0.02],'LineWidth',2,'Color','k')
% line([1.3 1.5],[0.0 0.02],'LineWidth',2,'Color','k')
% line([1.2 1.4],[0.08 0.1],'LineWidth',2,'Color','k')
% line([1.3 1.5],[0.08 0.1],'LineWidth',2,'Color','k')
set(L2,'String',{'ctr.' 'sim.' 'th1.' 'expt.' })
title('mNF-BL, linear CV');

%plot observed & predicted values for the CV before & after invasion,
%linear landscape theory
figure;hold on;
plot([doxRange(1:end-1) 10],cvrCtr,'ko','LineWidth',2);
errorbar([doxRange(1:end-1) 10],mean(cvrInv'),std(cvrInv'),'ms','LineWidth',2);
plot([doxRange(1:end-1) 10],sqrt(Nvar2)./Nmean2,'rs','MarkerSize',12);
plot([doxPL(1:end-1)' 10],GFPcvPL,'o','MarkerSize',12,'Color',[0.53 0 1]);
xlabel ( 'Dox (ng/ml)', 'FontSize', 24)
ylabel ( 'CV of lg(GFP)', 'FontSize', 24);
grid on;
%set(gca,'FontSize', 24,'XLim',[0 2.5],'XTick',[0 1 2],'XTickLabel',[0 1 10])
set(gca,'FontSize', 24)
L2=legend('ctr.','sim.','th2.','expt.','Location','EastOutside');
% line([1.2 1.4],[0.0 0.02],'LineWidth',2,'Color','k')
% line([1.3 1.5],[0.0 0.02],'LineWidth',2,'Color','k')
% line([1.2 1.4],[0.08 0.1],'LineWidth',2,'Color','k')
% line([1.3 1.5],[0.08 0.1],'LineWidth',2,'Color','k')
set(L2,'String',{'ctr.' 'sim.' 'th2.' 'expt.' })
title('mNF-BL, quad CV');

% figure;hold on;
% plot ( t, NcurrentNF, 'b-', 'LineWidth', 2 );
% plot ( t, NinvNF, 'c-', 'LineWidth', 2 );
% xlabel ( 't (hours)', 'FontSize', 20 )
% ylabel ( 'N(t)', 'FontSize', 20, 'HorizontalAlignment', 'right' )
% legend('NF','NFinv');
% set(gca,'FontSize',20);
% 
% figure;
% plot ( t, xNF(1:100,:), 'k-', 'LineWidth', 1 )
% xlabel ( 't (hours)', 'FontSize', 16 )
% ylabel ( 'X(t)', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right' )
% title ( 'Exact simulation of O-U SDE', 'FontSize', 16 )
% grid on;set(gca,'FontSize',20);

figure;hold on;
plot(Nmean1-aveCtr,mean(aveInv')-aveCtr,'ms','LineWidth',2);   %linear theory, simulation
plot(Nmean1-aveCtr,GFPavePL'-aveCtr,'ro','LineWidth',2); %linear theory, expt.
plot(Nmean2-aveCtr,mean(aveInv')-aveCtr,'cs','LineWidth',2);   %quad theory, simulation
plot(Nmean2-aveCtr,GFPavePL'-aveCtr,'bo','LineWidth',2); %quad theory, expt.
legend('lin-sim','lin-expt','quad-sim','quad-expt','Location','NorthWest');
xlabel('Theory mean');ylabel('Simulation & expt. means');
set(gca,'FontSize',20);

figure;hold on;
plot(sqrt(Nvar1)./Nmean1-cvrCtr,mean(cvrInv')-cvrCtr,'ms','LineWidth',2);
plot(sqrt(Nvar1)./Nmean1-cvrCtr,GFPcvPL'-cvrCtr,'ro','LineWidth',2);
plot(sqrt(Nvar2)./Nmean2-cvrCtr,mean(cvrInv')-cvrCtr,'cs','LineWidth',2);
plot(sqrt(Nvar2)./Nmean2-cvrCtr,GFPcvPL'-cvrCtr,'bo','LineWidth',2);
legend('lin-sim','lin-expt','quad-sim','quad-expt','Location','NorthWest');
xlabel('Theory CV');ylabel('Simulation & expt. CVs');
set(gca,'FontSize',20);
