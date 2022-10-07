% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% 2022/10/06
% Yusuke Mukuhira
% Institue of fluid science, Tohoku University
% tohoku@tohoku.ac.jp
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

close all;
clear;

fnc_pdir = './function/';
path(path,fnc_pdir);

wgs84 = wgs84Ellipsoid('kilometer');

%----------[Color map]------------------------------------------------------------
cm(1,:)=  [     0    0.4470    0.7410];
cm(2,:)=  [0.8500    0.3250    0.0980];
cm(3,:)=  [0.9290    0.6940    0.1250];
cm(4,:)=  [0.4940    0.1840    0.5560];
cm(5,:)=  [0.4660    0.6740    0.1880];
cm(6,:)=  [0.3010    0.7450    0.9330];
cm(7,:)=  [0.6350    0.0780    0.1840];

bluesmap = brewermap(16, 'Blues');
redsmap  = brewermap(16, 'Reds');
grnsmap  = brewermap(16, 'Greens');
prplsmap = brewermap(16, 'Purples');


% ---------[Read entire catalog]------------------
filename = '../Data/post.out';
delimiterIn = ' ';
headerlinesIn = 0;
Data = importdata(filename,delimiterIn,headerlinesIn);
year=Data(:,1);month=Data(:,2);day=Data(:,3);
hour=Data(:,4);minut=Data(:,5);sec=Data(:,6);
ElpsdayfromTohoku=Data(:,12);
% Elpsday=ElpsdayfromTohoku;
Elpsday=ElpsdayfromTohoku-ElpsdayfromTohoku(1);

Lon=Data(:,7);
Lat=Data(:,8);
Dep=Data(:,9);
refLon=Lon(1);refLat=Lat(1);refDep=Dep(1);
Mag=Data(:,10);
MagID = Mag<9.9;


% ---------[Read cluster catalog]------------------
filename = ['../Data/AizuE.out';'../Data/AizuN.out';'../Data/AizuS.out';'../Data/AizuW.out'; ]
delimiterIn = ' ';
headerlinesIn = 0;

% ---------[Read catalog cluster1]------------------
DataCL1 = importdata(filename(1,:),delimiterIn,headerlinesIn);
yearCL1=DataCL1(:,1);monthCL1=DataCL1(:,2);dayCL1=DataCL1(:,3);
hourCL1=DataCL1(:,4);minutCL1=DataCL1(:,5);secCL1=DataCL1(:,6);
ElpsdayfromTohokuCL1=DataCL1(:,12);
ElpsdayCL1=ElpsdayfromTohokuCL1-ElpsdayfromTohokuCL1(1);
LonCL1=DataCL1(:,7);LatCL1=DataCL1(:,8);DepCL1=DataCL1(:,9);MagCL1=DataCL1(:,10);
MagCL1ID = MagCL1<9.9;
CL1num=numel(LonCL1)
% [CL1x,CL1y,CL1z]=geodetic2ecef(wgs84,LatCL1,LonCL1,DepCL1*(-1));
% [CL1az,CL1elev,CL1slantRange] = geodetic2aer(LatCL1(1),LonCL1(1),DepCL1(1)*(-1),LatCL1,LonCL1,DepCL1*(-1),wgs84);
% [xtCL1,ytCL1,ztCL1] = geodetic2enu(LatCL1,LonCL1,DepCL1*(-1),LatCL1(1),LonCL1(1),DepCL1(1)*(-1),wgs84);

% ---------[Read catalog cluster2]------------------
DataCL2 = importdata(filename(2,:),delimiterIn,headerlinesIn);
yearCL2=DataCL2(:,1);monthCL2=DataCL2(:,2);dayCL2=DataCL2(:,3);
hourCL2=DataCL2(:,4);minutCL2=DataCL2(:,5);secCL2=DataCL2(:,6);
ElpsdayfromTohokuCL2=DataCL2(:,12);
ElpsdayCL2=ElpsdayfromTohokuCL2-ElpsdayfromTohokuCL2(1);
LonCL2=DataCL2(:,7);LatCL2=DataCL2(:,8);DepCL2=DataCL2(:,9);MagCL2=DataCL2(:,10);
MagCL2ID = MagCL2<9.9;
CL2num=numel(LonCL2)
% [CL2x,CL2y,CL2z]=geodetic2ecef(wgs84,LatCL2,LonCL2,DepCL2*(-1));
% [CL2az,CL2elev,CL2slantRange] = geodetic2aer(LatCL2(1),LonCL2(1),DepCL2(1)*(-1),LatCL2,LonCL2,DepCL2*(-1),wgs84);
% [xtCL2,ytCL2,ztCL2] = geodetic2enu(LatCL2,LonCL2,DepCL2*(-1),LatCL2(1),LonCL2(1),DepCL2(1)*(-1),wgs84);

% ---------[Read catalog cluster3]------------------
DataCL3 = importdata(filename(3,:),delimiterIn,headerlinesIn);
yearCL3=DataCL3(:,1);monthCL3=DataCL3(:,2);dayCL3=DataCL3(:,3);
hourCL3=DataCL3(:,4);minutCL3=DataCL3(:,5);secCL3=DataCL3(:,6);
ElpsdayfromTohokuCL3=DataCL3(:,12);
ElpsdayCL3=ElpsdayfromTohokuCL3-ElpsdayfromTohokuCL3(1);
LonCL3=DataCL3(:,7);LatCL3=DataCL3(:,8);DepCL3=DataCL3(:,9);MagCL3=DataCL3(:,10);
MagCL3ID = MagCL3<9.9;
CL3num=numel(LonCL3)
% [CL3x,CL3y,CL3z]=geodetic2ecef(wgs84,LatCL3,LonCL3,DepCL3*(-1));
% [CL3az,CL3elev,CL3slantRange] = geodetic2aer(LatCL3(1),LonCL3(1),DepCL3(1)*(-1),LatCL3,LonCL3,DepCL3*(-1),wgs84);
% [xtCL3,ytCL3,ztCL3] = geodetic2enu(LatCL3,LonCL3,DepCL3*(-1),LatCL3(1),LonCL3(1),DepCL3(1)*(-1),wgs84);

% ---------[Read catalog cluster4]------------------
DataCL4 = importdata(filename(4,:),delimiterIn,headerlinesIn);
yearCL4=DataCL4(:,1);monthCL4=DataCL4(:,2);dayCL4=DataCL4(:,3);
hourCL4=DataCL4(:,4);minutCL4=DataCL4(:,5);secCL4=DataCL4(:,6);
ElpsdayfromTohokuCL4=DataCL4(:,12);
ElpsdayCL4=ElpsdayfromTohokuCL4-ElpsdayfromTohokuCL4(1);
LonCL4=DataCL4(:,7);LatCL4=DataCL4(:,8);DepCL4=DataCL4(:,9);MagCL4=DataCL4(:,10);
MagCL4ID = MagCL4<9.9;
CL4num=numel(LonCL4)
% [CL4x,CL4y,CL4z]=geodetic2ecef(wgs84,LatCL4,LonCL4,DepCL4*(-1));
% [CL4az,CL4elev,CL4slantRange] = geodetic2aer(LatCL4(1),LonCL4(1),DepCL4(1)*(-1),LatCL4,LonCL4,DepCL4*(-1),wgs84);
% [xtCL4,ytCL4,ztCL4] = geodetic2enu(LatCL4,LonCL4,DepCL4*(-1),LatCL4(1),LonCL4(1),DepCL4(1)*(-1),wgs84);


%% GR Parameters 
% GR parameters from Yoshida et al., 2017 (JGR)
Mc = 2.0
% Mc = Mc_maxcurv_corrected
dM = 0.1

%% Seismogenic index setting
% SI = [-10 -5 -1 0 0.1 0.5 1];
% SI = [-1:0.2:1];
% SI = [-3:0.5:1];
SI = [-2:0.5:1];



%% Seismogenic index Fluid Inversion
timestep = [10:10:800];
for i=1:numel(timestep)
    % whole catalog
    ID=find(Mag<9.9 & Mag>=Mc & ElpsdayfromTohoku<=i*10);
    GRnum(i) =numel(ID);
    if GRnum(i)>100
        bvalue(i)=1/log(10)/(mean( Mag(ID) )-(Mc-dM/2));
        avalue(i)=calc_avalue(Mc,bvalue(i),GRnum(i));
    end
    clear ID
    
    % Cluster1
    ID=find(MagCL1<9.9 & MagCL1>=Mc & ElpsdayfromTohokuCL1<=i*10);
    GRnumCL1(i) =numel(ID);
    avalueCL1(i)=0;
    if GRnumCL1(i)>100
        bvalueCL1(i)=1/log(10)/(mean( MagCL1(ID) )-(Mc-dM/2));
        avalueCL1(i)=calc_avalue(Mc,bvalueCL1(i),GRnumCL1(i));
    end
    clear ID
    
    % Cluster2
    ID=find(MagCL2<9.9 & MagCL2>=Mc & ElpsdayfromTohokuCL2<=i*10);
    GRnumCL2(i) =numel(ID);
    avalueCL2(i)=0;
    if GRnumCL2(i)>100
        bvalueCL2(i)=1/log(10)/(mean( MagCL2(ID) )-(Mc-dM/2));
        avalueCL2(i)=calc_avalue(Mc,bvalueCL2(i),GRnumCL2(i));
    end
    clear ID
    
    % Cluster3
    ID=find(MagCL3<9.9 & MagCL3>=Mc & ElpsdayfromTohokuCL3<=i*10);
    GRnumCL3(i) =numel(ID);
    avalueCL3(i)=0;
    if GRnumCL3(i)>100
        bvalueCL3(i)=1/log(10)/(mean( MagCL3(ID) )-(Mc-dM/2));
        avalueCL3(i)=calc_avalue(Mc,bvalueCL3(i),GRnumCL3(i));
    end
    clear ID

    % Cluster4    
    ID=find(MagCL4<9.9 & MagCL4>=Mc & ElpsdayfromTohokuCL4<=i*10);
    GRnumCL4(i) =numel(ID);
    avalueCL4(i)=0;
    if GRnumCL4(i)>100
        bvalueCL4(i)=1/log(10)/(mean( MagCL4(ID) )-(Mc-dM/2));
        avalueCL4(i)=calc_avalue(Mc,bvalueCL4(i),GRnumCL4(i));
    end
    clear ID

    for j=1:numel(SI)
        Vall(i,j) = 10 ^ (avalue(i)-SI(j));
        VCL1(i,j) = 10 ^ (avalueCL1(i)-SI(j));
        VCL2(i,j) = 10 ^ (avalueCL2(i)-SI(j));
        VCL3(i,j) = 10 ^ (avalueCL3(i)-SI(j));
        VCL4(i,j) = 10 ^ (avalueCL4(i)-SI(j));
    end
    clear ID
end


%% ---[Data output]---
flag=0;
% flag=1 --> .mat format, flag=2 --> .dat (ASCII) format, flag=3 --> .nc format,  
if flag==1
    tickday=timestep;
    save('FluidInv_SI_Mc2_v4.mat','tickday','Vall','VCL1','VCL2','VCL3','VCL4');

elseif flag==2
    fileID = fopen('FluidInvResult_SeismogenicIndex.dat','w');
    fileID2 = fopen('FluidInvResult_SeismogenicIndex_clusterE.dat','w');
    fileID3 = fopen('FluidInvResult_SeismogenicIndex_clusterN.dat','w');
    fileID4 = fopen('FluidInvResult_SeismogenicIndex_clusterS.dat','w');
    fileID5 = fopen('FluidInvResult_SeismogenicIndex_clusterW.dat','w');
    
    fprintf(fileID,'%20s  %20s%20s %20s%20s%20s%20s%20s%20s%20s\n',...
        'Days from Tohoku event','GR a-value','GR b-value',...
        'FluidV SI:-2.0(m3)','FluidV SI:-1.5(m3)','FluidV SI:-1.0(m3)','FluidV SI:-0.5(m3)','FluidV SI:0.0(m3)',...
        'FluidV SI:0.5(m3)','FluidV SI:1.0(m3)');
    fprintf(fileID2,'%20s  %20s%20s %20s%20s%20s%20s%20s%20s%20s\n',...
        'Days from Tohoku event','GR a-value','GR b-value',...
        'FluidV SI:-2.0(m3)','FluidV SI:-1.5(m3)','FluidV SI:-1.0(m3)','FluidV SI:-0.5(m3)','FluidV SI:0.0(m3)',...
        'FluidV SI:0.5(m3)','FluidV SI:1.0(m3)');
    fprintf(fileID3,'%20s  %20s%20s %20s%20s%20s%20s%20s%20s%20s\n',...
        'Days from Tohoku event','GR a-value','GR b-value',...
        'FluidV SI:-2.0(m3)','FluidV SI:-1.5(m3)','FluidV SI:-1.0(m3)','FluidV SI:-0.5(m3)','FluidV SI:0.0(m3)',...
        'FluidV SI:0.5(m3)','FluidV SI:1.0(m3)');
    fprintf(fileID4,'%20s  %20s%20s %20s%20s%20s%20s%20s%20s%20s\n',...
        'Days from Tohoku event','GR a-value','GR b-value',...
        'FluidV SI:-2.0(m3)','FluidV SI:-1.5(m3)','FluidV SI:-1.0(m3)','FluidV SI:-0.5(m3)','FluidV SI:0.0(m3)',...
        'FluidV SI:0.5(m3)','FluidV SI:1.0(m3)');
    fprintf(fileID5,'%20s  %20s%20s %20s%20s%20s%20s%20s%20s%20s\n',...
        'Days from Tohoku event','GR a-value','GR b-value',...
        'FluidV SI:-2.0(m3)','FluidV SI:-1.5(m3)','FluidV SI:-1.0(m3)','FluidV SI:-0.5(m3)','FluidV SI:0.0(m3)',...
        'FluidV SI:0.5(m3)','FluidV SI:1.0(m3)');
    for i=1:numel(timestep)
        fprintf(fileID,'%20.3f  %20.3f%20.3f %20.3f%20.3f%20.3f%20.3f%20.3f%20.3f%20.3f\n',...
            timestep(i),avalue(i),bvalue(i),Vall(i,1),Vall(i,2),Vall(i,3),Vall(i,4),Vall(i,5),Vall(i,6),Vall(i,7));
        fprintf(fileID2,'%20.3f  %20.3f%20.3f %20.3f%20.3f%20.3f%20.3f%20.3f%20.3f%20.3f\n',...
            timestep(i),avalueCL1(i),bvalueCL1(i),VCL1(i,1),VCL1(i,2),VCL1(i,3),VCL1(i,4),VCL1(i,5),VCL1(i,6),VCL1(i,7));
        fprintf(fileID3,'%20.3f  %20.3f%20.3f %20.3f%20.3f%20.3f%20.3f%20.3f%20.3f%20.3f\n',...
            timestep(i),avalueCL2(i),bvalueCL2(i),VCL2(i,1),VCL2(i,2),VCL2(i,3),VCL2(i,4),VCL2(i,5),VCL2(i,6),VCL2(i,7));
        fprintf(fileID4,'%20.3f  %20.3f%20.3f %20.3f%20.3f%20.3f%20.3f%20.3f%20.3f%20.3f\n',...
            timestep(i),avalueCL3(i),bvalueCL3(i),VCL3(i,1),VCL3(i,2),VCL3(i,3),VCL3(i,4),VCL3(i,5),VCL3(i,6),VCL3(i,7));
        fprintf(fileID5,'%20.3f  %20.3f%20.3f %20.3f%20.3f%20.3f%20.3f%20.3f%20.3f%20.3f\n',...
            timestep(i),avalueCL4(i),bvalueCL4(i),VCL4(i,1),VCL4(i,2),VCL4(i,3),VCL4(i,4),VCL4(i,5),VCL4(i,6),VCL4(i,7));
    end
    fclose(fileID);fclose(fileID2);fclose(fileID3);fclose(fileID4);fclose(fileID5);

elseif flag==3
    nccreate('FluidInvResult_SI.nc','ElapsedTime(day)','Dimensions',{'Time',80});
    ncwrite('FluidInvResult_SI.nc','ElapsedTime(day)',timestep(1:80));
    
    nccreate('FluidInvResult_SI.nc','a-value','Dimensions',{'Time',80});
    ncwrite ('FluidInvResult_SI.nc','a-value',avalue);
    nccreate('FluidInvResult_SI.nc','b-value','Dimensions',{'Time',80});
    ncwrite ('FluidInvResult_SI.nc','b-value',bvalue);
    nccreate('FluidInvResult_SI.nc','Vall','Dimensions',{'Time',80,'SI',numel(SI)});
    ncwrite ('FluidInvResult_SI.nc','Vall',Vall);
    
    nccreate('FluidInvResult_SI.nc','a-valueCL','Dimensions',{'Time',80,'CL',4});
    ncwrite ('FluidInvResult_SI.nc','a-valueCL',[avalueCL1;avalueCL2;avalueCL3;avalueCL4]');
    nccreate('FluidInvResult_SI.nc','b-valueCL','Dimensions',{'Time',80,'CL',4});
    ncwrite ('FluidInvResult_SI.nc','b-valueCL',[bvalueCL1;bvalueCL2;bvalueCL3;bvalueCL4]');
    
    nccreate('FluidInvResult_SI.nc','VCL','Dimensions',{'Time',80,'SI',numel(SI),'CL',4,});
    VCL=zeros(80,numel(SI),4);
    VCL(:,:,1)=VCL1;VCL(:,:,2)=VCL2;VCL(:,:,3)=VCL3;VCL(:,:,4)=VCL4;
    ncwrite ('FluidInvResult_SI.nc','VCL',VCL);

    ncdisp('FluidInvResult_SI.nc');
end
% return


%% ---------[Plot]------------------
% --- Fluid volume for whole catalog ---
figure
subplot(5,1,1)
plot(ElpsdayfromTohoku(MagID),Mag(MagID),'.','MarkerSize',5,'Color',cm(1,:));
xlim([0 750]);
ylabel('Magnitude');

subplot(5,1,2)
plot(timestep,avalue,'.','MarkerSize',15,'Color',cm(1,:));hold on;
xlim([0 750]);
ylabel('a-value');

subplot(5,1,3)
plot(timestep,bvalue,'.','MarkerSize',15,'Color',cm(2,:));hold on;
xlim([0 750]);
ylabel('b-value');

subplot(5,1,[4 5])
cmtemp=colormap(copper(numel(SI)));
for j=1:numel(SI)
    plot(timestep,Vall(:,j),'.-','MarkerSize',10,'Color',cmtemp(j,:));hold on;
end
xlabel('Days after Tohoku event');
ylabel('Fluid volume (m^3)');
xlim([0 750]);
% ylim([10^4 10^7]);
set(gca,'YScale','log');

c=colorbar('eastoutside');
caxis([-2 1]);
c.Label.String = 'Seismogenic Index (-)';


% --- Fluid volume by SI---
figure
plot(timestep,Vall(:,6),'Color',cm(1,:),       'LineWidth',3);hold on;
plot(timestep,VCL1(:,6),'Color',bluesmap(10,:),'LineWidth',2);hold on;
plot(timestep,VCL2(:,6),'Color',redsmap(10,:) ,'LineWidth',2);hold on;
plot(timestep,VCL3(:,6),'Color',grnsmap(10,:) ,'LineWidth',2);hold on;
plot(timestep,VCL4(:,6),'Color',prplsmap(10,:),'LineWidth',2);hold on;
plot(timestep,VCL1(:,6)+VCL2(:,6)+VCL3(:,6)+VCL4(:,6),'Color',cm(1,:),'LineStyle','--','LineWidth',1);hold on;
xlabel('Days after Tohoku event');
ylabel('Fluid Volume (m^3)');box on;
title("Fluid volume inverted with Seismogenic index");
legend('All events','Cluster E','Cluster N','Cluster S','Cluster W','Sum of all clusters');
set(gca,'YScale','log');
ylim([10^3 10^7]);
title(["Seismogenic Index: ",num2str(SI(6))]);
hold on;


% --- Fluid volume for each cluster ---
figure
subplot(8,1,1)
plot(ElpsdayfromTohokuCL1(MagCL1ID),MagCL1(MagCL1ID),'.','MarkerSize',5,'Color',bluesmap(10,:));
hold on;
xlim([0 750]);ylabel('Magnitude');
subplot(8,1,2)
plot(ElpsdayfromTohokuCL2(MagCL2ID),MagCL2(MagCL2ID),'.','MarkerSize',5,'Color',redsmap(10,:));
hold on;
xlim([0 750]);ylabel('Magnitude');
subplot(8,1,3)
plot(ElpsdayfromTohokuCL3(MagCL3ID),MagCL3(MagCL3ID),'.','MarkerSize',5,'Color',grnsmap(10,:));
hold on;
xlim([0 750]);ylabel('Magnitude');
subplot(8,1,4)
plot(ElpsdayfromTohokuCL4(MagCL4ID),MagCL4(MagCL4ID),'.','MarkerSize',5,'Color',prplsmap(10,:));
hold on;
xlim([0 750]);
ylabel('Magnitude');

subplot(8,1,5)
plot(timestep,avalueCL1,'.','MarkerSize',7,'Color',bluesmap(10,:));hold on;
plot(timestep,avalueCL2,'.','MarkerSize',7,'Color',redsmap(10,:));hold on;
plot(timestep,avalueCL3,'.','MarkerSize',7,'Color',grnsmap(10,:));hold on;
plot(timestep,avalueCL4,'.','MarkerSize',7,'Color',prplsmap(10,:));hold on;
ylabel('a-value');
xlim([0 750]);ylim([3.5 6.5]);

subplot(8,1,6)
plot(timestep,bvalueCL1,'.','MarkerSize',7,'Color',bluesmap(10,:));hold on;
plot(timestep,bvalueCL2,'.','MarkerSize',7,'Color',redsmap(10,:));hold on;
plot(timestep,bvalueCL3,'.','MarkerSize',7,'Color',grnsmap(10,:));hold on;
plot(timestep,bvalueCL4,'.','MarkerSize',7,'Color',prplsmap(10,:));hold on;
ylabel('b-value');
xlim([0 750]);ylim([0.7 1.8]);

subplot(4,1,4)
for j=numel(SI)-1:numel(SI)-1
    plot(timestep,VCL1(:,j),'.-','MarkerSize',10,'Color',bluesmap(10,:));hold on;
    plot(timestep,VCL2(:,j),'.-','MarkerSize',10,'Color',redsmap(10,:));hold on;
    plot(timestep,VCL3(:,j),'.-','MarkerSize',10,'Color',grnsmap(10,:));hold on;
    plot(timestep,VCL4(:,j),'.-','MarkerSize',10,'Color',prplsmap(10,:));hold on;
end
xlabel('Days after Tohoku event');
ylabel('Fluid volume (m^3)');
xlim([0 750]);

set(gca,'YScale','log');

%% Function
% Following part is taken from Zmap, calc_MaxLikelihoodA.m ----              
function avalue=calc_avalue(Mc,bvalue,GRnum)
        % Set up the magnitude range for computing the loglikelihood
        fMc = Mc;
        vMagnitudes = fMc:0.1:10.1;  % 10 is the maximum magnitude. Add an additional bin for diff
        % Compute the unscaled number of expected events for a given b-value
        fBValue = bvalue;
        vExpectation = 10.^(-fBValue * vMagnitudes);
        vExpectation = -diff(vExpectation);
        fSum = sum(vExpectation);
        % Get the a-value as the maximum likelihood solution
        nNumber = GRnum;
        fAValue = log10(nNumber/fSum);
        avalue = fAValue;
end