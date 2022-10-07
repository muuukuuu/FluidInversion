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
Elpsday=ElpsdayfromTohoku-ElpsdayfromTohoku(1);

Lon=Data(:,7);
Lat=Data(:,8);
Dep=Data(:,9);
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

% ---------[Read catalog cluster2]------------------
DataCL2 = importdata(filename(2,:),delimiterIn,headerlinesIn);
yearCL2=DataCL2(:,1);monthCL2=DataCL2(:,2);dayCL2=DataCL2(:,3);
hourCL2=DataCL2(:,4);minutCL2=DataCL2(:,5);secCL2=DataCL2(:,6);
ElpsdayfromTohokuCL2=DataCL2(:,12);
ElpsdayCL2=ElpsdayfromTohokuCL2-ElpsdayfromTohokuCL2(1);
LonCL2=DataCL2(:,7);LatCL2=DataCL2(:,8);DepCL2=DataCL2(:,9);MagCL2=DataCL2(:,10);
MagCL2ID = MagCL2<9.9;
CL2num=numel(LonCL2)

% ---------[Read catalog cluster3]------------------
DataCL3 = importdata(filename(3,:),delimiterIn,headerlinesIn);
yearCL3=DataCL3(:,1);monthCL3=DataCL3(:,2);dayCL3=DataCL3(:,3);
hourCL3=DataCL3(:,4);minutCL3=DataCL3(:,5);secCL3=DataCL3(:,6);
ElpsdayfromTohokuCL3=DataCL3(:,12);
ElpsdayCL3=ElpsdayfromTohokuCL3-ElpsdayfromTohokuCL3(1);
LonCL3=DataCL3(:,7);LatCL3=DataCL3(:,8);DepCL3=DataCL3(:,9);MagCL3=DataCL3(:,10);
MagCL3ID = MagCL3<9.9;
CL3num=numel(LonCL3)

% ---------[Read catalog cluster4]------------------
DataCL4 = importdata(filename(4,:),delimiterIn,headerlinesIn);
yearCL4=DataCL4(:,1);monthCL4=DataCL4(:,2);dayCL4=DataCL4(:,3);
hourCL4=DataCL4(:,4);minutCL4=DataCL4(:,5);secCL4=DataCL4(:,6);
ElpsdayfromTohokuCL4=DataCL4(:,12);
ElpsdayCL4=ElpsdayfromTohokuCL4-ElpsdayfromTohokuCL4(1);
LonCL4=DataCL4(:,7);LatCL4=DataCL4(:,8);DepCL4=DataCL4(:,9);MagCL4=DataCL4(:,10);
MagCL4ID = MagCL4<9.9;
CL4num=numel(LonCL4)


%% ---[McGarr inversion based　on McGarr, 2014 JGR]
Mo = Mw2Mo(Mag(MagID));
G  = 3e10;        % Shear modulus: representative value from McGarr 2014
% G  = 3.3e10;      % Shear modulus: representative value from Nakajima
V = sum(Mo)/2/G    % [m3]
% V = max(Mo)/G;     % [m3]
mu=0.6;
% mu=0.85;
% lamda=G;
lamda=26.49;    % [m3]

V2=3/(2*mu*(3*lamda+2*G))*sum(Mo);

%% Figure GR and Mo ratio
Mo_descend = sort(Mo,'descend');
Mag_descend = sort(Mag(MagID),'descend');
labels = num2str(Mag_descend(1:10));

subplot(1,2,1)
dM=0.1;
edges=[-0.6-dM/2:0.1:4.6+dM/2];
h=histogram(Mag(MagID),edges);
h.BinWidth = 0.1;
set(gca,'YScale','log');
xlabel('Magnitude');ylabel('Frequency');
axis square

subplot(1,2,2)
pie(Mo_descend(1:10),labels);


%% McGarr's model Fluid Inversion
timestep=[10:10:800];
for i=1:numel(timestep)
    % whole catalog
    ID=find(ElpsdayfromTohoku<=timestep(i) & Mag<9.9);    
    Mo = Mw2Mo(Mag(ID));
    MoAll=Mo;
    Vall(i)=sum(Mo)/2/G;
    clear ID

    % Cluster1
    ID=find(ElpsdayfromTohokuCL1<=timestep(i) & MagCL1<9.9);    
    Mo = Mw2Mo(MagCL1(ID));
    VCL1(i)=sum(Mo)/2/G;
    clear ID

    % Cluster2    
    ID=find(ElpsdayfromTohokuCL2<=timestep(i) & MagCL2<9.9);    
    Mo = Mw2Mo(MagCL2(ID));
    VCL2(i)=sum(Mo)/2/G;
    clear ID
    
    % Cluster3
    ID=find(ElpsdayfromTohokuCL3<=timestep(i) & MagCL3<9.9);    
    Mo = Mw2Mo(MagCL3(ID));
    VCL3(i)=sum(Mo)/2/G;
    clear ID

    % Cluster4
    ID=find(ElpsdayfromTohokuCL4<=timestep(i) & MagCL4<9.9);    
    Mo = Mw2Mo(MagCL4(ID));
    VCL4(i)=sum(Mo)/2/G;
    clear ID
end

%% ---------[Plot]------------------
figure
plot(timestep,Vall,'Color',cm(1,:),       'LineWidth',3);hold on;
plot(timestep,VCL1,'Color',bluesmap(10,:),'LineWidth',2);hold on;
plot(timestep,VCL2,'Color',redsmap(10,:) ,'LineWidth',2);hold on;
plot(timestep,VCL3,'Color',grnsmap(10,:) ,'LineWidth',2);hold on;
plot(timestep,VCL4,'Color',prplsmap(10,:),'LineWidth',2);hold on;
xlabel('Days after Tohoku event');
ylabel('Fluid Volume (m^3)');box on;
title("Fluid volume inverted with McGarr's theory")
legend('All events','Cluster E','Cluster N','Cluster S','Cluster W');
hold on;
ylim([10^3 10^7]);
set(gca,'YScale','log');


%% ---[Data output]---
flag=0;
% flag=1 --> .mat format, flag=2 --> .dat (ASCII) format, flag=3 --> .nc format,  
if flag==1
    tickday=timestep;
    save('FluidInv_McGarr.mat','tickday','Vall','VCL1','VCL2','VCL3','VCL4')

elseif flag==2
    fileID = fopen('FluidInvResult_McGarr.dat','w');
    OutData=[timestep;Vall;VCL1;VCL2;VCL3;VCL4];
    fprintf(fileID,'%20s %20s %20s %20s %20s %20s\n',...
        'Days from Tohoku event','Entire FluidV(m3)','FluidV clusterE(m3)','FluidV clusterN(m3)','FluidV clusterS(m3)','FluidV clusterW(m3)');
    fprintf(fileID,'%20.3f %20.3f %20.3f %20.3f %20.3f %20.3f\n',OutData);
    fclose(fileID);

elseif flag==3
    nccreate('FluidInvResult_McGarr.nc','ElapsedTime(day)','Dimensions',{'Time',80});
    ncwrite('FluidInvResult_McGarr.nc','ElapsedTime(day)',timestep(1:80));

    nccreate('FluidInvResult_McGarr.nc','Vall(m3)','Dimensions',{'time',80});
    ncwrite('FluidInvResult_McGarr.nc','Vall(m3)',Vall(1:80));

    nccreate('FluidInvResult_McGarr.nc','VCL1(m3)','Dimensions',{'time',80});
    ncwrite('FluidInvResult_McGarr.nc','VCL1(m3)',VCL1(1:80));

    nccreate('FluidInvResult_McGarr.nc','VCL2(m3)','Dimensions',{'time',80});
    ncwrite('FluidInvResult_McGarr.nc','VCL2(m3)',VCL2(1:80));

    nccreate('FluidInvResult_McGarr.nc','VCL3(m3)','Dimensions',{'time',80});
    ncwrite('FluidInvResult_McGarr.nc','VCL3(m3)',VCL3(1:80));

    nccreate('FluidInvResult_McGarr.nc','VCL4(m3)','Dimensions',{'time',80});
    ncwrite('FluidInvResult_McGarr.nc','VCL4(m3)',VCL4(1:80));
end

return

