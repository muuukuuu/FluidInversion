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

a=81;
bluesmap=brewermap(a, 'Blues');
redsmap =brewermap(a, 'Reds');
grnsmap =brewermap(a, 'Greens');
prplsmap=brewermap(a, 'Purples');

alphaindex=0.5;

clr_prl=parula();

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
EVnum=numel(Lon)

[x,y,z]=geodetic2ecef(wgs84,Lat,Lon,Dep*(-1));
[az,elev,slantRange] = geodetic2aer(Lat(1),Lon(1),Dep(1)*(-1),Lat,Lon,Dep*(-1),wgs84);
[xt,yt,zt] = geodetic2enu(Lat,Lon,Dep*(-1),Lat(1),Lon(1),Dep(1)*(-1),wgs84);
Devents=(slantRange*1000).^2./(4*pi*Elpsday*24*3600);% Day-->Sec 

% ---------[Calc Diffusivity, permeability]------------------
% percent=96;
% [BestD]=Find_bestD(Devents,percent);
% Bestk=diff2perm(BestD);


% ---------[Read catalog]------------------
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
[CL1x,CL1y,CL1z]=geodetic2ecef(wgs84,LatCL1,LonCL1,DepCL1*(-1));
[CL1az,CL1elev,CL1slantRange] = geodetic2aer(LatCL1(1),LonCL1(1),DepCL1(1)*(-1),LatCL1,LonCL1,DepCL1*(-1),wgs84);
[xtCL1,ytCL1,ztCL1] = geodetic2enu(LatCL1,LonCL1,DepCL1*(-1),LatCL1(1),LonCL1(1),DepCL1(1)*(-1),wgs84);
refLonCL1=mean(LonCL1(1:10));refLatCL1=mean(LatCL1(1:10));refDepCL1=mean(DepCL1(1:10));

DeventsCL1=(CL1slantRange*1000).^2./(4*pi*ElpsdayCL1*24*3600);% Day-->Sec 
DeventsCL1(1)=0;


% ---------[Read catalog cluster2]------------------
DataCL2 = importdata(filename(2,:),delimiterIn,headerlinesIn);
yearCL2=DataCL2(:,1);monthCL2=DataCL2(:,2);dayCL2=DataCL2(:,3);
hourCL2=DataCL2(:,4);minutCL2=DataCL2(:,5);secCL2=DataCL2(:,6);
ElpsdayfromTohokuCL2=DataCL2(:,12);
ElpsdayCL2=ElpsdayfromTohokuCL2-ElpsdayfromTohokuCL2(1);
LonCL2=DataCL2(:,7);LatCL2=DataCL2(:,8);DepCL2=DataCL2(:,9);MagCL2=DataCL2(:,10);
MagCL2ID = MagCL2<9.9;
CL2num=numel(LonCL2)
[CL2x,CL2y,CL2z]=geodetic2ecef(wgs84,LatCL2,LonCL2,DepCL2*(-1));
[CL2az,CL2elev,CL2slantRange] = geodetic2aer(LatCL2(1),LonCL2(1),DepCL2(1)*(-1),LatCL2,LonCL2,DepCL2*(-1),wgs84);
[xtCL2,ytCL2,ztCL2] = geodetic2enu(LatCL2,LonCL2,DepCL2*(-1),LatCL2(1),LonCL2(1),DepCL2(1)*(-1),wgs84);

DeventsCL2=(CL2slantRange*1000).^2./(4*pi*ElpsdayCL2*24*3600);% Day-->Sec 
DeventsCL2(1)=0;


% ---------[Read catalog cluster3]------------------
DataCL3 = importdata(filename(3,:),delimiterIn,headerlinesIn);
yearCL3=DataCL3(:,1);monthCL3=DataCL3(:,2);dayCL3=DataCL3(:,3);
hourCL3=DataCL3(:,4);minutCL3=DataCL3(:,5);secCL3=DataCL3(:,6);
ElpsdayfromTohokuCL3=DataCL3(:,12);
ElpsdayCL3=ElpsdayfromTohokuCL3-ElpsdayfromTohokuCL3(1);
LonCL3=DataCL3(:,7);LatCL3=DataCL3(:,8);DepCL3=DataCL3(:,9);MagCL3=DataCL3(:,10);
MagCL3ID = MagCL3<9.9;
CL3num=numel(LonCL3)
[CL3x,CL3y,CL3z]=geodetic2ecef(wgs84,LatCL3,LonCL3,DepCL3*(-1));
[CL3az,CL3elev,CL3slantRange] = geodetic2aer(LatCL3(1),LonCL3(1),DepCL3(1)*(-1),LatCL3,LonCL3,DepCL3*(-1),wgs84);
[xtCL3,ytCL3,ztCL3] = geodetic2enu(LatCL3,LonCL3,DepCL3*(-1),LatCL3(1),LonCL3(1),DepCL3(1)*(-1),wgs84);

DeventsCL3=(CL3slantRange*1000).^2./(4*pi*ElpsdayCL3*24*3600);% Day-->Sec 
DeventsCL3(1)=0;


% ---------[Read catalog cluster4]------------------
DataCL4 = importdata(filename(4,:),delimiterIn,headerlinesIn);
yearCL4=DataCL4(:,1);monthCL4=DataCL4(:,2);dayCL4=DataCL4(:,3);
hourCL4=DataCL4(:,4);minutCL4=DataCL4(:,5);secCL4=DataCL4(:,6);
ElpsdayfromTohokuCL4=DataCL4(:,12);
ElpsdayCL4=ElpsdayfromTohokuCL4-ElpsdayfromTohokuCL4(1);
LonCL4=DataCL4(:,7);LatCL4=DataCL4(:,8);DepCL4=DataCL4(:,9);MagCL4=DataCL4(:,10);
MagCL4ID = MagCL4<9.9;
CL4num=numel(LonCL4)
[CL4x,CL4y,CL4z]=geodetic2ecef(wgs84,LatCL4,LonCL4,DepCL4*(-1));
[CL4az,CL4elev,CL4slantRange] = geodetic2aer(LatCL4(1),LonCL4(1),DepCL4(1)*(-1),LatCL4,LonCL4,DepCL4*(-1),wgs84);
[xtCL4,ytCL4,ztCL4] = geodetic2enu(LatCL4,LonCL4,DepCL4*(-1),LatCL4(1),LonCL4(1),DepCL4(1)*(-1),wgs84);

DeventsCL4=(CL4slantRange*1000).^2./(4*pi*ElpsdayCL4*24*3600);% Day-->Sec 
DeventsCL4(1)=0;


%% ---[Time variable Diffusivity estimation & Plot]-------------------------------
RefD=[0.1 0.2 0.5 1 2];

% Analysis condition setting
timestep=[0:10:800];
tickday=timestep;

temp=diff(timestep);
window=temp(1);


% Loop for four Clusters
for j=1:4
    if j==1
        ElpsdayTemp=ElpsdayCL1;
        TempslantRange=CL1slantRange;
        DeventsTemp=DeventsCL1;
        xtTemp=xtCL1;ytTemp=ytCL1;DepTemp=DepCL1;
        cmap_temp=bluesmap(a,:);
        temp=ElpsdayfromTohoku(1);
    elseif j==2
        ElpsdayTemp=ElpsdayCL2;
        TempslantRange=CL2slantRange;
        DeventsTemp=DeventsCL2;
        xtTemp=xtCL2;ytTemp=ytCL2;DepTemp=DepCL2;
        cmap_temp=redsmap(a,:);
        temp=ElpsdayfromTohoku(2);
    elseif j==3
        ElpsdayTemp=ElpsdayCL3;
        TempslantRange=CL3slantRange;
        DeventsTemp=DeventsCL3;
        xtTemp=xtCL3;ytTemp=ytCL3;DepTemp=DepCL3;
        cmap_temp=grnsmap(a,:);
        temp=ElpsdayfromTohoku(3);
    elseif j==4
        ElpsdayTemp=ElpsdayCL4;
        TempslantRange=CL4slantRange;
        DeventsTemp=DeventsCL4;
        xtTemp=xtCL4;ytTemp=ytCL4;DepTemp=DepCL4;
        cmap_temp=prplsmap(a,:);
        temp=ElpsdayfromTohoku(4);
    end
    
    % Loop for timestep
    for i=2:numel(timestep)
        tickday(i)=timestep(i);
        ID=find(ElpsdayTemp<=tickday(i));        
        
        % Cluster dimension analysis
        lengthAxis(j,i,1:3)=[0 0 0];
        if numel(ID)>=3
            src1=[xtTemp(ID),ytTemp(ID),DepTemp(ID)];
            [V,D,M]=shuseibun_num(src1);
            Ds(1:3)=[D(1,1),D(2,2),D(3,3)];
            G(j,i,:)=M;
            lengthAxis(j,i,1:3)=2*sqrt(Ds);
            [azimuth,elevation]=cart2sph(V(1,:),V(2,:),V(3,:));
            azim_cart(j,i,:) = azimuth;
            elv_cart(j,i,:)  = elevation;
        end
        clear ID
    end    
end


%% ---[Inversion with cubic law ]---
%  --- Pore pressure distribution setting ---  
% Pore pressure at outlet
Pout=[0:10:260];
% Pore pressure at inlet
IniPin=270;

% Most plausible scinario
targetPP=170;

% --- Plot pore pressure distribution in time series ----------------------
figure
subplot(2,1,1)
for i=1:numel(Pout)
    temp=round(i/(numel(Pout))*256)-1;
    Pinlet_time(i,:)=IniPin-(IniPin-Pout(i))/max(tickday).*tickday;
    dP_time(i,:)=Pinlet_time(i,:)-Pout(i);
        
    plot(tickday,Pinlet_time(i,:),'Color',clr_prl(temp,:));hold on;
    if dP_time(i,1)==200
        plot(tickday,Pinlet_time(i,:),'Color',clr_prl(temp,:),'LineWidth',3);hold on;
    end
end
xlim([0 800]);
xlabel('Days after Tohoku event');
ylabel('Pressure (MPa)');
title('inlet pressure (MPa)')

subplot(2,1,2)
for i=1:numel(Pout)
    temp=round(i/(numel(Pout))*256)-1;
    plot(tickday,dP_time(i,:),'Color',clr_prl(temp,:));hold on;
    if dP_time(i,1)==200
        plot(tickday,dP_time(i,:),'Color',clr_prl(temp,:),'LineWidth',3);hold on;
    end
end
xlim([0 800]);
xlabel('Days after Tohoku event');
ylabel('Pressure (MPa)');
title('Pressure difference (MPa)')


% ---[Flow zone definition]------------------------------------------------
w(1,:)=lengthAxis(1,:,3)*2;
L(1,:)=lengthAxis(1,:,2)*2; % km

w(2,:)=lengthAxis(2,:,3)*2;
L(2,:)=lengthAxis(2,:,2)*2; % km

w(3,:)=lengthAxis(3,:,2)*2;
L(3,:)=lengthAxis(3,:,1)*2; % km

w(4,:)=lengthAxis(4,:,2)*2;
L(4,:)=lengthAxis(4,:,3)*2; % km


% ---[Hydraulic aperture to permeability]---
d=0.1;    %[mm] Norbeck et al., 
% 59.	Norbeck, J. H. & Shelly, D. R. Exploring the Role of Mixed-Mechanism Fracturing and Fluid-Faulting Interactions During the 2014 Long Valley Caldera, California, Earthquake Swarm. 
% Proc. Geotherm. Reserv. Eng. Stanford Univ. 1–16 (2018).
% 60.	Watanabe, N., Hirano, N. & Tsuchiya, N. Diversity of channeling flow in heterogeneous aperture distribution inferred from integrated experimental-numerical analysis on flow through shear fracture in granite. 
% J. Geophys. Res. Solid Earth 114, 1–17 (2009). 10.1029/2008JB005959

d=d/10^3; %[m]
k_d=d^2/12; %[m^2]

ratio=0.2;    % ratio of flow path to the fault area
% 61.	Ishibashi, T., Watanabe, N., Hirano, N., Okamoto, A. & Tsuchiya, N. Beyond‐laboratory‐scale prediction for channeling flows through subsurface rock fractures with heterogeneous aperture distributions revealed by laboratory evaluation. 
% J. Geophys. Res.. Solid Earth 120, 106–124 (2015). 10.1002/2014JB011555

% % [Q and V calc]_______________________________________________________________________
for j=1:4 %Cluster   
    for i=1:numel(tickday)
        for k=1:numel(Pout)
            Q(j,i,k)=CalcQ_Cubic(d,w(j,i)*10^3,L(j,i)*10^3,dP_time(k,i)*10^6,5)*ratio;
            V(j,i,k)=Q(j,i,k)* window *24*60*60; % 100 day-->sec
            if i==1
                Vsum(j,i,k)=V(j,i,k);
            else
                Vsum(j,i,k)=Vsum(j,i-1,k)+V(j,i,k);
            end
        end
    end
end
VCL1(:,:)=Vsum(1,:,:);VCL2(:,:)=Vsum(2,:,:);VCL3(:,:)=Vsum(3,:,:);VCL4(:,:)=Vsum(4,:,:);
% % _______________________________________________________________________


figure
subplot(2,1,1)
a=40;
k=11; % --> dP=170 MPa
for j=1:4 %Cluster
    if j==1
        cmap_temp=bluesmap(:,:);
%         plot(tickday,Q_cl1*60,'--','Color',cmap_temp(a,:),'LineWidth',2);hold on;
    elseif j==2
        cmap_temp=redsmap(:,:);
    elseif j==3
        cmap_temp=grnsmap(:,:);
    elseif j==4
        cmap_temp=prplsmap(:,:);
    end
    plot(tickday,Q(j,:,k)*60,'Color',cmap_temp(a,:),'LineWidth',2);hold on;
end
xlim([0 800]);
xlabel('Days after Tohoku event');
ylabel('Flow rate (m^3/min)');
%     set(gca,'YScale','log');
legend('Cluster E','Cluster N','Cluster S','Cluster W');
title(['Pore pressure difference:',num2str(dP_time(k,1)),' MPa']);
set(gca,'YScale','log');

subplot(2,1,2)
% subplot(3,1,2)
for kk=1:numel(Pout)
    for i=1:numel(tickday)
        AllV(i,kk)=sum( Vsum(:,i,kk));
    end
end
plot(tickday,AllV(:,k),'Color',cm(1,:),'LineWidth',3);hold on;
for j=1:4 %Cluster
    if j==1
        cmap_temp=bluesmap(:,:);
        plot(tickday,VCL1(:,k),'-','Color',cmap_temp(a,:),'LineWidth',2);hold on;
    elseif j==2
        cmap_temp=redsmap(:,:);
        plot(tickday,VCL2(:,k),'-','Color',cmap_temp(a,:),'LineWidth',2);hold on;
    elseif j==3
        cmap_temp=grnsmap(:,:);
        plot(tickday,VCL3(:,k),'-','Color',cmap_temp(a,:),'LineWidth',2);hold on;
    elseif j==4
        cmap_temp=prplsmap(:,:);
        plot(tickday,VCL4(:,k),'-','Color',cmap_temp(a,:),'LineWidth',2);hold on;
    end    
end
legend('All event','Cluster E','Cluster N','Cluster S','Cluster W');
xlim([0 800]);
%     xlabel('Days from 1st event');
xlabel('Days after Tohoku event');
ylabel('Volume (m^3)');
%     set(gca,'YScale','log');
title(['Pore pressure difference:',num2str(dP_time(k,1)),' MPa']);


%% ---[Data Output]---
flag=0;
% flag=1 --> .mat format, flag=2 --> .dat (ASCII) format, flag=3 --> .nc format,  
if flag==1
    tickday=timestep;
    save('FluidInv_cubic_v4.mat','tickday','AllV','VCL1','VCL2','VCL3','VCL4');

elseif flag==2

elseif flag==3
    tickday=timestep;
    nccreate('FluidInvResult_cubic.nc','ElapsedTime(day)','Dimensions',{'Time',numel(timestep)});
    ncwrite('FluidInvResult_cubic.nc','ElapsedTime(day)',tickday);
    nccreate('FluidInvResult_cubic.nc','Vall','Dimensions',{'Time',numel(timestep),'dP',numel(Pout)});
    ncwrite('FluidInvResult_cubic.nc','Vall',AllV);

    nccreate('FluidInvResult_cubic.nc','VCL','Dimensions',{'time',numel(timestep),'dP',numel(Pout),'CL',4});
    VCL=zeros(numel(timestep),numel(Pout),4);
    VCL(:,:,1)=VCL1;VCL(:,:,2)=VCL2;VCL(:,:,3)=VCL3;VCL(:,:,4)=VCL4;
    ncwrite('FluidInvResult_cubic.nc','VCL',VCL);

    ncdisp('FluidInvResult_cubic.nc');
end

return


%% Local Functions
% Function Cubic's law
function Q=CalcQ_Cubic(d,w,L,dP,flag)
    if flag==1
        myu= 105 * 10^(-6); %300 degC  8 km depth (80 MPa) [Pa*sec]
    elseif flag==2
        myu= 109 * 10^(-6); %300 degC 10 km depth (100 MPa) [Pa*sec]
    elseif flag==3
        myu= 91 * 10^(-6); %350 degC 8 km depth (80 MPa) [Pa*sec]
    elseif flag==4
        myu= 95 * 10^(-6); %350 degC 10 km depth (100 MPa) [Pa*sec]
    elseif flag==5    
        eta=1.24*10^(-4); %Pa s, 
        % IAPWS95 formulation
        myu=eta;
    end
    Q=d^3/12/myu*w/L*dP; % [m3/sec]
    if w==0
        Q=0;
    end
end



%% Function PCA(Shuseibun) %%
function [V,D,M]=shuseibun_num(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  x : ( N x 3 )行列
%  V(:,1)　: 第三固有ベクトル(Pole)
%  V(:,2)  : 第二固有ベクトル
%  V(:,3)  : 第一固有ベクトル
%  D  : [第三、第二、第一固有値]
%  M  : [Mean of X, Mean of Y, Mean of Z]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=mean(x);S=std(x);

Ds=[S(1) 0    0   ; ...
    0    S(2) 0   ; ...
    0    0    S(3)];

cor=Ds*corrcoef(x)*Ds;
[V,D]=eig(cor);

D2=D-D;V2=V-V;
for i=1:3;
%     abs( diag(D) )
%     abs( max(diag(D)) )
  id=find( diag(D)==max( diag(D) ) );
  if (numel(id)==1);
      D2(4-i,4-i)=D(id,id);D(id,id)=0;
      V2(1:3,4-i)=V(1:3,id);
  else
      D=0;
  end;
end;
D=D2;
V=V2;
%if (D(2,2)<D(1,1)) & (D(2,2)<D(3,3));
%   DD=D(1,1);D(1,1)=D(2,2);D(2,2)=DD;
%   VV=V(:,1);V(:,1)=V(:,2);V(:,2)=VV;
%elseif (D(3,3)<D(1,1)) & (D(3,3)<D(2,2));
%   DD=D(1,1);D(1,1)=D(3,3);D(3,3)=DD;
%   VV=V(:,1);V(:,1)=V(:,3);V(:,3)=VV;
%end
clear S cor DD VV Ds N V2 D2
end
