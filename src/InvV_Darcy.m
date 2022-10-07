% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% 2022/10/06
% Yusuke Mukuhira
% Institue of fluid science, Tohoku University
% tohoku@tohoku.ac.jp
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% To model time decreasing source pore pressure
% Estimate Q, V, and Vall with decreasing Pp
% We estimate D more subjective way
% D for cluster D=5 ~ 70 days
% and W are estimated from 100 days data and used
% deterministic value for later

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
Devents=(slantRange*1000).^2./(4*pi*Elpsday*24*3600);           % Day-->Sec 

% ---------[Calc Diffusivity for whole catalog]------------------
percent=96;
[BestD]=Find_bestD(Devents,percent);


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
[CL1x,CL1y,CL1z]=geodetic2ecef(wgs84,LatCL1,LonCL1,DepCL1*(-1));
[CL1az,CL1elev,CL1slantRange] = geodetic2aer(LatCL1(1),LonCL1(1),DepCL1(1)*(-1),LatCL1,LonCL1,DepCL1*(-1),wgs84);
[xtCL1,ytCL1,ztCL1] = geodetic2enu(LatCL1,LonCL1,DepCL1*(-1),LatCL1(1),LonCL1(1),DepCL1(1)*(-1),wgs84);
refLonCL1=mean(LonCL1(1:10));refLatCL1=mean(LatCL1(1:10));refDepCL1=mean(DepCL1(1:10));

DeventsCL1=(CL1slantRange*1000).^2./(4*pi*ElpsdayCL1*24*3600);% Day-->Sec 
DeventsCL1(1)=0;
[BestDCL1]=Find_bestD(DeventsCL1,percent);


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
[BestDCL2]=Find_bestD(DeventsCL2,percent);


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
[BestDCL3]=Find_bestD(DeventsCL3,percent);


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
[BestDCL4]=Find_bestD(DeventsCL4,percent);


%% ---[Time variable Diffusivity estimation & Plot]-------------------------------
% Reference Diffusivity setting
RefD=[0.1 0.2 0.5 1 2 5];

% Figure axis defined
tempax=figure;

% Analysis condition setting
timestep=[0:10:800];
temp=diff(timestep);
window=temp(1);

a=40;
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
    
    k=1;
    tickday(1)=0;
    figure;
    % Loop for timestep
    for i=2:numel(timestep)
        tickday(i)=timestep(i);
        ID=find(ElpsdayTemp<=tickday(i));
        [BestD_time(j,i)]=Find_bestD(DeventsTemp(ID),percent);
        
        if rem(tickday(i),100)==0
            subplot(3,3,k)
            plot(ElpsdayTemp(ID)+temp,TempslantRange(ID),'.','MarkerSize',5,'Color',cmap_temp);hold on;
            rtheo=sqrt(4*pi*BestD_time(j,i)*ElpsdayTemp*24*3600);hold on;
            rref=sqrt(4*pi.*RefD.*ElpsdayTemp*24*3600);hold on;
            plot(ElpsdayTemp(ID)+temp,rtheo(ID)/1000,'-','LineWidth',2,'Color',cmap_temp);hold on;
            if i~=1 & j==1 || j==4
                plot(ElpsdayTemp(ID)+temp,rref(ID,1)/1000,'--','LineWidth',1,'Color',cm(1,:));hold on;
                plot(ElpsdayTemp(ID)+temp,rref(ID,2)/1000,'--','LineWidth',1,'Color',cm(1,:));hold on;
                plot(ElpsdayTemp(ID)+temp,rref(ID,3)/1000,'--','LineWidth',1,'Color',cm(1,:));hold on;
                plot(ElpsdayTemp(ID)+temp,rref(ID,4)/1000,'--','LineWidth',1,'Color',cm(1,:));hold on;
                plot(ElpsdayTemp(ID)+temp,rref(ID,5)/1000,'--','LineWidth',1,'Color',cm(1,:));hold on;
                plot(ElpsdayTemp(ID)+temp,rref(ID,6)/1000,'--','LineWidth',1,'Color',cm(1,:));hold on;
            end

            title(['Diffusivity = ',num2str(BestD_time(j,i)),' (m^2/sec)']);
            xlabel('Days after Tohoku event');
            ylabel('Distance from 1st event (km)');
            k=k+1;
        end
        ylim([0 10])
        
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
    
    subplot(3,3,9)
    plot(tickday,BestD_time(j,:),'.','MarkerSize',20,'Color',cmap_temp);hold on;
    if j==1
        SuperBestD_time(j,:)=BestD_time(j,:);
        ID=find(tickday<=50);
        SuperBestD_time(j,ID(end)+1:end)=RefD(2);
        plot(tickday,SuperBestD_time(j,:),'.','MarkerSize',20,'Color',cm(2,:));hold on;
    elseif j==4
        SuperBestD_time(j,:)=BestD_time(j,:);
        ID=find(tickday<=50);
        SuperBestD_time(j,ID(end)+1:end)=RefD(2);
        plot(tickday,SuperBestD_time(j,:),'.','MarkerSize',20,'Color',cm(2,:));hold on;
    end
    
    % Plot referenceD r-t plot
    if j==1 | j==4
        line([0,tickday(end)],[RefD(1) RefD(1)],'LineWidth',1,'Color',cm(1,:));hold on;
        line([0,tickday(end)],[RefD(2) RefD(2)],'LineWidth',1,'Color',cm(2,:));hold on;
        line([0,tickday(end)],[RefD(3) RefD(3)],'LineWidth',1,'Color',cm(3,:));hold on;
        line([0,tickday(end)],[RefD(4) RefD(4)],'LineWidth',1,'Color',cm(4,:));hold on;
        line([0,tickday(end)],[RefD(5) RefD(5)],'LineWidth',1,'Color',cm(5,:));hold on;
        line([0,tickday(end)],[RefD(6) RefD(6)],'LineWidth',1,'Color',cm(6,:));hold on;
    end
    xlabel('Days after Tohoku event');
%     xlabel('Days from 1st event');
    ylabel('Diffusivity (m^2/sec)');
    set(gcf,'Position',[159   0   860*0.9   771*0.9]); % By Hotta Dec.18

    
    if true
        figure(tempax);
        subplot(4,4,[4*j-3 4*j-2 4*j-1])
        plot(ElpsdayTemp+temp,TempslantRange,'.','MarkerSize',5,'Color',cmap_temp);hold on;
        rref=sqrt(4*pi.*RefD.*ElpsdayTemp*24*3600);hold on;
        
        plot(ElpsdayTemp+temp,rref(:,1)/1000,'--','LineWidth',0.5,'Color',cm(1,:));hold on;
        plot(ElpsdayTemp+temp,rref(:,2)/1000,'--','LineWidth',0.5,'Color',cm(1,:));hold on;
        plot(ElpsdayTemp+temp,rref(:,3)/1000,'--','LineWidth',0.5,'Color',cm(1,:));hold on;
        plot(ElpsdayTemp+temp,rref(:,4)/1000,'--','LineWidth',0.5,'Color',cm(1,:));hold on;
        plot(ElpsdayTemp+temp,rref(:,5)/1000,'--','LineWidth',0.5,'Color',cm(1,:));hold on;
        plot(ElpsdayTemp+temp,rref(:,6)/1000,'--','LineWidth',0.5,'Color',cm(1,:));hold on;
        ylabel('Distance from 1st event (km)');
                    
        if j==4
            xlabel('Days after Tohoku event');
        end
        ylim([0 10]);
        
        subplot(4,4,4*j)    
        plot(tickday,BestD_time(j,:),'.','MarkerSize',10,'Color',cmap_temp);hold on;
        if j==1
            plot(tickday,SuperBestD_time(j,:),'.','MarkerSize',15,'Color',cm(2,:));hold on;
        end
        if j==4
            plot(tickday,SuperBestD_time(j,:),'.','MarkerSize',15,'Color',cm(2,:));hold on;
            xlabel('Days after Tohoku event');
        end
        ylabel('Diffusivity (m^2/sec)');
    end
end

% Update BestD
BestD_time(1,:)=SuperBestD_time(1,:);
BestD_time(4,:)=SuperBestD_time(4,:);


%% ---[Plot 3D distribution, Diffusivity,and permeability]---
a=81;
for j=1:4
    figure
    if j==1
        ElpsdayTemp=ElpsdayfromTohokuCL1;
        TempslantRange=CL1slantRange;
        DeventsTemp=DeventsCL1;
        xtTemp=xtCL1;ytTemp=ytCL1;DepTemp=DepCL1;
        cmap_temp=bluesmap(:,:);
        temp=ElpsdayfromTohoku(1);
    elseif j==2
        ElpsdayTemp=ElpsdayfromTohokuCL2;
        TempslantRange=CL2slantRange;
        DeventsTemp=DeventsCL2;
        xtTemp=xtCL2;ytTemp=ytCL2;DepTemp=DepCL2;
        cmap_temp=redsmap(:,:);
        temp=ElpsdayfromTohoku(2);
    elseif j==3
        ElpsdayTemp=ElpsdayfromTohokuCL3;
        TempslantRange=CL3slantRange;
        DeventsTemp=DeventsCL3;
        xtTemp=xtCL3;ytTemp=ytCL3;DepTemp=DepCL3;
        cmap_temp=grnsmap(:,:);
        temp=ElpsdayfromTohoku(3);
    elseif j==4
        ElpsdayTemp=ElpsdayfromTohokuCL4;
        TempslantRange=CL4slantRange;
        DeventsTemp=DeventsCL4;
        xtTemp=xtCL4;ytTemp=ytCL4;DepTemp=DepCL4;
        cmap_temp=prplsmap(:,:);
        temp=ElpsdayfromTohoku(4);
    end
    
    % Plot 3D distribution    
    subplot(3,3,[1 2 4 5])
    for i=1:numel(timestep)
        ID=find(ElpsdayTemp<=i*window & ElpsdayTemp>(i-1)*window);
        if ID~=0
            scatter3(xtTemp(ID),ytTemp(ID),DepTemp(ID),10,cmap_temp(a-i,:),'filled');hold on;
        end
        clear ID
    end
    alpha(alphaindex);
   
    c = colorbar('eastoutside');
    colormap(flipud(cmap_temp))
    c.Label.String = 'Days after Tohoku event';
    caxis([0 800])

    % Cuboid to model seismogenic zone
    i=numel(timestep)-1;
    h2=pltCuboid(lengthAxis(j,i,3),lengthAxis(j,i,2),lengthAxis(j,i,1),G(j,i,:),...
        rad2deg(azim_cart(j,i,3)),rad2deg(elv_cart(1,i,3)));hold on;
    
    box on;grid on;axis equal;
    set(gca,'ZDir','reverse');
    xlabel('Easting (km)');
    ylabel('Northing (km)');
    zlabel('Depth (km)');
    zlim([5 12]);
    
    % Plot time-diffusivity     
    subplot(3,3,3)
    for i=1:numel(timestep)-1
        plot(tickday(i),BestD_time(j,i),'.','MarkerSize',20,'Color',cmap_temp(a-1,:));hold on;
    end
    xlabel('Days after Tohoku event');
    ylabel('Diffusivity (m^2/sec)');

    % Plot time-permeability
    subplot(3,3,6)
    k_time(j,:)=diff2perm(BestD_time(j,:));
    plot(tickday,k_time(j,:),'.','MarkerSize',20,'Color',cmap_temp(a-1,:));
    xlabel('Days from 1st event');
    ylabel('Permeability (m^2)');

    % Plot poles to PCA axises     
    subplot(3,3,7)
    MatlabPoleplot2(0,0,1);hold on;
    for i=1:numel(timestep)-1
        if lengthAxis(j,i,:)==[0 0 0]
            continue
        end
        if elv_cart(j,i,3)<0
            theta  = azim_cart(j,i,3)+pi;
            rho    = sqrt(2)*sin( (pi/2-abs(elv_cart(j,i,3)))/2);
        else
            theta  = azim_cart(j,i,3);
            rho    = sqrt(2)*sin( (pi/2-abs(elv_cart(j,i,3)))/2);
        end
        xp     = rho .* cos(theta);
        yp     = rho .* sin(theta);
        plot(xp,yp,'ok','MarkerSize',10,'MarkerFaceColor',cmap_temp(a-i,:));hold on;
        title('Major axis');
    end
    
    subplot(3,3,8)
    MatlabPoleplot2(0,0,1);hold on;
    for i=1:numel(timestep)-1
        if lengthAxis(j,i,:)==[0 0 0]
           continue
        end
        if elv_cart(j,i,2)<0
            theta  = azim_cart(j,i,2)+pi;
            rho    = sqrt(2)*sin( (pi/2-abs(elv_cart(j,i,2)))/2);
        else
            theta  = azim_cart(j,i,2);
            rho    = sqrt(2)*sin( (pi/2-abs(elv_cart(j,i,2)))/2);
        end
        xp     = rho .* cos(theta);
        yp     = rho .* sin(theta);
        plot(xp,yp,'ok','MarkerSize',10,'MarkerFaceColor',cmap_temp(a-i,:));hold on;
        title('Intermediate axis');
    end

    subplot(3,3,9)
    MatlabPoleplot2(0,0,1);hold on;
    for i=1:numel(timestep)-1
        if lengthAxis(j,i,:)==[0 0 0]
            continue
        end
        if elv_cart(j,i,1)<0
            theta  = azim_cart(j,i,1)+pi;
            rho    = sqrt(2)*sin( (pi/2-abs(elv_cart(j,i,1)))/2);
        else
            theta  = azim_cart(j,i,1);
            rho    = sqrt(2)*sin( (pi/2-abs(elv_cart(j,i,1)))/2);
        end
        xp     = rho .* cos(theta);
        yp     = rho .* sin(theta);
        plot(xp,yp,'ok','MarkerSize',10,'MarkerFaceColor',cmap_temp(a-i,:));hold on;
        title('Minor axis');
    end
    set(gcf,'Position',[159   0   860*0.9   771*0.9]); % By Hotta Dec.18
end


% figure
cmap_temp=bluesmap(:,:);
temp_Elpsday(1)=0;
temp_Elpsday(2:numel(tickday)+1)=tickday;


%% ---[Inversion with Darcy law ]---
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
    if dP_time(i,1)==targetPP
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
    if dP_time(i,1)==targetPP
        plot(tickday,dP_time(i,:),'Color',clr_prl(temp,:),'LineWidth',3);hold on;
    end
end
xlim([0 800]);
xlabel('Days after Tohoku event');
ylabel('Pressure (MPa)');
title('Pressure difference (MPa)')


% ---[Flow zone definition]------------------------------------------------
A(1,:)=lengthAxis(1,:,3)*2.*lengthAxis(1,:,1)*2;  % km2
L(1,:)=lengthAxis(1,:,2)*2; % km

A(2,:)=lengthAxis(2,:,3)*2.*lengthAxis(2,:,1)*2;  % km2
L(2,:)=lengthAxis(2,:,2)*2; % km

A(3,:)=lengthAxis(3,:,2)*2.*lengthAxis(3,:,3)*2;  % km2
L(3,:)=lengthAxis(3,:,1)*2; % km

A(4,:)=lengthAxis(4,:,2)*2.*lengthAxis(4,:,1)*2;  % km2
L(4,:)=lengthAxis(4,:,3)*2; % km


% ---[Flow rate and Volume calclation]-------------------------------------
for j=1:4 %Cluster
    for i=1:numel(tickday)
        for k=1:numel(Pout)
            Q(j,i,k)=CalcQ_Darcy(k_time(j,i),A(j,i)*10^6,L(j,i)*10^3,dP_time(k,i)*10^6,2);
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


% --- Plot Fluid volume inversion result in time series -------------------
figure
a=40;
k=21; % --> dP=170 MPa
subplot(2,1,1)
for j=1:4 %Cluster
    if j==1
        cmap_temp=bluesmap(:,:);
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
xlim([0 max(tickday)]);
xlabel('Days after Tohoku event');
ylabel('Flow rate (m^3/min)');
%     set(gca,'YScale','log');
legend('Cluster E','Cluster N','Cluster S','Cluster W');
title(['Differential pressure :',num2str(dP_time(k,1)),' MPa']);
set(gca,'YScale','log');


subplot(2,1,2)
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
xlim([0 max(tickday)]);
%     xlabel('Days from 1st event');
xlabel('Days after Tohoku event');
ylabel('Volume (m^3)');
%     set(gca,'YScale','log');
title(['Pore pressure difference:',num2str(dP_time(k,1)),' MPa']);
set(gcf,'Position',[159   0   860*0.9   771*0.9]); % By Hotta Dec.18



%% ---[Data Output]---
flag=0;
% flag=1 --> .mat format, flag=2 --> .dat (ASCII) format, flag=3 --> .nc format,  

if flag==1
    save('FluidInv_darcy_v5.mat','tickday','AllV','VCL1','VCL2','VCL3','VCL4');

elseif flag==2

elseif flag==3
    nccreate('FluidInvResult_Darcy.nc','tickday','Dimensions',{'Time',80});
    ncwrite('FluidInvResult_Darcy.nc','tickday',tickday(1:80));

    nccreate('FluidInvResult_Darcy.nc','BestDiffusivity(m2s-1)','Dimensions',{'cluster',4,'time',80});
    ncwrite('FluidInvResult_Darcy.nc','BestDiffusivity(m2s-1)',BestD_time(1:4,2:81));

    nccreate('FluidInvResult_Darcy.nc','Permeability(m2)','Dimensions',{'cluster',4,'time',80});
    ncwrite('FluidInvResult_Darcy.nc','Permeability(m2)',k_time(1:4,2:81));

    nccreate('FluidInvResult_Darcy.nc','LengthAxis(m)','Dimensions',{'cluster',4,'time',80,'Axis',3});
    ncwrite('FluidInvResult_Darcy.nc','LengthAxis(m)',lengthAxis(1:4,2:81,1:3));

    nccreate('FluidInvResult_Darcy.nc','Area of flow region(m2)','Dimensions',{'cluster',4,'time',80});
    ncwrite('FluidInvResult_Darcy.nc','Area of flow region(m2)',A(1:4,2:81));

    nccreate('FluidInvResult_Darcy.nc','Length of flow region(m2)','Dimensions',{'cluster',4,'time',80});
    ncwrite('FluidInvResult_Darcy.nc','Length of flow region(m2)',L(1:4,2:81));

    nccreate('FluidInvResult_Darcy.nc','dP(MPa)','Dimensions',{'dP',27,'time',80});
    ncwrite('FluidInvResult_Darcy.nc','dP(MPa)',dP_time(1:27,2:81));

    nccreate('FluidInvResult_Darcy.nc','Q(m3s-1)','Dimensions',{'cluster',4,'time',80,'dP',27,});
    ncwrite('FluidInvResult_Darcy.nc','Q(m3s-1)',Q(1:4,2:80,1:27));

    nccreate('FluidInvResult_Darcy.nc','AllV(m3)','Dimensions',{'time',80,'dP',27,});
    ncwrite('FluidInvResult_Darcy.nc','AllV(m3)',AllV(2:81,1:27));

    nccreate('FluidInvResult_Darcy.nc','VCL1(m3)','Dimensions',{'time',80,'dP',27,});
    ncwrite('FluidInvResult_Darcy.nc','VCL1(m3)',VCL1(2:81,1:27));

    nccreate('FluidInvResult_Darcy.nc','VCL2(m3)','Dimensions',{'time',80,'dP',27,});
    ncwrite('FluidInvResult_Darcy.nc','VCL2(m3)',VCL2(2:81,1:27));

    nccreate('FluidInvResult_Darcy.nc','VCL3(m3)','Dimensions',{'time',80,'dP',27,});
    ncwrite('FluidInvResult_Darcy.nc','VCL3(m3)',VCL3(2:81,1:27));

    nccreate('FluidInvResult_Darcy.nc','VCL4(m3)','Dimensions',{'time',80,'dP',27,});
    ncwrite('FluidInvResult_Darcy.nc','VCL4(m3)',VCL4(2:81,1:27));
end

return

%% Local Functions
% Function Darcy's law
function Q=CalcQ_Darcy(k,A,L,dP,flag)
    if flag==1
        myu= 105 * 10^(-6); %300 degC  8 km depth (80 MPa) [Pa*sec]
    elseif flag==2
        myu= 109 * 10^(-6); %300 degC 10 km depth (100 MPa) [Pa*sec]
    elseif flag==3
        myu= 91 * 10^(-6); %350 degC 8 km depth (80 MPa) [Pa*sec]
    elseif flag==4
        myu= 95 * 10^(-6); %350 degC 10 km depth (100 MPa) [Pa*sec]
        
    end
    Q=k*A/myu/L*dP; % [m3/sec]
    if L==0
        Q=0;
    end
end


%%------------------------------------------- 
function  [X,Y]=MatlabPoleplot2(azi,inc,flag)
% figure
    % input azi (Azimuth of FPS, map cordinate, clockwise from North)
    % input inc (Inclination of FPS, map cordinate, clockwise from horizontal)

    % % ---Convert Pole of Azim and inc 
    azi=azi+90;
    % azi=azi-90;
    inc=90-inc;


    % ---Convert xyz cordinate-----
    azi=90-azi;
    circle=[cos((0:360)*pi/180) ; sin((0:360)*pi/180) ; (0:360)-(0:360)]';

    if(flag==1)
      plot(circle(:,1),circle(:,2),'-k','Linewidth',1);hold on;
      xlim([-1.1,1.1]);
      ylim([-1.1,1.1]);
      zlim([-1.1,1.1]);
      view(0,90);
      axis off;
      axis equal;
      text(0,1.1,0,'N','HorizontalAlignment','center','FontSize',15);
      text(1.1,0,0,'E','HorizontalAlignment','center','FontSize',15);
      text(0,-1.1,0,'S','HorizontalAlignment','center','FontSize',15);
      text(-1.1,0,0,'W','HorizontalAlignment','center','FontSize',15);
    else
        if azi==90||azi==270
           X=0;
           Y=sin(azi.*pi./180).*sqrt(2)*sin((90-inc).*pi./(180*2));
        elseif azi==180||azi==360
           X=cos(azi.*pi./180).*sqrt(2)*sin((90-inc).*pi./(180*2));
           Y=0;
        else
           X=cos(azi.*pi./180).*sqrt(2)*sin((90-inc).*pi./(180*2));
           Y=sin(azi.*pi./180).*sqrt(2)*sin((90-inc).*pi./(180*2));
        end
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
