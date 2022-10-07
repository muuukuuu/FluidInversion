%% Function get best Diffusivity %%
% Original version is from Free University Berlin, Dr. Serge Shapiro group
% Mainly written by PhD student Carsten Dinske

function [BestD]=Find_bestD(Devents,percent)

steps=1000;
Dvec=linspace(0,max(Devents),steps);
percent1(steps)=0;
BestD=0;

for i=1:steps;
    a=find(Devents < Dvec(i) );
    percent1(i)=numel(a)/numel(Devents)*100;
    if(percent1(i)>percent);
       BestD=Dvec(i);
       disp(['Best fitting hydraulic diffusivity:  D = ',num2str(BestD),'m^2/s']);
       break;
    end; 
end