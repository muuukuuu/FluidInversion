function k=diff2perm(D)

hoge=2;
if hoge ==1
    k=diff2perm_local(D);
elseif hoge==2
    k=diff2perm_shapiro(D);
elseif hoge==3
    k=diff2perm_TZ(D);
end

end

%% Functions: Permeability estimation from difusivity
% Various parameter set and various theory

% Prameter sets (KTB) Shapiro et al., 1997 GJI
function k=diff2perm_local(D)
    phi=0.023*10^(-3);
    Kf = 2.3*10^9;
    Kd = 5.0*10^10;
    Kg = 36.4*10^9;
    alpha = 1-Kd/Kg;
    if alpha<0
        alpha=0;
    end
    eta=10^(-3);
    
    N=(phi/Kf + alpha/Kg)^(-1);
%     N=(phi/Kf + (alpha-phi)/Kg)^(-1);

    k=D/N*eta;

    disp(['Apparent permeability =',num2str(k),' m2'])
end


% Prameter sets for this field:
% Equation: Method; Shapiro et al., 1997 GJI
function k=diff2perm_shapiro(D)
%     phi=0.023*10^(-3); % Porosity
    phi=0.62*10^(-2); % Porosity

    Kf = 2.692*10^9; % Pa, 2692 MPa
    % IAPWS95 formulation    
       
    Kd = 5.0*10^10; % Pa, 50GPa 
    % Iwamori, et al., 2021; Nakajima et al., 2001
    
    Kg = 4.87*10^10; % Pa, 48.7GPa
    % Granite mineral composition: quartz 20, plagioclase 54, potash feldspar 23, biotite 3 vol.%
    % 318 MPa, 375degC (12 km assumed)
    
    eta=1.24*10^(-4); %Pa s, 
    % IAPWS95 formulation
    
    % Biot system conversion Shapiro,1999
    alpha = 1-Kd/Kg;
    
    N=(phi/Kf + (alpha-phi)/Kg)^(-1);
    k=D/N*eta;
    disp(['Apparent permeability (Shapiro) =',num2str(k),' m2']);
end


% Prameter sets for this field:
% Equation: Method; Townend and Zoback, 2000 Geology
function k=diff2perm_TZ(D)
%     phi=0.023*10^(-3); % Porosity
    phi=0.62*10^(-2); % Porosity
    
    Kf = 2.692*10^9; % Pa, 2692 MPa
    % IAPWS95 formulation    
       
    Kd = 5.0*10^10; % Pa, 50GPa 
    % Iwamori, et al., 2021; Nakajima et al., 2001
    
    Kg = 4.87*10^10; % Pa, 48.7GPa
    % Granite mineral composition: quartz 20, plagioclase 54, potash feldspar 23, biotite 3 vol.%
    % 318 MPa, 375degC (12 km assumed)
    
    eta=1.24*10^(-4); %Pa s, 
    % IAPWS95 formulation
    
    Kr=Kd;
    
    N=(phi/Kf + 1/Kr)^(-1);
    k=D/N*eta;
    disp(['Apparent permeability by Townend & Zoback=',num2str(k),' m2']);
end