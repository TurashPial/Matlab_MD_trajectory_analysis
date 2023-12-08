%------------pure-pressure - bulk------------
% this code takes input from MD ion distribution and that info to employ it in a FDM to generate continuum flow distribution 
clear all
load('monomer_volume_fraction_gd05_f1_70_82_ns')  % gd0.05 dp1
% load('monomer_volume_fraction_gd05_f2_55_67_ns')  % gd0.05 dp2
% load('monomer_volume_fraction_gd03_f1_17_29_ns')  % gd0.03 dp1
nsteps = 3000;

x = vol_frac';
y = block_cen';
X = vertcat(0,x,0);

% % Y = vertcat(0,y,120.4509);
% Y = vertcat(62.6215,y,180.94);
Y = vertcat(46.5341,y,166.9850);
Y1 = Y - 46.5341;   % gd 0.05
% Y1 = Y - 62.6215; %gd 0.03
N = nsteps+1;
yq = linspace(0,120.4509,nsteps+1)';  %gd 0.05
% yq = linspace(0,180.94-62.6215,nsteps+1)'; % gd 0.03
M = interp1(Y1,X,yq,'linear','extrap');
h = max(Y1)/2;
ak=1.53*1e-10;
dy=(yq(2)-yq(1))*1e-10;

%%


% bulk-region
M1 = round((89.1938-46.5341)*(1e-10)/dy); %gd 0.05 70_82ns
M2 = round((119.3065-46.5341)*(1e-10)/dy); %gd 0.05 70_82ns
% % % 
% M1 = round(42.6597*(1e-10)/dy); %gd 0.05 55_67_ns
% M2 = round((72.7724)*(1e-10)/dy); %gd 0.05 55_67_ns
% 
% M1 = round(41.9045*(1e-10)/dy); %gd 0.03 17_29ns
% M2 = round((76.414)*(1e-10)/dy); %gd 0.03 17_29ns

% inter-region
% M1 = round((76.3108-46.5341)*(1e-10)/dy); %gd 0.05 70_82ns
% M2 = round((134.9742-46.5341)*(1e-10)/dy); %gd 0.05 70_82ns
% % 
% M1 = round(27.944*(1e-10)/dy); %gd 0.05 55_67_ns
% M2 = round((2*h-32.162)*(1e-10)/dy); %gd 0.05 55_67_ns

% M1 = round(26.854*(1e-10)/dy); %gd 0.03 17_29ns
% M2 = round((2*h-26.206)*(1e-10)/dy); %gd 0.03 17_29ns


% Constants
ere0=8.8*79.8*10^(-12); % permittivity 
NA=6.023*10^23;         % Avogadro constant
kB=8.314/NA;            % Boltzmann constant J/K
T=300;%298;                  % temperature, k
e=1.6*10^(-19);         % electric charge, c 
eta=8.90*10^(-4);       % viscosity, Pa*s

A = zeros(M2-M1+2,M2-M1+2);    
A(1,1) = 1;
A(M2-M1+2,M2-M1+2) = 1;

for i = 2:(M2-M1+1)
    A(i,i)=-2/dy^2;%- M(i+M1-2)^2/ak^2;
    A(i,i-1)=1/dy^2;
    A(i,i+1)=1/dy^2;
end

dp_dx = -1e15;
% dp_dx = -2e15;
B = linspace(0,0,M2-M1+2)';
B(2:M2-M1+1) = dp_dx/eta;

% % % % B(1) = 1.54164; B(M2-M1+2)=1.92086; % gd 0.05 f1_70_82_ns

B(1) = 1.54535; B(M2-M1+2)=1.94432; % gd 0.05 f1_70_82_ns
% B(1) = 3.60915; B(M2-M1+2)=3.81705; % gd 0.05 f2_55_67_ns
% B(1) = 2.07213; B(M2-M1+2)=2.50219; % gd 0.03 f1_17_29_ns
u = A\B;
plot(u,yq(M1:M2+1))
Q_in=trapz(yq(M1:M2+1),u)*1e-10*93.9e-10;
%  Q_in=trapz(yq(M1:M2+1),u)*1e-10*101.007e-10;
%%
%%%%%%%%%%%%%%%%%%%%Electro-osmotic bulk%%%%%%%%%%%%%%%%%%%%%%%%%%
%  clear all
% load('monomer_volume_fraction_gd05_f1_70_82_ns')  % gd0.05 dp1
% load('n_Na_Cl_gd05_f1_70_82_ns')


% load('monomer_volume_fraction_gd05_f2_55_67_ns')  % gd0.05 dp2
% load('n_Na_Cl_gd05_f2_55_67_ns')

load('monomer_volume_fraction_gd03_f1_17_29_ns')  % gd0.03 dp1
load('n_Na_Cl_gd03_f1_17_29_ns')


x = vol_frac';
y = block_cen';
n_Na_q = vertcat(0,n_Na',0);
% y_Na_q = vertcat(0,y_Na',120.4509);
y_Na_q = vertcat(0,y_Na',180.94-62.6215);
n_Cl_q = vertcat(0,n_Cl',0);
% y_Cl_q = vertcat(0,y_Cl',120.4509);
y_Cl_q = vertcat(0,y_Na',180.94-62.6215);
X = vertcat(0,x,0);
% % Y = vertcat(0,y,120.4509);
Y = vertcat(62.6215,y,180.94);
% Y = vertcat(46.5341,y,166.9850);
% Y1 = Y - 46.5341;   % gd 0.05
Y1 = Y - 62.6215; %gd 0.03
nsteps=3000;
N = nsteps+1;
% yq = linspace(0,120.4509,nsteps+1)';  %gd 0.05
yq = linspace(0,180.94-62.6215,nsteps+1)'; % gd 0.03
M = interp1(Y1,X,yq,'linear','extrap');
N_Na_q = interp1(y_Na_q,n_Na_q,yq,'linear','extrap');
N_Cl_q = interp1(y_Cl_q,n_Cl_q,yq,'linear','extrap');
h = max(Y1)/2;
ak=1.53*1e-10;
dy=(yq(2)-yq(1))*1e-10;


% bulk-region
% % % % M1 = round(42.6597*(1e-10)/dy); %gd 0.05 70_82ns
% % % % M2 = round((72.7724)*(1e-10)/dy); %gd 0.05 70_82ns
% % 

% M1 = round((89.1938-46.5341)*(1e-10)/dy); %gd 0.05 76_86ns
% M2 = round((119.3065-46.5341)*(1e-10)/dy); %gd 0.05 76_86ns

% M1 = round(42.6597*(1e-10)/dy); %gd 0.05 55_67_ns
% M2 = round((72.7724)*(1e-10)/dy); %gd 0.05 55_67_ns

M1 = round(41.9045*(1e-10)/dy); %gd 0.03 17_29ns
M2 = round((76.414)*(1e-10)/dy); %gd 0.03 17_29ns

% inter-region
% % % % M1 = round(29.883*(1e-10)/dy); %gd 0.05 70_82ns
% % % % M2 = round((2*h-32.045)*(1e-10)/dy); %gd 0.05 70_82ns

% M1 = round((76.3108-46.5341)*(1e-10)/dy); %gd 0.05 70_82ns
% M2 = round((134.9742-46.5341)*(1e-10)/dy); %gd 0.05 70_82ns
% 
% M1 = round(27.944*(1e-10)/dy); %gd 0.05 55_67_ns
% M2 = round((2*h-32.162)*(1e-10)/dy); %gd 0.05 55_67_ns

% M1 = round(26.854*(1e-10)/dy); %gd 0.03 17_29ns
% M2 = round((2*h-26.206)*(1e-10)/dy); %gd 0.03 17_29ns

n_p = N_Na_q*1e30;
n_m = N_Cl_q*1e30;
% Constants
    ere0=8.8*79.8*10^(-12); % permittivity 
    NA=6.023*10^23;         % Avogadro constant
    kB=8.314/NA;            % Boltzmann constant J/K
    T=300;%298;                  % temperature, k
    e=1.6*10^(-19);         % electric charge, c 
    eta=8.90*10^(-4);       % viscosity, Pa*s

    %Elecroosmotic flow field
    
       A=zeros(M2-M1+2,M2-M1+2);
        %B.C. u=0   @ left y=-h
        A(1,1)=1;
        
        A(M2-M1+2,M2-M1+2) = 1;
        % -1<y<-1+d in PE
        for i=2:(M2-M1+1)
            A(i,i)=-2/dy^2;%- M(i+M1-2)^2/ak^2; %M(i)^2 *(h*1e-10)^2/(ak^2);
            A(i,i-1)=1/dy^2;
            A(i,i+1)=1/dy^2;
        end
        
%         
%         E = -2.38*1e-3*1e9; % Case 1 - 70_82
%         E = -3.68*1e-3*1e9; % case 2 - 55_67
        E = -4.34*1e-3*1e9; % Case 3 - 17_29
        
%% [B]  
        B=linspace(0,0,M2-M1+2)';
        B(2:(M2-M1+1))= e*E/eta*(-n_p(M1+1:M2) + n_m(M1+1:M2));%1;
        
%         B(1) = 0.0147; B(M2-M1+2)=0.0147; % gd 0.05 f1_70_82_ns
%         B(1) = 0.0227; B(M2-M1+2)=0.0227; % gd 0.05 f2_55_67_ns
        B(1) = 0.0268; B(M2-M1+2)=0.0268; % gd 0.03 f1_17_29_ns
%         B(MM)=0;
        u_EOS= A\B;
        
%         dpdx = -1e15;
% %         dpdx = -2e15;
%         u_E0=((h*1e-10)^2/eta)*dpdx;
%         u_EOS_dim=u_EOS*u_E0;
        plot(u_EOS,yq(M1:M2+1))
%         Q_in=trapz(yq(M1:M2+1),u_EOS)*1e-10*93.9e-10;
         Q_in=trapz(yq(M1:M2+1),u_EOS)*1e-10*101.007e-10;
