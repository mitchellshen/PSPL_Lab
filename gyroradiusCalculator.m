%% Simple Gyroradius calculator and simple deflection
E = 1E3; %10 MeV in eV
mp = 16*1.67262192E-27; % proton mass * 16 = ox
L = 1; %m
BE = 10*5E-5; % T

[rg_HC,rgHC_km,rgHC_au,dXHC_m,dXHC_cm] = computeRg_dX(L,BE,mp,E);

%% A loop of Gyroradius calculator
Brange = [1:1:4000].*BE;
L = 0.0508; %m
for i = 1:1:length(Brange)
    B = Brange(i);
    [rg(i,1),rg_km(i,1),rg_au(i,1),dX_m(i,1),dX_cm(i,1)] = computeRg_dX(L,B,mp,E);
end

%% Detailed deflection 
% Universal parameters
eV2J = 1.60218E-19; %1eV to J
q = 1.60218E-19; mp = 1.67262192E-27; 
R = 0.1651; 
Lfull = 3.77; 

% Compute the magnetic field per current that we imposed...
I = 1.2; % HC current (reverse magnetic field)
n = 85; mu = 1.26E-6;
Bz2 = (4/5)^(3/2)*mu*n.*I/R+5.11E-5

% User input
amu = 4; % AMU of the particle of interest. 
L = 0.89; % Magnet length

% Compute the dy assume Vx = constant and Vy0 ~0
m = mp*amu;
L1 = (Lfull - L)./2; L3 = L1; L2 = L;

vperp = sqrt(E*eV2J*2./m);
t1 = L1/vperp; t2 = L2/vperp; t3 = L3/vperp;
Bz1 = 5E-5; Bz3 = 5E-5;

dy_m = q/m*Bz1*(vperp*t1^2/2+L1*t2+L1*t3)+q/m*Bz2*(vperp*t2^2/2-L2*t3)...
    + q/m*Bz3*(t3^2/2)
dy_cm = dy_m*1E2

%% Plot
figure();
yyaxis left;
plot(Brange,rg,'k','LineWidth',2); hold on;
xlabel('B field Factor (B/B_E)'); ylabel('Gyroradius (m)');
set(gca,'FontSize',12,'YScale','log');

yyaxis right;
plot(Brange,dX_cm,'r','LineWidth',2); hold on;
plot([Brange(1) Brange(end)],[8 8],':','LineWidth',1);

xlabel('B field Factor (B/B_E)'); ylabel('deflection (cm)');
set(gca,'FontSize',12);

%% Create function
function [rg,rg_km,rg_au,dX_m,dX_cm] = computeRg_dX(L,B,mp,E)
    eV2J = 1.60218E-19; %1eV to J
    q = 1.60218E-19;
    
    vperp = sqrt(E*eV2J*2./mp);
    
    rg = mp.*vperp./(q.*B); %m
    rg_km = mp.*vperp./(q.*B)./1000; %km
    rg_au = rg.*6.68459e-9; % AU
    
    %% dX
    dX_m = rg*(1-cos(asin(L/rg)))
    dX_cm = dX_m/1E-2
end
