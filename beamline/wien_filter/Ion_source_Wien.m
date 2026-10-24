close all; clear all; clc;
format long

q = 1.6E-19;              % C
mp= 1.67262192369E-27;    % kg
M_amu = 1.66E-27;
d = 0.0090 %86%875;      % m
B = 0.04;       % T
E = 90;    % kV

Ma_Xe = 131.290;
Ma_Ar = 39.9480;    % amu, Ar+
Ma_N2 = 28.0134;    % amu, N2+
Ma_Ne = 20.1797;    % amu, Ne+
Ma_N  = 14.0067;    % amu, N+
Ma_He =  4.0026;    % amu, He+
Ma_DH =  3;
Ma_H2 =  2.0156;    % amu, H2+
Ma_H  =  1.0073;    % amu, H

U = (d * sqrt((2*B^2*E*q) ./ ...
    ([Ma_Xe Ma_Ar Ma_N2 Ma_Ne Ma_N Ma_He Ma_DH Ma_H2 Ma_H]*M_amu)))'

% d = [188 270.4 389.5 480.2 559.6 629.6]./sqrt(2*0.04^2*[0.5 1 2 3 4 5]*1e3*q./(2*1.66E-27)) % H2, d = 0.022
% d = [134.1 193.4 281.5 345 407 454.7  ]./sqrt(2*0.04^2*[0.5 1 2 3 4 5]*1e3*q./(4*1.66E-27)) % He, d = 0.022
% d = [111 157 192]./sqrt(2*0.04^2*[1 2 3]*1e3*q./(2*1.66E-27)) % H2


%% PLOT WIEN FILTER PLATE POTENTIAL SCAN
format long

path = '/Users/mitchellshen/Downloads/';
num  = 2;
yyyy = 2023;
mm   = 02;
dd   = 14;
name = [path 'Wien-Filter-Plate-Potential-Scan__' ...
        num2str(yyyy) '.' num2str(mm,'%02d') '.' num2str(dd,'%02d') '_' num2str(num) ...
        '.txt'];
% data = readtable(name);

% folder = 'C:\Users\ms3648\Dropbox (Princeton)\Space-Lab-Starter-Kit\Ion-Beam-Calibration\WienFilter Scan\';
% name = '/Users/mitchellshen/Downloads/Wien-Filter-Plate-Potential-Scan__2021.10.29_2.txt';
% name = [folder 'Wien-Filter-Plate-Potential-Scan__2021.10.27_1.txt'];
data = dlmread(name);

V_wien = data(:,1);
I_pico = data(:,2);

% figure(1)
% plot(V_wien,(I_pico-min(I_pico))*1E12,'+')

% V_wien_seq = V_wien; %min(V_wien):1:max(V_wien);
V_wien_seq = [0:1:256 258:2:300]';
for i = 1:1:length(V_wien_seq)
% I_pico_seq(i) = (mean( I_pico(V_wien==V_wien_seq(i)) ) -min(I_pico) )*1E12; 
I_pico_seq(i) = mean(I_pico(V_wien==V_wien_seq(i)))*1E12; 
end

I_pico_seq = I_pico_seq';
I_pico_seq = I_pico_seq-min(I_pico_seq);

figure
set(gcf,'Position', [50, 150, 1250, 600 ])
subplot(2,1,1)
% plot(V_wien_seq,I_pico_seq(isfinite(I_pico_seq)),'-','linewidth', 1.6,'color', [0.8500 0.3250 0.0980])
plot(V_wien_seq,I_pico_seq,'-','linewidth', 1.6,'color', [0.8500 0.3250 0.0980])
hold on; grid on;
xlabel('Wien filter potential (V)')
ylabel('Beam current (pA)')
title('Ion source Wien filter scanning (Energy: 1kV)')
xticks(0:20:300)
% ylim([0 0.1])
set(gca,'lineWidth',1,'FontSize',14,'FontWeight','bold');


amu = (d./V_wien_seq).^2 * (2*B^2*E*q./M_amu);
subplot(2,1,2)
plot(amu,(I_pico_seq),'+','linewidth', 1.6,'color', [0.8500 0.3250 0.0980],'markersize',16)
hold on; grid on;
xlabel('Atomic mass (amu)')
ylabel('Beam current (pA)')
title('Ion source Wien filter scanning (Energy: 1kV)')
xticks([1 2:2:60]); xlim([1 60])
% ylim([0 0.06])
set(gca,'lineWidth',1,'FontSize',14,'FontWeight','bold');

if 0
ax2 = axes('Position',[.170 .276 .130 .120]);
box on;
plot(amu,(I_pico_seq),'+','linewidth', 1.6,'color', [0.8500 0.3250 0.0980],'markersize',10)
hold on; grid on;
xlim([0 5]);    xticks(0:1:5)
ylim([0 0.08]); yticks([0:0.02:0.08])
set(gca,'lineWidth',1,'FontSize',12,'FontWeight','bold');
end

set(0, 'DefaultFigurePaperPositionMode', 'auto')
print('-dpng'  ,'-r300', [name(1:end-4) '_v2.png'])
% print('-dtiff' ,'-r300',savename(1:end-4))
disp('Mission Completed')