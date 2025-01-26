clc; clear all; close all; format long;

yr  = 24;
mon = 10;
day = 24;
filenum_max = 2;

% dir = '/Users/mitchellshen/Dropbox (Personal)/ABM_JG-MMS/Hi-ABM log/';
% dir = '/Users/mitchellshen/OneDrive - Princeton University/Beam_CAL_database/ABM/';
% dir = '/Users/mitchellshen/OneDrive - Princeton University/IMAP-Lo_CAL_database/FM Cal/ABM_reference/T101_TOF_Cal/';
% dir = '/Users/mitchellshen/OneDrive - Princeton University/IMAP-Lo_CAL_database/FM Cal/ABM_reference/T102_Pre-Cal1/';
dir = '/Users/mitchellshen/OneDrive - Princeton University/IMAP-Lo_CAL_database/FM Cal/ABM_reference/T104_Pre-Cal2/';
dir = ['/Users/mitchellshen/Library/CloudStorage/OneDrive-PrincetonUniversity/' ...
       'IMAP-Lo_CAL_database/FM Cal/ABM_reference/T105_Final_Cal/'];
folder=['HVPS-ABM-Data-20' num2str(yr,'%02d') '.' num2str(mon,'%02d') '.' num2str(day,'%02d') '/'];

tbl = zeros(filenum_max,15);

for i = 1:filenum_max
data = readtable([dir folder 'ABM-Counts_20' num2str(yr,'%02d') '.' num2str(mon,'%02d') '.' ...
                  num2str(day,'%02d') '_' num2str(i) '.csv']);
t0   = datetime(table2array(  data(1,1)),'InputFormat','MM/dd/yyyy/HH:mm:ss.S');
t1   = datetime(table2array(data(end,1)),'InputFormat','MM/dd/yyyy/HH:mm:ss.S');
ttol = seconds(t1-t0);
corr = size(data,1)./ttol;

tbl(i, 1) = i;
tbl(i, 2) = ttol;
tbl(i, 3) = mean( table2array(data(:,5)) ) * corr; % FCEM_ave 
tbl(i, 4) = mean( table2array(data(:,4)) ) * corr; % ACEM_ave 
tbl(i, 5) = mean( table2array(data(:,6)) ) * corr; % COIN_ave
tbl(i, 6) = std(  table2array(data(:,5)) ) * corr ./ sqrt( size(data,1) ); % FCEM uncertainty
tbl(i, 7) = std(  table2array(data(:,4)) ) * corr ./ sqrt( size(data,1) ); % ACEM uncertainty
tbl(i, 8) = std(  table2array(data(:,6)) ) * corr ./ sqrt( size(data,1) ); % COIN uncertainty

% absolute detection efficiency
ADE_N_bin = table2array(data(:,6)).^2 ./ ( table2array(data(:,5)).*table2array(data(:,4)) ); 
ADE_Y_bin = sum( table2array(data(:,6)) )^2 ./ ... 
          ( sum( table2array(data(:,5)) )*sum( table2array(data(:,4)) ) ); 
cts_N_bin = table2array(data(:,5))    .*   table2array(data(:,4))./table2array(data(:,6)  ) * corr;
cts_Y_bin =(sum( table2array(data(:,5)) )*sum( table2array(data(:,4)) ) ) ./ ...
            sum( table2array(data(:,6)) ) ./ ttol;
      
tbl(i, 9) = mean(ADE_N_bin(isfinite(ADE_N_bin)));   % ADE seperately calculated with mean 
tbl(i,10) =  std(ADE_N_bin(isfinite(ADE_N_bin)))./...
            sum(isfinite(ADE_N_bin));               % ADE uncertainty, biased when more zero counts

tbl(i,11) = ADE_Y_bin;                              % ADE binned, single value
tbl(i,12) =  std(ADE_N_bin(isfinite(ADE_N_bin)))./...
            sum(isfinite(ADE_N_bin));               % ADE uncertainty

tbl(i,13) = mean(cts_N_bin(isfinite(cts_N_bin)));   % cts seperately calculated with mean of non-zeros
tbl(i,14) =  std(cts_N_bin(isfinite(cts_N_bin)))./...
            sum(isfinite(cts_N_bin));               % cts uncertainty, biased when more zero counts

tbl(i,15) = cts_Y_bin;                              % cts binned, single value
tbl(i,16) =  std(cts_N_bin(isfinite(cts_N_bin)))./...
            sum(isfinite(cts_N_bin));               % cts uncertainty

end

writematrix(tbl,[dir 'tbl_' folder(15:end-1) '.txt'],'Delimiter','tab');
disp('Process done');
