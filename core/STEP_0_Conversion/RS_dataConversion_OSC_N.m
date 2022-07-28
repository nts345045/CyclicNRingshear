%% Data for Oscillating load experiment 
% Author: Dougal Hansen (3/30/2022)
% Edited: Nathan Stevens (4/12/2022)
% Paths formatted for Linux/Unix operating systems

clear; close all;
load('../../../data/raw_measures/RS_RAWdata_OSC_N.mat')


%% variable  explanations

% raw_chamber_postExp = raw data for the weight of sample chamber + ice +
% bed at conclusion of experiment [mV]

% raw_postORING_torque = raw torque sensor measurement following drag test
% (chamber bolted to frame; no applied load) [mV]

% rawOringTorque = raw data for Oring drag, spinning platen rubbing against
% stationary, empty sample chamber walls [mV]

% rawPw1_zero = Best estimate for a water pressure zero for gauge 1 [V]

% rawPw2_zero = Best estimate for a water pressure zero for gauge 2 [V]


% rawTorque_zero = Torque reading when suspended in air on rack (no normal
% load) [mV]


%rawTime = time vector for experiment
%rawPw1 = raw experimental time series for water pressure from gauge 1 [V]
%rawPw2 = raw experimental time series for water pressure from gauge 2 [V]
%rawTorque = raw experimental time series for torque sensor [mV]
%rawRam = raw experimental time series for ram pressure transducer [mV]


%% Calibration coefficients

C_torque = -14.34; % Nm/mV
C_ram = 2000; % PSIG/mV
C_Pw = 50; % PSIG/V


%% Subtracted values

% zero torque (with ice-filled sample chamber installed in frame)
zero_torque = mean(rawTorque_zero); % [mV]

% mean raw transducer values
zero_Pw1 = mean(rawPw1_zero); % [V]
zero_Pw2 = mean(rawPw2_zero); % [V]


% raw sample chamber weight

raw_chamber_preEXP = 0.109; % single value recorded in lab notebook, no data file [mV]
raw_chamber_postExp = mean(raw_chamber_postExp); % mean value from data file [mV]
chamber_weight =  (raw_chamber_preEXP + raw_chamber_postExp)/2; % average of pre-and post-exp measurements [mV]

% raw o-ring drag
Oring = (mean(rawOringTorque) - mean(raw_postORING_torque)); % [mV]

%% converting time vector to UTC 

timedata=rawTime/86400;
timedata=timedata+datenum(1904,1,1); 
TimeUTC = datetime(timedata,'ConvertFrom','datenum');


%% Sample chamber dimensions

rout = .3; % outer radius [m]
rin =.1; % inner radius [m];
a_ring = pi*(.30^2 - .10^2); % area of the ring [m^2]


%% Converting raw Ram pressure data to kPa

% Pressure on the ram [PSIG]
Ram_pressure = C_ram * rawRam; % Pressure on the ram [PSIG]

% Convert pressure on ram to axial load in [kN] (Peter's calibration)
% Load = (Ram_pressure - chamber_weight)* 0.0231252; % Convert pressure on ram to axial load in [kN] (Peter's calibration)
%% NTS Edit based on DDH comment on April 10. 2022
Load = (Ram_pressure - C_ram*chamber_weight)* 0.0231252;

% Applied normal stress in [kPa]
SigmaN =  Load/a_ring; 


%% Converting water pressure data to kPa

psi2kpa = 6.895; % conversion factor for PSIG to kPa

%Corrected gauge records
Pw1 = C_Pw * psi2kpa * (rawPw1 - zero_Pw1); % [kPa]
Pw2 = C_Pw * psi2kpa * (rawPw2 - zero_Pw2); % [kPa]


%% Converting torque record to kPa

% Applied Torque (includes Peter's correction for sensor dependence on axial load)
Torque = (C_torque * (rawTorque - Oring - zero_torque)) - (0.243*Ram_pressure); % [Nm]

%Converting Torque to shear stress
Tau = (3/(2*pi) * (Torque) / (rout^3 - rin^3) /1000 ) ; %[kPA]


%% PLOTTING Normal stress, torque, and shear/stress torqe
figure; 

subplot(3, 1, 1);
plot(TimeUTC, SigmaN, 'k.')
xlim([TimeUTC(1) TimeUTC(end)])
ylabel('\sigma_N (kPa)')
set(gca, 'FontSize', 12)


subplot(3, 1, 2);
plot(TimeUTC , Tau, 'r.');
xlim([TimeUTC(1) TimeUTC(end)])
ylabel('\tau_b (kPa)')
set(gca, 'FontSize', 12)


subplot(3, 1, 3);
plot(TimeUTC, Pw1, 'c.'); 
hold on 
plot(TimeUTC, Pw2, 'b.'); 
legend('P_{w}1', 'P_{w}2');
ylabel('P_w (kPa)');
xlabel('Time (UTC)')

% Save data as output *.mat file
save('../../../data/physical_units/Pressures_Stresses/RS_procData_OSC_N.mat', 'TimeUTC', 'Tau', 'SigmaN', 'Pw1', 'Pw2')
% Save tables as *.csv files


% Compose tables
N_kPa = SigmaN;
T_kPa = Tau;
Pw1_kPa = Pw1;
Pw2_kPa = Pw2;
% Do UTC Time indexed data 
t_NTP_UTC = table(TimeUTC,N_kPa,T_kPa,Pw1_kPa,Pw2_kPa);

CSV_OUT = '../../../data/physical_units/Pressures_Stresses/UTC_N_T_P.csv'


%% Write tables to CSV
writetable(t_NTP_UTC,CSV_OUT);

