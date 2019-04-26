%% ********************** Important Variables **********************

clear all;
close all;
% save directory (remember trailing slash)
recDir = 'C:\Users\shayan\Desktop\data_collect\b5c1\';
% number of records to capture
n = 10816; % record every 7 secs for 24 hrs (two min of safeguard)

% *******************************************************************

%% ******************* Supply Voltage Sensor Section *******************

% Creating Analog Input Object for MCC Device
ai = analoginput('mcc');

% Adding Channels for Measuring Supply Voltage
addchannel(ai,0:1);

% Specifying Sample Rate (Hz)
ai.SampleRate = 500;

% Making Data Acqusition Continuous
ai.SamplesPerTrigger = Inf;

% The Sensors are Triggered Immediately by Start Function
ai.TriggerType = 'Immediate';

% Creating a Log File for Scans
ai.LogToDiskMode = 'Overwrite';
ai.LogFileName = 'volt.daq';
ai.LoggingMode = 'disk';

% *******************************************************************

%% ******************* Temperature Sensor Section *******************

% Session Variables:
% 1) DAQ Variables
fscan = 500; % Scan Rate (Sampling , Hz)
Vr = 1; % Input Voltage Range , (+/- volts)=1 gives -100 to +100 C
% 2) Temperature Sensor Variables
sf = 0.010; % Volts per Centrigrade Scaling Factor
Vs = 5.0; % Sensor Source Voltage
nsensor = 4; % Number of Sensors

% Average of Collected Data
global temp;

% DAQ Connection and Setup
global ni_temp; ni_temp = daq.createSession('ni'); % Creating Session
ni_temp.Rate = fscan; % Number of Scans per Second
ni_temp.IsContinuous = true;
ni_temp.NotifyWhenDataAvailableExceeds = 50;

% Configuring Input Channels
ni_temp.addAnalogInputChannel('Dev1',0:3,'Voltage'); %add channels ai0--3
for j = 1:nsensor
    ni_temp.Channels(j).Range=[-Vr Vr]; %set input range, all channels
end

% *******************************************************************

%% ********************** Oscilloscope Section **********************

% length of record
l = 1e6; % (??? or 500e3)

% CONNECTION VAR

% Find a VISA-GPIB object.
SCP = instrfind('Type', 'visa-gpib', 'RsrcName', 'GPIB8::1::0::INSTR', 'Tag', '');

% Create the VISA-GPIB object if it does not exist
% otherwise use the object that was found.

if isempty(SCP)
    SCP = visa('TEK', 'GPIB8::1::0::INSTR');
else
    fclose(SCP);
    SCP = SCP(1);
end

% Configure instrument object, SCP

set(SCP, 'InputBufferSize', l);

% Configure instrument object, SCP

% set(SCP, 'OutputBufferSize', 512);

if isempty(SCP)
    disp('Scope not present!');
    return;
end

% Connect to instrument object, SCP.
fopen(SCP);

% Configure binary field data width for the waveform
fprintf(SCP,'WFMO:BYT_N 1');

% determine y-increment (voltage scale)
fprintf(SCP,'DATa:SOUrce CH1');
yinc = str2num(query(SCP,'WFMO:YMU?'));

fprintf(SCP,'DATa:SOUrce CH2');
if (yinc ~= str2num(query(SCP,'WFMO:YMU?')))
    disp('Vertical increments are not matched!');
    fclose(SCP);
    return;
end

% sample rate
s = str2num(query(SCP,'HORizontal:MODE:SAMPLERate?'));

% set data collection points
fprintf(SCP,'DATa:STARt 1'); %beginning of record, in buffer of 1*10^6
fprintf(SCP,'DATa:STOP 1e6'); %end of record, in buffer % (??? or 500e6)

i = 0; %current rec cnt
ib = 0; %current bad rec cnt

% begin timer
tic;
% instruct scope to take single measurement
fprintf(SCP,'ACQuire:STATE ON');

% *******************************************************************

%% *********************** Associated Section ***********************

disp('Run initiated...');
while (i < n)
    
    % * Start Supply Voltage Measurements*
    start(ai);
    
    % * Start Temperature Measurements *
    lh_s = ni_temp.addlistener('DataAvailable',@UpdateStatus);
    temp = [];
    ni_temp.startBackground();

    % wait for scope to finish measurement (check scope state every 1/10
    % of a second)
    disp('Waiting for trigger...')
    while str2num(query(SCP,'ACQuire:STATE?')) == 1
            pause(0.1);
    end
    disp('...triggered!')
    
    % get channel data
    fprintf(SCP,'DATa:SOUrce CH1');
    fprintf(SCP,'CURVe?');
    ch1 = binblockread(SCP,'int8')';
    
    fprintf(SCP,'DAT:SOURCE CH2');
    fprintf(SCP,'CURVe?');
    ch2 = binblockread(SCP,'int8')';
    
    % to save time, trigger new measurement on scope
    fprintf(SCP,'ACQuire:STATE ON');
    
    % * Stop Supply Voltage Measurements *
    stop(ai);
    temp_volt = daqread('volt.daq');
    volt_dat = mean(temp_volt);

    % * Stop Temperature Measurements *
	stop(ni_temp);
	clear lh_s;

    % recover differential signal
    rec = yinc*(ch2-ch1);
    
    ts = now;

    % check if record is good
    if (max(rec(50e3:60e3)) > 3 ... % sync signal
            && min(rec(150e3:160e3)) < -2.5 ... % payload
            && max(rec(990e3:1000e3)) < 1) % end of record

        % increment good recs counter
        i = i + 1;

        is = num2str(i);
        
        % * Storing Desired Variables *
        f = [recDir 'sample' ...
            strrep(num2str(zeros(1,5-length(is))), ' ', '') is '.mat'];
        save(f,'ts','s','yinc','rec','ch1','ch2','temp','volt_dat');
    else
        % increment bad recs counter
        ib = ib + 1;
        
        ibs = num2str(ib);
        
        % * Storing Desired Variables *
        f = [recDir 'bsample' ...
            strrep(num2str(zeros(1,5-length(ibs))), ' ', '') ibs '.mat'];
        save(f,'ts','s','yinc','rec','ch1','ch2','temp','volt_dat');
    end
    
    t = round(toc); % get current running time
    disp(['Record: ' num2str(i) '/' num2str(n) ...
        ' (%' num2str((i)/n*100) '); ' ...
        'Bad record count: ' num2str(ib) ';' ...
        ' Elapsed time: ' num2str(floor(t/3600)) ' hr. ' ...
        num2str(mod(floor(t/60),60)) ' min. ' ...
        num2str(mod(t,60)) ' sec.']);
    
    % Showing the Measurements to User
    plot(rec);title(num2str([temp,volt_dat]));

end

disp('...run complete');

% Releasing Supply Voltage Sensors
delete(ai);

% Releasing Temperature Sensors
ni_temp.release();

% close connection to scope
fclose(SCP);

% *******************************************************************
