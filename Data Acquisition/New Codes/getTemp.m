%% session vars

% daq vars
fscan = 200; %scan rate (sampling , Hz)
Vr = 1; %input voltage range (+/- volts); =1 gives -100 to +100 C

% temp sensor vars
sf = 0.010; %volts per centrigrade scaling factor
Vs = 5.0; %sensor source voltage
nsensor = 4; %num of sensors

%% data: num of samples = length of session / length of sample

global temp;
global cond;

%% daq connection and setup
global s; s = daq.createSession('ni');
s.Rate = fscan;
s.IsContinuous = true;

%% input: use four inputs in differential mode
s.addAnalogInputChannel('Dev1',0:3,'Voltage'); %add channels ai0--3
for j = 1:nsensor
    s.Channels(j).Range=[-Vr Vr]; %set input range, all channels
end

%% output: used to power sensors (generate dc voltage of 'Vs')

% Start Measurements (in Loop: Its Beginning)
lh_s = s.addlistener('DataAvailable',@UpdateStatus);
temp = zeros(1,nsensor);
cond = 0;
s.startBackground();

% Stop Measurements (in Loop: Its End)
% cond = 1;
% clear lh_s;

% Releasing the Hardware Resources (Outside of Loop)
% s.release();