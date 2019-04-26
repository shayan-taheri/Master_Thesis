%% session vars
T = 54000;%1*60; %session length: number of seconds to take data for

% daq vars
fscan = 100; %scan rate (sampling , Hz)
tscan = 1; %scan time (sec): length of single sample (will average fscan*tscan samples)
Vr = 1; %input voltage range (+/- volts); =1 gives -100 to +100 C

% temp sensor vars
sf = 0.010; %volts per centrigrade scaling factor
Vs = 5.0; %sensor source voltage
nsensor = 4; %num of sensors

%% data: num of samples = length of session / length of sample
global temp; temp = zeros(ceil(T/tscan),nsensor); %temperature data
global time; time = zeros(ceil(T/tscan),1); %time data acquired
global i; i=1; %sample number

%% daq connection and setup
s = daq.createSession('ni');
s.Rate = fscan;
s.DurationInSeconds = T; %session time should be rounded up so we always get same number of points per sample
s.NotifyWhenDataAvailableExceeds = fscan*tscan; %listener called after this many samples acquired
% sOut = daq.createSession('ni');

%% input: use four inputs in differential mode
s.addAnalogInputChannel('Dev1',0:3,'Voltage'); %add channels ai0--3
for j = 1:nsensor
    s.Channels(j).Range=[-Vr Vr]; %set input range, all channels
end
lh_s = s.addlistener('DataAvailable',@getDataCont); %after fscan*tscan samples we call 'getDataCont'


%% output: used to power sensors (generate dc voltage of 'Vs')
% sOut.outputSingleScan(Vs);

% start measurements
s.startBackground(); 

% clean-up
% delete(lh_s);
% s.stop();
% s.release();
