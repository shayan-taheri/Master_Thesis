% In all the bad data calculation after the bad data is detected, it is
% kept in the record
% recDir = 'C:\records\nl\25\';%'C:\Documents and Settings\rgerdes\Desktop\2014-SPRING-01-NOLOG\01\';%'Z:\keyboard\1MS\2013-FALL-02-LOG\12\';
% mkdir(recDir);

%% keyboard daq setup
AI=analoginput('mcc',0);
out = daqhwinfo(AI)
chan = addchannel(AI,0);
% The total sample rate per channel
Sample_Rate = 1e6;

flag = 0;

set(AI,'SampleRate',Sample_Rate); % total sample rate is 1000000 and for each channel it is
% 500000 as there are two channels
ActualRate = get(AI,'SampleRate');
set(chan, {'InputRange'}, {[-5 5] });
% Samples per Trigger is defined
Samples_Per_Trigger = 350000;

%Amount of samples to be collected before triggering occurs
PreTrigger_Data = 2000;

%Counting bad records
bad_rec = 0;

% number of records to acquire
n = 10;

% Counter for counting the good records
j = 0;

set(AI,'SamplesPerTrigger',Samples_Per_Trigger) % 2 seconds of data after triggering
% so samples per trigger is kept double of sample rate
 
set(AI,'TriggerChannel',chan) % triggering is done at channel 1 where the square   
% wave from the function generator goes in

set(AI,'TriggerType','Software')
set(AI,'TriggerCondition','Falling')  
% set(AI,'TriggerConditionValue',4) %logic level shifter
set(AI,'TriggerConditionValue',1) %no logic level shifter
set(AI,'TriggerDelayUnits','Samples')
set(AI,'TriggerDelay',-PreTrigger_Data) 

%% temperature daq setup
fscan = 100; %scan rate (sampling , Hz)
tscan = 6000; %max scan time (sec): after tscan secs we average temp readings
Vr = 1; %input voltage range (+/- volts); =1 gives -100 to +100 C

% temp sensor vars
sf = 0.010; %volts per centrigrade scaling factor
Vs = 5.0; %sensor source voltage
nsensor = 4; %num of sensors

% data: num of samples = length of session / length of sample
global temp; 
temp = [];

% daq connection and setup
s = daq.createSession('ni');
s.Rate = fscan;
s.IsContinuous = true; %sample until stopped, unless we exceed tscan sec
s.NotifyWhenDataAvailableExceeds = fscan*tscan; %stop scanning after tscan sec

% input: use four inputs in differential mode
s.addAnalogInputChannel('Dev1',0:3,'Voltage'); %add channels ai0--3
for j = 1:nsensor
    s.Channels(j).Range=[-Vr Vr]; %set input range, all channels
end
lh_s = addlistener(s,'DataAvailable',@getDataCont); %after s.stop() or fscan*tscan samples we call 'getDataCont'

%% Directory where data is to be recorded
recDir = 'C:\Documents and Settings\saptarshi.mallick\Desktop\getTemp\temp_data\';%'C:\Documents and Settings\rgerdes\Desktop\2014-SPRING-01-NOLOG\01\';%'Z:\keyboard\1MS\2013-FALL-02-LOG\12\';
mkdir(recDir);

%% setup to control AWG
vu = visa('agilent', 'USB0::2391::1031::my44007623::0::INSTR'); %from instrhwinfo('visa','agilent')
pause(1);             
fopen(vu);   
fprintf(vu,'OUTP ON'); %enable output of awg

%% a while loop for getting the data until the amount of good data is captured
disp('Acquiring data...')
close all;
tic;
while j < n
    s.startBackground();
    start(AI); % Starting to Acquire the data
    
    %AI.TriggersExecuted; %Keeping track of the number of triggers executed
    %while AI.TriggersExecuted == 0
    %pause(2e-6); %assume daq fs=500e3 Hz to reduce the CPU usage
    %end   
    %disp('Waiting...');
    wait(AI,6000);
    %disp(['...triggered ' num2str(j)]);
    [S,time] = getdata(AI,Samples_Per_Trigger); % Getting the samples for 2 seconds with the data
    %stop(AI);
    
    s.stop(); %get average temperature over daq period for keypress
    
    d = datevec(now);
    
    clk = S(:,1); % sampling on the firs t input and recording the clock signal
    %dat = s(:,2); % sampling on the second input and recording the data
    
    plot(clk);
    % Checking whether there is bad data by checking whether the data and the
    % clock value goes negative once after the triggering occurs
    
    if chkData1(clk,PreTrigger_Data) == 1
        bad_rec = bad_rec + 1;
        badRecStr = num2str(bad_rec);
        f = [recDir 'bsample' badRecStr '.mat'];
        save(f, 'clk', 'd','temp','recDir');
    else
        j = j + 1;
        js = num2str(j);
        f = [recDir 'sample' ...
            strrep(num2str(zeros(1,5-length(js))), ' ', '') js '.mat'];
        save(f, 'clk', 'd','recDir','temp','recDir');
    end
    
    t = round(toc); % get current running time
    disp(['Record: ' num2str(j) '/' num2str(n) ... 
        ' (%' num2str(j/n*100) ...
        ', num. bad: ' num2str(bad_rec) '); ' ...
        'Elapsed time: ' num2str(floor(t/3600)) ' hr. ' ...
        num2str(mod(floor(t/60),60)) ' min. ' ...
        num2str(mod(t,60)) ' sec.']);
    
    % clear data
    clk = 0;
end

AI
delete(lh_s);
% turn off awg
fprintf(vu,'OUTP OFF'); %disable output on awg
fclose(vu);

% release temp daq
s.release();