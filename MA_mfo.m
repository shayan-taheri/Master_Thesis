% ***** Finding Moving Average of Matched Filter Output *****
clear all;
close all;
% ---- Preparing Environment ----
wind_size = 2000; % Window Size for Moving Average
test_dir = '/media/SHAYAN_HDD/Results/Collection_3/test_res/1/';
fig_dir = '/media/SHAYAN_HDD/Results/Analysis/fig_dir/col3/';
% *****************************************************
load ([test_dir,'b5c1_b5c1.mat']);

MA_out = tsmovavg(m0.gen,'s',wind_size);
nan_ind = isnan(MA_out);
MA_out = MA_out(~nan_ind);
x_first = median(1:wind_size);
x_axis = x_first:(x_first+length(MA_out)-1);
fig_id = figure();
plot(m0.gen);
hold on;
plot(x_axis,MA_out,'r');

saveas(fig_id,[fig_dir,'b5c1_',num2str(wind_size),'.pdf']);
% *****************************************************
load ([test_dir,'b5c2_b5c2.mat']);

MA_out = tsmovavg(m0.gen,'s',wind_size);
nan_ind = isnan(MA_out);
MA_out = MA_out(~nan_ind);
x_first = median(1:wind_size);
x_axis = x_first:(x_first+length(MA_out)-1);
fig_id = figure();
plot(m0.gen);
hold on;
plot(x_axis,MA_out,'r');

saveas(fig_id,[fig_dir,'b5c2_',num2str(wind_size),'.pdf']);
% *****************************************************
load ([test_dir,'b5c3_b5c3.mat']);

MA_out = tsmovavg(m0.gen,'s',wind_size);
nan_ind = isnan(MA_out);
MA_out = MA_out(~nan_ind);
x_first = median(1:wind_size);
x_axis = x_first:(x_first+length(MA_out)-1);
fig_id = figure();
plot(m0.gen);
hold on;
plot(x_axis,MA_out,'r');

saveas(fig_id,[fig_dir,'b5c3_',num2str(wind_size),'.pdf']);
% *****************************************************
load ([test_dir,'b5c4_b5c4.mat']);

MA_out = tsmovavg(m0.gen,'s',wind_size);
nan_ind = isnan(MA_out);
MA_out = MA_out(~nan_ind);
x_first = median(1:wind_size);
x_axis = x_first:(x_first+length(MA_out)-1);
fig_id = figure();
plot(m0.gen);
hold on;
plot(x_axis,MA_out,'r');

saveas(fig_id,[fig_dir,'b5c4_',num2str(wind_size),'.pdf']);
% *****************************************************
load ([test_dir,'b5c5_b5c5.mat']);

MA_out = tsmovavg(m0.gen,'s',wind_size);
nan_ind = isnan(MA_out);
MA_out = MA_out(~nan_ind);
x_first = median(1:wind_size);
x_axis = x_first:(x_first+length(MA_out)-1);
fig_id = figure();
plot(m0.gen);
hold on;
plot(x_axis,MA_out,'r');

saveas(fig_id,[fig_dir,'b5c5_',num2str(wind_size),'.pdf']);
% -----------------------------------------------------