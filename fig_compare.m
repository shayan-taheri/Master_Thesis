test_dir1 = '/media/SHAYAN_HDD/Results/Collection_1/test_res_fig/1/';
test_dir2 = '/media/SHAYAN_HDD/Results/Collection_2/test_res/1/';
fig_dir = '/media/SHAYAN_HDD/Results/Collection_2/fig_dir/';
figx = 0;
% *****************************************************
load ([test_dir1,'b5c1_b5c1.mat']);
figx=figx+1;
fig_id=figure(figx);
plot(m0.gen,'b');
hold on;
load ([test_dir2,'b5c1_b5c1.mat']);
plot(m0.gen,'r');
legend('Dataset 1','Dataset 2');
hold off;
saveas(fig_id,[fig_dir,'b5c1_comp.pdf']);
% *****************************************************
load ([test_dir1,'b5c2_b5c2.mat']);
figx=figx+1;
fig_id=figure(figx);
plot(m0.gen,'b');
hold on;
load ([test_dir2,'b5c2_b5c2.mat']);
plot(m0.gen,'r');
legend('Dataset 1','Dataset 2');
hold off;
saveas(fig_id,[fig_dir,'b5c2_comp.pdf']);
% *****************************************************
load ([test_dir1,'b5c3_b5c3.mat']);
figx=figx+1;
fig_id=figure(figx);
plot(m0.gen,'b');
hold on;
load ([test_dir2,'b5c3_b5c3.mat']);
plot(m0.gen,'r');
legend('Dataset 1','Dataset 2');
hold off;
saveas(fig_id,[fig_dir,'b5c3_comp.pdf']);
% *****************************************************
load ([test_dir1,'b5c4_b5c4.mat']);
figx=figx+1;
fig_id=figure(figx);
plot(m0.gen,'b');
hold on;
load ([test_dir2,'b5c4_b5c4.mat']);
plot(m0.gen,'r');
legend('Dataset 1','Dataset 2');
hold off;
saveas(fig_id,[fig_dir,'b5c4_comp.pdf']);
% *****************************************************
load ([test_dir1,'b5c5_b5c5.mat']);
figx=figx+1;
fig_id=figure(figx);
plot(m0.gen,'b');
hold on;
load ([test_dir2,'b5c5_b5c5.mat']);
plot(m0.gen,'r');
legend('Dataset 1','Dataset 2');
hold off;
saveas(fig_id,[fig_dir,'b5c5_comp.pdf']);