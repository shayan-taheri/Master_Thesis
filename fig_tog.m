test_dir1 = '/media/SHAYAN_HDD/Results/Collection_1/test_res_fig/1/';
test_dir2 = '/media/SHAYAN_HDD/Results/Collection_2/test_res/1/';
fig_dir = '/media/SHAYAN_HDD/Results/Collection_2/fig_dir/';
figx = 0;
% *****************************************************
load ([test_dir1,'b5c1_b5c1.mat']);
figx=figx+1;
fig_id=figure(figx);
subplot(2,1,1);
plot(m0.gen,'b');
title('Dataset 1','FontSize',18);
load ([test_dir2,'b5c1_b5c1.mat']);
subplot(2,1,2);
plot(m0.gen,'r');
title('Dataset 2','FontSize',18);
saveas(fig_id,[fig_dir,'b5c1_tog.pdf']);
% *****************************************************
load ([test_dir1,'b5c2_b5c2.mat']);
figx=figx+1;
fig_id=figure(figx);
subplot(2,1,1);
plot(m0.gen,'b');
title('Dataset 1','FontSize',18);
load ([test_dir2,'b5c2_b5c2.mat']);
subplot(2,1,2);
plot(m0.gen,'r');
title('Dataset 2','FontSize',18);
saveas(fig_id,[fig_dir,'b5c2_tog.pdf']);
% *****************************************************
load ([test_dir1,'b5c3_b5c3.mat']);
figx=figx+1;
fig_id=figure(figx);
subplot(2,1,1);
plot(m0.gen,'b');
title('Dataset 1','FontSize',18);
load ([test_dir2,'b5c3_b5c3.mat']);
subplot(2,1,2);
plot(m0.gen,'r');
title('Dataset 2','FontSize',18);
saveas(fig_id,[fig_dir,'b5c3_tog.pdf']);
% *****************************************************
load ([test_dir1,'b5c4_b5c4.mat']);
figx=figx+1;
fig_id=figure(figx);
subplot(2,1,1);
plot(m0.gen,'b');
title('Dataset 1','FontSize',18);
load ([test_dir2,'b5c4_b5c4.mat']);
subplot(2,1,2);
plot(m0.gen,'r');
title('Dataset 2','FontSize',18);
saveas(fig_id,[fig_dir,'b5c4_tog.pdf']);
% *****************************************************
load ([test_dir1,'b5c5_b5c5.mat']);
figx=figx+1;
fig_id=figure(figx);
subplot(2,1,1);
plot(m0.gen,'b');
title('Dataset 1','FontSize',18);
load ([test_dir2,'b5c5_b5c5.mat']);
subplot(2,1,2);
plot(m0.gen,'r');
title('Dataset 2','FontSize',18);
saveas(fig_id,[fig_dir,'b5c5_tog.pdf']);