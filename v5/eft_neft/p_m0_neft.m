% plot filter output of ineffecitve tests

neft = tind(a_t/(nc-ntrn)<.9)

for i = 1:length(neft);
    figure;
    tst = neft(i);
    plot(m0_fp(:,tst));
    hold on;
    plot(m0_fn(:,tst),'r');
    title(['test number ' num2str(tst)]);
end