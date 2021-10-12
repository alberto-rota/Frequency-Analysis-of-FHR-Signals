m = ar(h_clean{1},100);
sys = idtf(1,m.A);
fpe(sys)
wn = randn(length(h_clean{1}),1);
sigest = filter(1, m.A,wn);
close
hhfigure
plot(h_clean{1});
hold on
plot(sigest);
hold off