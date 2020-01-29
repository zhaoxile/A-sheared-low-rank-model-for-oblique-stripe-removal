function  power_spectrum(signal)

[sx sy]=size(signal);
P=0;
for i=1:sy
[Px,w] = periodogram(double(signal(:,i)));
P=P+Px;
end;
P=P/sy;
% ws=w;
% figure(1);plot((w./(2*max2(w))),log10(P1));hold;plot((w./(2*max2(w))),log10(P2),'red')
figure;plot((w./(2*max(w))),log10(P+1));
% ylim([0,10]);
xlabel('Normalized frequency');
ylabel('Power spectrum ');
grid on;
set(gca,'linewidth',1);
