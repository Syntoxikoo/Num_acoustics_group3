function plotSPLff(SPLff,f);

% plotSPLff(SPLff);
% plot the SPL relative to free field versus frequency

figure;
for pp=1:size(SPLff,2)
   semilogx(f,SPLff(:,pp),'-o');
   hold on
end
xlabel('Frequency');
ylabel('SPL relative to free field');
grid;
