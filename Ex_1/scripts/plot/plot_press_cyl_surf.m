figure;
plot(the,abs(p_s).^2, 'LineWidth',1.); grid on
legend('ka= 1','ka = 2','ka = 3','ka = 4','ka = 5');
xlabel('angle /°')
ylabel('|p|^2')
title('Scattering pressure field along the cylinder surface');
