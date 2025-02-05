figure;
polarplot(the,abs(p_s).^2, 'LineWidth',1.); grid on
ax = gca;
% ax.ThetaTick = 0:45:315;
subtitle(ax,'\theta /° re. \theta_0','VerticalAlignment','bottom')

ax.RAxis.Label.String = '|p| (Pa)';

legend('k = 1','k = 2','k = 3','k = 4','k = 5');
title('Scattering pressure field along the cylinder surface');
