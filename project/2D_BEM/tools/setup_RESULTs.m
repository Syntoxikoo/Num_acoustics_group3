% field points
fp_theta_plot = fp_theta*pi/180;
Rmed = 1;
fpxy=[Rmed*sin(fp_theta_plot) Rmed*cos(fp_theta_plot)];


nrows = 1;
ncols = 3;
heightScale = 0.5; % Adjust height scaling if needed
Nleg = 1;

[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

tile1 = nexttile;

for ii =1:size(topologyUF,1)
    Leg(1) = plot(xybUF(topologyUF(ii,1:3),1),xybUF(topologyUF(ii,1:3),2),"Color","black","Linewidth",1.5); 
    hold on 
end
Leg(2) = plot(fpxy(:,1),fpxy(:,2),"Color", "#810100","Marker",".","LineStyle","none", "MarkerSize",10);

xlim([-1.2 1.2]);
ylim([-1.2 1.2]);

grid on ; hold off
title( "Recessed 10 cm")


tile2 = nexttile;

for ii =1:size(topologyR,1)
    plot(xybR(topologyR(ii,1:3),1),xybR(topologyR(ii,1:3),2),"Color","black","Linewidth",1.5); 
    hold on;
end
plot(fpxy(:,1),fpxy(:,2),"Color", "#810100","Marker",".","LineStyle","none", "MarkerSize",10);
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);



xlabel(tiled, "x (m)")
ylabel(tiled, "y (m)")
title( "Rounded")
grid on;
hold off;

tile3 = nexttile;

for ii =1:size(topologySl,1)
    plot(xybSl(topologySl(ii,1:3),1),xybSl(topologySl(ii,1:3),2),"Color","black","Linewidth",1.5); 
    hold on;
end
plot(fpxy(:,1),fpxy(:,2),"Color", "#810100","Marker",".","LineStyle","none", "MarkerSize",10);
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);



xlabel(tiled, "x (m)")
ylabel(tiled, "y (m)")
title( "Slit")
grid on;
hold off;

leg = legend(Leg, {"Meshed body", "Observation point"}, 'NumColumns', 3); 
leg.Layout.Tile = 'north'; 

saveas(gcf,'figures/RESULTs_setup.svg')



