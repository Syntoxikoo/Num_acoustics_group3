nrows = 1;
ncols = 2;
heightScale = 0.55; % Adjust height scaling if needed
Nleg = 1;

[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

tile1 = nexttile;

for ii =1:size(topology,1)
    Leg(1) = plot(xyb(topology(ii,1:3),1),xyb(topology(ii,1:3),2),"Color","black","Linewidth",1.5); 
    hold on 
end
Leg(2) = plot(fpxy(1),fpxy(2),"Color", "#810100","Marker",".","LineStyle","none", "MarkerSize",10);

xlim([-1.2 1.2]);
ylim([-1.2 1.2]);

grid on ; hold off
title( "Recessed 10 cm")


tile2 = nexttile;

for ii =1:size(topology2,1)
    plot(xyb2(topology2(ii,1:3),1),xyb2(topology2(ii,1:3),2),"Color","black","Linewidth",1.5); 
    hold on;
end
plot(fpxy(1),fpxy(2),"Color", "#810100","Marker",".","LineStyle","none", "MarkerSize",10);
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);



xlabel(tiled, "x (m)")
ylabel(tiled, "y (m)")
title( "Recessed 15 cm")
grid on;
hold off;

leg = legend(Leg, {"Meshed body", "Observation point"}, 'NumColumns', 3); 
leg.Layout.Tile = 'north'; 

saveas(gcf,'figures/geometry_sensitivity_study_recessing.svg')



