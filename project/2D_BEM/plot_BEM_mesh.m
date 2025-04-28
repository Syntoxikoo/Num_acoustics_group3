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

for ii =1:size(topology_low,1)
    plot(xyb_low(topology_low(ii,1:3),1),xyb_low(topology_low(ii,1:3),2),"Color","black","Linewidth",1.5); 
    hold on;
    Leg(1) = plot(xyb_low(topology_low(ii,[1 3]),1),xyb_low(topology_low(ii,[1 3]),2),"Color", "#181748","Marker",".","LineStyle","none","MarkerSize",10);
    midpoint_x =xyb_low(topology_low(ii,2),1); midpoint_y = xyb_low(topology_low(ii,2),2); 
    % plot(midpoint_x,midpoint_y,"Color","#181748","Marker","*","LineStyle","none", "MarkerSize",4);
end
Leg(2) =plot(fpxy_low(:,1),fpxy_low(:,2),"Color", "#810100","Marker",".","LineStyle","none", "MarkerSize",10);
Leg(3) = plot(xyb_low(nn_low,1),xyb_low(nn_low,2), "Color", "#81819B", "Linewidth",2);
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
% Leg(4) =plot(xyb_chief(:,1),xyb_chief(:,2),"Color","black","Marker","square","Linestyle","none","MarkerSize", 4);
grid on ; hold off
title( M_low + " nodes, "+ N_low +" elements, " + Nfield_low + " Field p.")


tile2 = nexttile;

for ii =1:size(topology_high,1)
    plot(xyb_high(topology_high(ii,1:3),1),xyb_high(topology_high(ii,1:3),2),"Color","black","Linewidth",1.5); 
    hold on;
    plot(xyb_high(topology_high(ii,[1 3]),1),xyb_high(topology_high(ii,[1 3]),2),"Color", "#181748","Marker",".","LineStyle","none","MarkerSize",10);
    midpoint_x =xyb_high(topology_high(ii,2),1); midpoint_y = xyb_high(topology_high(ii,2),2); 
    % plot(midpoint_x,midpoint_y,"Color","#181748","Marker","*","LineStyle","none", "MarkerSize",4);
end
plot(fpxy_high(:,1),fpxy_high(:,2),"Color", "#810100","Marker",".","LineStyle","none", "MarkerSize",10);
plot(xyb_high(nn_high,1),xyb_high(nn_high,2), "Color", "#81819B", "Linewidth",2)
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);



plot(xyb_chief(:,1),xyb_chief(:,2),"Color","black","Marker","square","Linestyle","none","MarkerSize", 5);
xlabel(tiled, "x (m)")
ylabel(tiled, "y (m)")
title( M_high + " nodes, "+ N_high +" elements, "+ Nfield_high + " Field p.")
grid on;
hold off;

leg = legend(Leg, {"Mesh body", "Field point", "Velocity"}, 'NumColumns', 3); 
leg.Layout.Tile = 'north'; 

saveas(gcf,'project/figures/mesh_size.svg')



