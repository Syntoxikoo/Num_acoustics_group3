function FEM_ex8_1_DrawElementTopology(ex,ey,elnum)

%-------------------------------------------------------------
% PURPOSE 
%   Draw the undeformed 2D mesh
% INPUT  
%    ex,ey:.......... nen:   number of element nodes
%                     nel:   number of elements          
%    elnum=edof(:,1) ; i.e. the first column in the topology matrix
%         
%-------------------------------------------------------------


[nel nen]=size(ex);

 
%  plot coordinates 
x0=sum(ex')/nen; 
y0=sum(ey')/nen;
    
x=ex' ;
y=ey';   
xc=[x ; x(1,:)]; yc=[y ; y(1,:)];

%% plot commands 

 hold on
 axis equal 
 plot(xc,yc,'b') 
 plot(x,y,'b')

%% Element Numbering
for i=1:nel
    h=text(x0(i),y0(i),int2str(elnum(i)));
    set(h,'fontsize',9);
end

ss = get(0,'screensize');
set(gcf,'color','w','position',[30 200 ss(3)/1.5 ss(4)/1.5]);
title(['Element topology'],'FontName','Arial','FontSize',12)    
 hold off 

