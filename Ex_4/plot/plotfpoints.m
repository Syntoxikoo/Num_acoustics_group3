function plotfpoints(XX,YY,q,varargin);

% plotfpoints(XX,YY,pfield,animate);
%
% calculation of the pressure on field points
%
% Input variables:
%     -XX:       matrix of x-coordinates, output from 'meshgrid' of the script.
%     -YY:       matrix of y-coordinates, output from 'meshgrid' of the script.
%     -q:        vector with the variable to represent, for the same coordinates
%                as above, but extended with (1:end).
%     -animate:  If 1, the plot is animated using the phase. The values of
%                q must therefore be complex. If 0 (default), a regular
%                plot is made and q should be real. The animation may be 
%                stopped by pressing CTRL-C.
%

% Vicente Cutanda Henriquez 2-2015

[XYr,XYc]=size(XX);

qM=zeros(XYr,XYc);
for cc=1:XYc
   for rr=1:XYr
      qM(rr,cc)=q((cc-1)*XYr+rr);
   end
end

if nargin>3
    % Animated plot of the modulus of q
    hh=figure  % screen animation
    phase=0;
    colorrange=[-max(abs(q)) max(abs(q))];
    while 1 % Execution may be stopped using Ctr-C
        clf;figure(hh)
        phase=phase+pi/180*10;    if phase>=2*pi-eps, phase=0;end;
        surf(XX,YY,real(qM*exp(-j*phase)));hold on
        shading interp
        caxis(colorrange)
        colorbar;%('vert');
        view(0,90);
        title('Press Ctrl-C to stop')
        drawnow;  %   pause(1)
    end
    close(hh)
else
    % Static plot of the modulus of q
    surf(XX,YY,qM);
    shading interp
    colorbar;%('vert');
    view(0,90);
end

end
