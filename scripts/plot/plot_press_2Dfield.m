figure(1);
surf(X,Y,abs(p_tot)); hold on;
% 3D cylinder
[Xc, Yc, Zc] = cylinder(a, 100);
Zc = Zc * max(max(abs(p_tot)));
Zc(Zc==0) = min(min(abs(p_tot)));
surf(Xc, Yc, Zc); hold off;