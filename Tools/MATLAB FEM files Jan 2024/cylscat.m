function [ptot,pinc,pscat]=cylscat(kp,a,xy,nterm)
% [ptot,pinc,pscat]=cylscat(kp,a,xy,nterm);

% Calculates formula 29.2 in Morse
% for a list of field points around the cylinder.
% kp=wavenumber
% a=cylinder radius
% xy=field points, col 1 (x), col 2 (y)
% nterm=number of terms in the series
%
% ptot,pinc,pscat = total, incident and scattered pressures

% 31265 Numerical Acoustics - Vicente Cutanda Henriquez

% Change to cylindrical (polar) coordinates
R=abs(xy(:,1)+j*xy(:,2))';
phi=angle(xy(:,1)+j*xy(:,2))';

pinc=besselj(0,kp*R);
pscat=zeros(1,size(xy,1));

% Incident pressure
for m=1:nterm-1
   if m==0;em=1;else;em=2;end;
   pinc=pinc+2*i^(m)*besselj(m,kp*R).*cos(m*(phi));
end

% Scattered pressure
for m=0:nterm-1
   if m==0;em=1;else;em=2;end;
   pscat=pscat-em*i^(m+1)*exp(-i*pa(m,kp*a)).*sin(pa(m,kp*a)).*(besselj(m,kp*R)+i*bessely(m,kp*R)).*cos(m*(phi));
end

ptot=pinc+pscat;
end

function p=pa(mn,ka);
% Calculates formula 8.1.2b in Morse-Ingard

if mn==0
   p=atan(-besselj(1,ka)/bessely(1,ka));
else
   p=atan((besselj(mn-1,ka)-besselj(mn+1,ka))/(bessely(mn+1,ka)-bessely(mn-1,ka)));
end
end
   
