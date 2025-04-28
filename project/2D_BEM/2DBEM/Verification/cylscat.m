function Int=cylscat(ka,ndiv,nterm,see);
% [Int]=cylscat(ka,ndiv,nterm);
% Calculates formula 8.1.2a in Morse-Ingard
% ka=wavenumber*radius
% ndiv=angular divisions in 2*pi (40?)
% nterm=number of terms in the series (nterm^2) (10?)
% see: if 'y', the result is plotted

% Vicente Cutanda 3-2001.

%ndiv=40;
%nterm=10

phi=0:2*pi/ndiv:2*pi;

IntI=ones(1,ndiv+1)*besselj(0,ka);
IntS=zeros(1,ndiv+1);

for m=1:nterm-1
   if m==0;em=1;else;em=2;end;
   IntI=IntI+2*i^(m)*besselj(m,ka).*cos(m*(phi-pi));
end

for m=0:nterm-1
   if m==0;em=1;else;em=2;end;
   IntS=IntS-em*i^(m+1)*exp(-i*pa(m,ka))*sin(pa(m,ka))*(besselj(m,ka)+i*bessely(m,ka)).*cos(m*(phi-pi));
end

Int=IntI+IntS;

if see=='y'|see=='Y'
    figure;
    plot(phi,abs(Int));
    title(['ka = ' num2str(ka)]);
    axis([0 2*pi 0 2])
    grid;
end

function p=pa(mn,ka);
% Calculates formula 8.1.2b in Morse-Ingard

if mn==0
   p=atan(-besselj(1,ka)/bessely(1,ka));
else
   p=atan((besselj(mn-1,ka)-besselj(mn+1,ka))/(bessely(mn+1,ka)-bessely(mn-1,ka)));
end
   