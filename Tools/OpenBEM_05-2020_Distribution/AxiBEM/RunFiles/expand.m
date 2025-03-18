function [a,b]=expand(x,numterms,varargin)

% [a,b]=expand(x,numterms,see)
% 
% Performs FFT on the sequence x and formats the output
% as coefficients of the Fourier series as defined in 
% Peter Juhl's thesis (p.31).
% 
% The 'a' are the coefficients of the cosine terms, and
% the 'b' are the coefficients of the sine terms. AxiBEM 
% uses cosine expansions only, but the sine part can be computed
% as cosine and later shifted by pi/2 and superposed.
% 
% If the input x should be even, x(i)=x(N-i), only
% cosine terms are obtained. N is x's length.
% 
% The number of terms in the series is 'numterms' and can be
% up to half the lenght of the x sequence.
% 
% Some plots are presented in a figure to check the
% accuracy of the terms requested (Set 'see' to 'yes').

% Vicente Cutanda 2000

% VCH 2009: produces full Fourier expansions in cosines and sines and
% takes complex sequences.
% http://en.wikipedia.org/wiki/Fourier_series


if nargin>2
   see=varargin{1};
else
   see='n';
end

if size(x,1)==1
   x=x';
end

% normalized FFT
y=fft(x)/length(x);

if ~exist('numterms')|numterms>=length(y)
   numterms=length(y)/2-1;
end

% Convert the FFT coefficients into sine/cosine expansion coefficients:
fplot=[1:numterms]';
a=[y(1) ; (y(fplot+1) + y(end-fplot+1))];
b=[0 ; j*(y(fplot+1) - y(end-fplot+1))];


if see(1)=='y' | see(1)=='Y'
   phi=0:2*pi/(length(x)-1):2*pi;

   xF=a(1) + a(2:end).'*cos(fplot*phi) + b(2:end).'*sin(fplot*phi);
   xC=a(1) + a(2:end).'*cos(fplot*phi);
   
   figure;
   subplot(2,3,1)
   plot(180*phi/pi,abs(xF),'b')
   hold on;
   plot(180*phi/pi,abs(x),'r')
   title(['Fourier expansion, ' num2str(numterms) ' terms']);xlabel('Theta (deg.)');ylabel('Modulus');
   hold off;
   legend('Reconstructed','Original')
   
   subplot(2,3,4)
   plot(180*phi/pi,abs(xC),'b')
   hold on;
   plot(180*phi/pi,abs(x),'r')
   title(['Cosine expansion, ' num2str(numterms) ' terms']);xlabel('Theta (deg.)');ylabel('Modulus');
   hold off;
   
   subplot(2,3,2)
   plot(180*phi/pi,180/pi*unwrap(angle(xF)),'b')
   hold on;
   plot(180*phi/pi,180/pi*unwrap(angle(x)),'r')
   title(['Fourier expansion, ' num2str(numterms) ' terms']);xlabel('Theta (deg.)');ylabel('Phase (deg.)');
   hold off;
   
   subplot(2,3,5)
   plot(phi/pi,180/pi*unwrap(angle(xC)),'b')
   hold on;
   plot(phi/pi,180/pi*unwrap(angle(x)),'r')
   title(['Cosine expansion, ' num2str(numterms) ' terms']);xlabel('Theta (deg.)');ylabel('Phase (deg.)');
   hold off;

   subplot(2,3,3)
   plot(real(a),'b')
   hold on;
   plot(imag(a),'r')
   title('a coefficients');
   legend('Real','Imag')
   hold off;

   subplot(2,3,6)
   plot(real(b),'b')
   hold on;
   plot(imag(b),'r')
   title('b coefficients');
   legend('Real','Imag')
   hold off;
   
%    subplot(2,3,3)
%    plot(real(y(fplot)),'b')
%    hold on;
%    plot(imag(y(fplot)),'r')
%    title('Real and imaginary coeff.');
%    hold off;
%    
%    subplot(2,3,6)
%    plot(abs(y(fplot)),'b')
%    title('Absolute value coeff.');

end
