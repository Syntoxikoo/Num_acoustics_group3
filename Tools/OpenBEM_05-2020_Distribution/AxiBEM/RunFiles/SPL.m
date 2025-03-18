function pout=SPL(pin,varargin);

% pout=SPL(pin,direct_inv)
%
% Calculation of Sound Pressure Levels from Sound Pressure in Pascals, and
% viceversa.
%
% direct_inv=1: Obtains the Sound Pressure Levels of a matrix or vector of complex
%               harmonic pressures in Pascals.
% direct_inv=2: Obtains the absolute value of Sound Pressure in Pascals from a matrix
%               or vector of SPL values.

% Vicente Cutanda 2009

pref=20e-6; % Reference pressure, 20 microPa

if nargin>1
    if varargin{1}==1
        pout=20*log10(abs(pin)/sqrt(2)/pref);
    elseif varargin{1}==2
        pout=10.^(pin/20)*sqrt(2)*pref;
    else
        error('Invalid input in SPL function')
    end
else
    pout=20*log10(abs(pin)/sqrt(2)/pref);
end

