function SPLout=SPL(pres);

% Obtains the Sound Pressure Levels of a matrix or vector of complex
% harmonic pressures.
% plot the SPL relative to free field versus frequency

% Vicente Cutanda 2007

pref=20e-6; % Reference pressure, 20 microPa
SPLout=20*log10(abs(pres)/sqrt(2)/pref);
