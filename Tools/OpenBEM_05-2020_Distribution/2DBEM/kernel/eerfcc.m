function [output]= eerfcc(t)
%function [output]= eerfcc(t)
%INPUT
%t=iw=i(x+iy)
%OUTPUT
%output=erfc
%NOTE:it doesn't work with vectors

% Function til beregning af exp(-w^2)erfc(iw)
% for komplekse argumenter 
% efter forbillede i K.B. Rasmussen's "Lydudbredelse Udendørs" AT note 
% nr. 2206. APPENDIX: Pascal program på side 49. Funktionen vil sikkert
% kunne implementeres mere simpelt i MATLAB end gjort her, men det må 
% vente til en anden gang. Læg mærke til at t ofte vil være en vektor,
% hvilket betyder at man må passe på når man anvender if-else sætninger
%

output=[];
for z=t;
	A=abs(real(z));
	B=abs(imag(z));
	if ((A<=3.9)&(B<=3))
		AA=A*A;
		BB=B*B;
		PIN=2.5*pi;
		AB=2*A*B;
		A1=cos(AB);
		B1=sin(AB);
		C1=exp(-B*PIN)-cos(A*PIN);
		D1=sin(A*PIN);
		CC=C1*C1+D1*D1;
		DD=2*exp(-(AA+B*PIN-BB));
		P2=DD*(A1*C1-B1*D1)/CC;
		Q2=DD*(A1*D1+B1*C1)/CC;
		DD=0.8/(pi*(AA+BB));
		AH=DD*B;
		AK=DD*A;
		for N = 1:5;
			CC=0.64*N*N;
			AN=(BB-AA+CC)^2+4*AA*BB;
			DD=1.6/pi*exp(-CC)/AN;
			AH=AH+B*DD*(AA+BB+CC);
			AK=AK+A*DD*(AA+BB-CC);
		end
		if (B < 1.25*pi)
			AH=AH+P2;
			AK=AK-Q2;
		end;
		if (B==1.25*pi)
			AH=AH+P2/2;
			AK=AK-Q2/2;
		end;
		w=AH+i*AK;
	else
		w=0;
		C=0.4613135;
		D=0.1901635;
		w=suberf(A,B,C,D,w);
		C=0.09999216;
		D=1.7844927;
		w=suberf(A,B,C,D,w);
		C=0.002883894;
		D=5.5253437;
		w=suberf(A,B,C,D,w);
	end
	if real(z)<0
		w=conj(w);
		A=-A;
	end	
	if imag(z)<0
		DD=2*exp(B*B-A*A);
		w=DD*exp(2*A*B)-conj(w);
	end
	output=[output w];
end;

function [wo] = suberf (A,B,C,D,w)
% Underfunktion til rekursion i eerfcc.m

AA=A*A;
BB=B*B;
CC=AA-BB-D;
DD=2*A*B;
AN=CC*CC+DD*DD;
wr=C*B*(AA+BB+D)/AN;
wi=C*A*(AA+BB-D)/AN;
wo=(wr+i*wi)+w;
