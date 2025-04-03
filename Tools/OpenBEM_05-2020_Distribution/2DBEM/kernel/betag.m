function betaDB=betag(sigma,f)

% betaDB=betag(sigma,f);
%
% Calculates the admittance by Delany and Bazley formula
%
% Input:
%    -sigma: flow resistance in Nsm-4. It may be a vector.
%    -f:     frequency. It may be a vector.
%
% Output variables:
%    -betaDB:  normalised admittance. As many rows as the length of sigma.
%              As many columns as the length of f.

% Reference: Propagation of Sound in Porous Media: Modelling Sound Absorbing Materials,
% Second Edition, Jean F. Allard, Noureddine Atalla. 2009, John Wiley.
% Section 2.5.3

% Modified with new expression, Vicente Cutanda Henr�quez 4-2011.


% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (�C)
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

betaDB = zeros(length(sigma), length(f));
for ff=1:length(f)
    for rr=1:length(sigma)
        if  isinf(sigma(rr))
            betaDB(rr,ff)=0;
        elseif isnan(sigma(rr))
            betaDB(rr,ff)=NaN;
        else
            %betaDB(rr,ff)=rho*c./(1+9.08*(1000*f(ff)./sigma(rr)).^(-0.75)+11.9i.*(1000*f(ff)./sigma(rr)).^(-0.73));
            X=rho*f(ff)./sigma(rr);
            betaDB(rr,ff)=1./(1+0.0571*(X.^-0.754)+1j*0.087*(X.^-0.732)); % normalized admittance, exp(-jwt) convention
            % betaDB(rr,ff)=betaDB(rr,ff)/(rho*c); % de-normalize
        end
    end
end
end