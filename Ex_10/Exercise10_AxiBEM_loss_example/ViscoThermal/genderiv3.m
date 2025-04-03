function [D0,Dl] = genderiv3(rzb,rzline)

% [D0,D1] = genderiv3(rzb,rzline);
%
% Gradient and laplacian matrices including the whole generator and varying
% spacings. Based on the secant formula applied twice. It uses the node
% spacings along the generator (rzline) as produced by nodegen and the
% nodes matriz (rzb). It detects body numbers and performs backward and
% forward derivatives at the begining and end of each body.

% This version (3) uses mirror nodes at the ends and a simple formula.(VCH 04-2011)
% This is the simplest and best performing version so far. (VCH 22-4-2011)
% The same result and other possibilities can be done with "genderiv2"

% Vicente Cutanda Henriquez 1-2011

rzline=rzline(:,1);

M=size(rzb,1);

Dl=zeros(M,M);D0=zeros(M,M);

%rzdif=sqrt(diff(rzb(:,1)).^2+diff(rzb(:,2)).^2);
rzdif=diff(rzline);
for jj=2:M-1
    if rzb(jj,end)~=rzb(jj+1,end) % last node of a body
        h1=rzdif(jj);h2=h1;
        Dl(jj,jj-1:jj)=[4/(h1^2+h1*h2) -2/(h1*h2)]; 
        %D0: The gradient should be zero for rho=0.
    elseif rzb(jj,end)~=rzb(jj-1,end) % first node of the next body
        h1=rzdif(jj);h2=h1;
        Dl(jj,jj:jj+1)=[-2/(h1*h2) 4/(h2^2+h1*h2)]; % Forward difference
        %D0: The gradient should be zero for rho=0.
    else % remaining nodes
        h1=rzdif(jj-1);h2=rzdif(jj);
        Dl(jj,jj-1:jj+1)=[2/(h1^2+h1*h2) -2/(h1*h2) 2/(h2^2+h1*h2)]; % Centered difference
        D0(jj,jj-1:jj+1)=[-1/(h1+h2) 0 1/(h1+h2)];                   % Centered difference
    end
end

% first node of the setup
h1=rzdif(1);h2=h1;
Dl(1,1:2)=[-2/(h1*h2) 4/(h2^2+h1*h2)]; % Forward difference
%D0: The gradient should be zero for rho=0.

% last node of the setup
h1=rzdif(end);h2=h1;
Dl(end,end-1:end)=[4/(h1^2+h1*h2) -2/(h1*h2)]; % Backward difference
%D0: The gradient should be zero for rho=0.

end
