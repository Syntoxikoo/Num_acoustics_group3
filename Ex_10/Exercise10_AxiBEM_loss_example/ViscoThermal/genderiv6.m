function [D0,Dl] = genderiv6(rzb,topology,rzline)

% [D0,D1] = genderiv6(rzb,topology,rzline);
%
% Gradient and laplacian matrices including the whole generator and varying
% spacings. Based on the secant formula applied twice. It uses the node
% spacings along the generator (rzline) as produced by nodegen and the
% nodes matriz (rzb). It detects body numbers and performs backward and
% forward derivatives at the begining and end of each body.

% Version 3 uses mirror nodes at the ends and a simple formula.(VCH 04-2011)
% This is the simplest and best performing version so far. (VCH 22-4-2011)
% The same result and other possibilities can be done with "genderiv2"

% Vicente Cutanda Henriquez 1-2011

% VCH, 11-2012: Allows not assuming zero gradient at rho=0, and calculate
% forward and backward differences instead. Increases stability of the
% matrix.

% VCH, 5-2014: The possibility of bodies with a hole in the center and not
% connected to the z-axis should has been added. The D0 and Dl are not 
% calculated differently at the first and last nodes of such bodies. (VCH)

rho0D0=1;  % flag to indicate zero gradient at rho=0.

rzline=rzline(:,1);

M=size(rzb,1);

Dl=zeros(M,M);D0=zeros(M,M);
nbodies=max(abs(rzb(:,end))); % number of bodies


for bb=1:nbodies
    bbn=find(abs(rzb(:,end))==bb);bbnt=length(bbn); % find nodes in body bb
    bbe=find(abs(topology(:,end))==bb); % find nodes in body bb
    
    if abs(rzb(topology(bbe(1),1),1)-rzb(topology(bbe(end),end-1),1))<eps && abs(rzb(topology(bbe(1),1),2)-rzb(topology(bbe(end),end-1),2))<eps % Detect closed bodies (with a hole)
        rzdif=diff([0; rzline(bbn(2:end)); rzline(bbn(1))]); % Make the interval calculation similar to that in bodies with no holes
        for jj=1:bbnt
            iih=[(mod((jj)-2,length(rzdif))+1) (mod((jj)-1,length(rzdif))+1)];  %indexes of the spacings
            iin=[(mod((jj)-2,bbnt)+1) (mod((jj)-1,bbnt)+1) (mod((jj),bbnt)+1)]; %indexes of the nodes
            h1=rzdif(iih(1));h2=rzdif(iih(2));
            Dl(bbn(jj),bbn(iin))=[2/(h1^2+h1*h2) -2/(h1*h2) 2/(h2^2+h1*h2)]; % Centered difference
            D0(bbn(jj),bbn(iin))=[-1/(h1+h2) 0 1/(h1+h2)];                   % Centered difference
        end
    else % bodies with no holes, attached to the z-axis
        rzdif=diff(rzline);
        for jj=1:bbnt
            if jj==1 % first node of the body
                h1=rzdif(jj);h2=h1;
                Dl(bbn(jj),bbn(1:2))=[-2/(h1*h2) 4/(h2^2+h1*h2)];   % Centered difference using the symmetrical node
                %D0: The gradient should be zero for rho=0.
                if ~rho0D0, D0(bbn(jj),bbn(1:2))=[-1/h1 1/h1]; end  % Forward difference
                
            elseif jj==bbnt  % last node of the body
                h1=rzdif(jj-1);h2=h1;
                Dl(bbn(jj),bbn(end-1:end))=[4/(h1^2+h1*h2) -2/(h1*h2)];   % Centered difference using the symmetrical node
                %D0: The gradient should be zero for rho=0.
                if ~rho0D0, D0(bbn(jj),bbn(end-1:end))=[-1/h1 1/h1]; end  % Backward difference
            else
                h1=rzdif(jj-1);h2=rzdif(jj);
                Dl(bbn(jj),bbn(jj-1:jj+1))=[2/(h1^2+h1*h2) -2/(h1*h2) 2/(h2^2+h1*h2)]; % Centered difference
                D0(bbn(jj),bbn(jj-1:jj+1))=[-1/(h1+h2) 0 1/(h1+h2)];                   % Centered difference
            end
        end
    end
end
end
