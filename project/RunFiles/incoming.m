function pI=incoming(sources,Dum1,Dum2,Dum3,varargin)

% pI=incoming(sources,NodeRhoZ,k,m)
%
% pI=incoming(sources,RR,ZZ,k,m)
%
% Computes the incoming field on the nodes of the generator.
% 
%  Input parameters (in SI units):
%  
%    -sources  : Sources information matrix. Each row is a source.
%                 Column 1: angle from z+ axis (0 <= angle <= PI).
%                 Column 2: 0 if plane wave, else point source's
%                           distance from origin (>0).
%                 Column 3: wave amplitude.
%                 Column 4: initial phase (point source) or
%                           phase at origin (plane wave).
%    -NodeRhoZ : Node coordinates along the generator. Matrix with
%                M rows, first column is the node rho coodinate,
%                second column is its z coordinate and third column
%                is the number of the body it belongs to.
%    -RR,ZZ    : Node coordinates in matrix form, as produced by "meshgrid"
%    -k        : Wavenumber (k=2*pi*f/c)
%    -m        : Term in the cosine expansion to be calculated.
%
%  Output parameters:
%  
%    -pI       : Incoming field on the nodes. It is a complex vector if
%                the node coordinates are supplied as NodeRhoZ, or as matrix
%                if RR and ZZ are inputs.
%                A vector of zeros with lenght the number of nodes
%                if the "sources" variable is empty.

% Vicente Cutanda 2000
% See Ph.D thesis by Peter M. Juhl, equation 3.6

% VCH 02-2009: meshgrid input

% Future improvements:
% -sources outside the rho-z plane
% -allow k and m as vectors and output a 3D matrix

if nargin>4
    [ff,cc]=size(Dum1);
    NodeRhoZ=[Dum1(1:end)' Dum2(1:end)'];
    k=Dum3;
    m=varargin{1};
else
    NodeRhoZ=Dum1;
    k=Dum2;m=Dum3;
end    


% Select integration formula
if m<2
   [bp,wf]=gaussrule(10);
else
   % For large 'm' values it is better to use the midpoint formula
   nip=8*m+8;
   bp=linspace(-1+1/nip,1-1/nip,nip)';
   wf=2*ones(nip,1)./nip;
end

% Short names
srs=sources;
nrz=NodeRhoZ(:,1:2);

if isempty(srs)
   pI=zeros(size(nrz,1),1);
else
   % Integration parameters
   thmin=0;
   thmax=pi;
   d=(thmax-thmin)/2;
   th=d*bp'+thmin+d;

   % dummy variables
   th11=ones(size(th));
   nrz11=ones(size(nrz,1),1);

   % Integrand
   pItmp=zeros(size(nrz,1),length(th));
   for ss=1:size(srs,1)
      if srs(ss,2)==0
         dps=nrz(:,1)*cos(th)*sin(srs(ss,1))+nrz(:,2)*th11*cos(srs(ss,1));
      
         pItmp=pItmp+nrz11*cos(m*th)*srs(ss,3).*...
            exp(i*(k*dps+srs(ss,4)));
      else
         dps=sqrt((nrz(:,1)*th11).^2+(srs(ss,2)*sin(srs(ss,1)))^2+...
            (srs(ss,2)*cos(srs(ss,1))-nrz(:,2)*th11).^2-...
            2*nrz(:,1)*srs(ss,2)*sin(srs(ss,1))*cos(th));
      
         pItmp=pItmp+nrz11*cos(m*th)*srs(ss,3)./dps.*...
            exp(-i*(k*dps+srs(ss,4)));
      end
   end

   % Integration, times 1 if m=0 or times 2 for the rest
   pI=d*(pItmp*wf)*(1+sign(m))/pi;
end

if nargin>4
    pI=reshape(pI,ff,cc);
end
