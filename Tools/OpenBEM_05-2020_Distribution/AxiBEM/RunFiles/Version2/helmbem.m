function [phi,SVals,nchiefp,threshold]=helmbem(k,rzb,incoming,varargin)

%  [phi]=helmbem(k,rzb,pI)
%  [phi]=helmbem(k,rzb,pI,v)
%
%  Axisymetric BEM calculation.
%  Scattering of an axisymmetric incident field.
%  Radiation in testfase
%  
%  Input parameters (in SI units):
%  
%    -k       : Wavenumbers (k=2*pi*f/c) (real vector)
%    -rzb     : Geometry matrix
%               one row for each node
%               rho-coordinate in column 1
%               z-coordinate in column 2
%               body number in column 3 (not necessary if there is
%               just one body)
%    -incoming: function that computes the incoming field from 'rzb'
%               and a wavenumber
%               Usually of type 'inline'
%               Alternatively the name of a function
%               [] or '' if there is no incoming field
%    -dphidn  : normal derivative of phi (complex column vector)
%               normal velocity of boundary = -1 * dphidn
%               Alternatively, one column can be given for each k-value
%
%  Output parameters:
%  
%    -phi     : complex matrix, the solution along the generator.
%               Each column contains the solution in the nodes
%               for the corresponding k-value

%  This program is based on:
%  Peter Juhl "The Boundary Element Method for Sound Field Calculations"
%  Report no. 55 1993, The Acoustics Laboratory,
%  Technical University of Denmark
%
%  It was written by:
%    pmj - Peter Møller Juhl, Institute of Applied Physics,
%          Odense University
%    vc  - Vicente Cutanda, Brüel & Kjær, Nærum and
%          Department of Acoustic Technology, Lyngby
%    msj - Morten Skaarup Jensen, Department of Acoustic Technology,
%          Technical University of Denmark

%  File history:
%  msj 990906 singular part of integral nolonger calculated for every frequency
%  pmj 990129 added some support for radiation
%  msj 981210
%  

%  Key points about input
[M, ncolrzb] = size(rzb);
if ncolrzb<3
	rzb=[rzb ones(M,1)];
end
NumBodies = rzb(M,3);

inode=0;
iel=0;
for ibody=1:NumBodies
  inode=inode+1;
  while inode~=M & rzb(inode+1,3)==ibody
     iel=iel+1;
     if rzb(inode,:)==rzb(inode+1,:)
        inode=inode+1;
     end
     ElemNodeNum(iel,1)=inode;
     ElemNodeNum(iel,2)=inode+1;
     ElemNodeNum(iel,3)=inode+2;
     ElemNodeNum(iel,4)=ibody;
     inode=inode+2;
  end
end
N=iel;

% If 'incoming' is not an inline function, turn it into one
if isempty(incoming)
   incoming='';
end
if ~isa(incoming,'inline')
   if isempty(incoming) | incoming==''
      incoming=inline('zeros(size(rzb,1),1)*wavenum');
   else
      incoming=inline([incoming '(rzb,k)']);
   end
end

% Circumferential mode number
m=0;
% Impedance (dummy at the moment)
beta=zeros(size(rzb,1));

% Neumann boundary condition
if nargin<4
   dphidn=zeros(size(rzb,1));
else
   dphidn=varargin{1};
end

%  Calculate C constants and the singular part of 'A' and 'B'
[A0,B0,CConst]=SteadEquat(rzb,ElemNodeNum,beta,m);
SVals0=svd(A0);
SVals0=SVals0(length(SVals0))

l_max=max(rzb(:,1))-min(rzb(:,1))+max(rzb(:,2))-min(rzb(:,2));
l_char=sqrt((max(rzb(:,1))-min(rzb(:,1)))*(max(rzb(:,2))-min(rzb(:,2))));
nchiefp=zeros(size(k));
threshold=zeros(size(k));
for ik=1:length(k)
   disp(sprintf('Frequency %3d of %4d',ik,length(k)));
   kk=k(ik);
   % Calculate 'A' and 'B'
   if kk==0
      A=A0; B=B0;
   else
      [A,B]=HarmEquat(rzb,ElemNodeNum,kk,beta,m);
      A=A+A0;  B=B+B0;
   end

   % Compute the incoming field in the nodes
   pI=incoming(rzb,kk);

   % Check for non-uniqueness problem
   % If necessary correct using CHIEF
   % See Juhl 1994 J. Sound & Vib. 75(1), 39-50
   % "A Numerical Study of the Coefficient Matrix of the Boundary
   %  Element Method Near Characteristic Frequencies"
   % N.B. Can be further improved as suggested by Juhl
   % Threshold values according to
   % Jensen 1999, to be submitted to J. Sound & Vib.
   % "Solving the BEM non-uniqueness Problem with CHIEF"
   ik_ref=find((kk-k)*l_max < 1 & kk > k);
   if ~isempty(ik_ref)
      weights=1-(kk-k(ik_ref))*l_max;
      weights=weights./sum(weights);
      threshold(ik)=10^(sum(log10(SVals(ik_ref))' .*weights)-0.1*(kk-sum(k(ik_ref).*weights))*l_max);
   elseif ik~=1
      ik_ref=ik-1;
      threshold(ik)=10^(log10(SVals(ik_ref))-0.1*(kk-k(ik_ref))*l_max);
   else
      threshold(ik)=SVals0/(0.1*kk*l_max+1.1);
   end
   while 1
      singvals=svd(A);
      SVals(1,ik)=singvals(length(singvals));
      disp(sprintf('A matrix lowest singular value = %g',SVals(ik)));
      nchiefp(ik)=size(A,1)-size(A,2);
      if (ik==1 & SVals(ik) > threshold(ik) & nchiefp(ik)>=k*l_char/pi ) ...
       | (ik~=1 & SVals(ik) > threshold(ik))
         break
      end
      % Generate a random CHIEF point
      [Ax,Bx,chiefpoint]=chief(rzb,ElemNodeNum,kk,beta,m);
      A=[A; Ax];
      B=[B; Bx];
      pI=[pI; incoming(chiefpoint,kk)];
   end

   %  Solve system of equations:
   if size(dphidn,2)==1
      phi(:,ik)=A\(B*dphidn-4*pi*pI);
   else
      phi(:,ik)=A\(B*dphidn(:,ik)-4*pi*pI);
   end
end





