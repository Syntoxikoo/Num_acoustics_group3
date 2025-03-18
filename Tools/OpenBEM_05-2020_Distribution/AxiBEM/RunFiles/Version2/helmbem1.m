function [pm,NodeRhoZ,CondNum,CConst]=helmbem1(fs,file,sources,mter,Y,v,varargin)

%  [pm,NodeRhoZ,CondNum,CConst]=helmbem1(fs,file,sources,mter,Y,v,filename)
%
%  BEM calculation for the scattering by an axisymmetric object
%  with non-axisymmetric sound sources.
%  
%  Input parameters (in SI units):
%  
%    -fs      : Frequencies (real vector) in Hz.
%    -file    : Output file from NodeGen.exe with node coordinates
%    -sources : Sources information matrix. Each row is a source.
%               Column 1: angle from z+ axis (0 <= angle <= PI).
%               Column 2: 0 if plane wave, else point source's
%                         distance from origin (>0).
%               Column 3: wave amplitude.
%               Column 4: initial phase (point source) or
%                         phase at origin (plane wave).
%               If there are no sources, it should be empty, []
%    -mter    : Terms in the cosine expansion to be calculated
%               (integer vector).
%    -Y       : Admittances on the nodes. Vector, length M (number
%               of nodes).
%    -v       : Velocities on the nodes. Matrix, M (number
%               of nodes) rows and m (number of expansion terms)
%               columns.
%    -filename: file name to be used for storing the BEM coefficient
%               matrices on the disk temporally.
%               The variables in the file are named Am$ and Bm$,
%               where $ is the index in the m vector.
%               The file also contains the variables 'm' and 'k'. 
%               If set to [], no file is used, just memory.
%
%  Output parameters:
%  
%    -pm      : Complex matrix, the solution along the generator.
%               Each column contains the solution in the nodes
%               for the corresponding frequency. It is expanded in
%               the third dimension to hold solutions for every term
%               of the cosine expansion requested in 'mter'.
%    -NodeRhoZ: Node coordinates along the generator. Matrix with
%               M rows, first column is the node rho coodinate,
%               second column is its z coordinate and third column
%               is the number of the body it belongs to, with a
%               minus sign if the interior domain is specified.
%    -CondNum : Condition numbers of the A matrix for every
%               frequency (rows)and expansion term (columns).
%    -CConst  : C constants on every node (real vector).

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

% Physical constants
c=345.773;  % Speed of sound in air
dens=1.2;   % Density of medium (kg/m^3)

[NodeRhoZ,ElemNodeNum,M]=readnods(file,'y');

if nargin<7
   filename=[];
else
   filename=varargin{1};
end

for ik=1:length(fs)
   kk=2*pi*fs(ik)/c;
   disp(' ');
   disp(['  Frequency = ' num2str(fs(ik)) '  Number ' num2str(ik) ...
         ' of ' num2str(length(fs))]);
   [A,B,CConst]=BEMEquat0(NodeRhoZ,ElemNodeNum,kk,mter,[],filename);
   
   for mm=1:length(mter)
      
      m=mter(mm);   
      disp(['  Expansion term = ' num2str(m) '  Number ' num2str(mm) ...
            ' of ' num2str(length(mter))]);
   
      if isempty(filename)
         Am=A(:,:,mm);
         Bm=B(:,:,mm);
      else
         [Am,Bm]=getAB(mm,filename,M);
      end

      % 0 = A phi + B v + 4 pi phi^I
      % multiply B with i*k*dens*c to get p instead of phi
      Bm=i*kk*dens*c*Bm;
      
      % Compute the incoming field on the nodes
      pI=incoming(sources,NodeRhoZ,kk,m);

      %  Solve system of equations:
      pm(:,ik,mm)=(Am+Bm*diag(Y))\(-Bm*v(:,mm)-4*pi*pI);
      CondNum(ik,mm)=cond(Am);
   end
end

   
function  [Am,Bm]=getAB(im,filename,M);

Am=zeros(M,M);
Bm=zeros(M,M);                                             %*********************
for ii=1:M
   load(filename,['Arow' int2str(ii) 'm' int2str(im)]...        % load filename Arow$ii$m$im$ Brow$ii$m$im$
       ,['Brow' int2str(ii) 'm' int2str(im)]);
   eval(['Am(ii,:)=Arow' int2str(ii) 'm' int2str(im) ';']);     % Am(ii,:)=Arow$ii$m$im$;
   eval(['Bm(ii,:)=Brow' int2str(ii) 'm' int2str(im) ';']);     % Bm(ii,:)=Brow$ii$m$im$;    *********************
end
