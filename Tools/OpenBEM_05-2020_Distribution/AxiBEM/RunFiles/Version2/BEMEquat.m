function [A,B,C]=BEMEquat(rzb,Topology,k,m,varargin)

%  Calculate C constants and A,B matrices for Helmholtz equation:
%
%  [A,B,C]=BEMEquat(rzb,Topology,k,m{,chiefpoints})
%  [A,B,C]=BEMEquat(rzb,Topology,k,m,chiefpoints,filename)
%
%  Axisymetric BEM calculation.
%  Scattering of an axisymmetric incident field.
%  Radiation in testfase
%  
%  Input parameters (in SI units):
%  
%    -rzb     : Geometry matrix
%               one row for each node
%               rho-coordinate in column 1
%               z-coordinate in column 2
%               body number in column 3 (not necessary if there is
%               just one body)
%    -Topology: Topology matrix
%               One row for each element
%               Each row contains the global node number
%               in the column no. corresponding to the local
%               node number. Last column contains the body number.
%               (If left empty it is generated automatically)
%    -k       : Wave number (k=2*pi*f/c) (real scalar)
%    -m       : Circumferential mode numbers (vector of integers>=0)
%    -chiefpoints: Like 'rzb', but contains CHIEF points instead
%               One row for each chief point
%    -filename : file name to be used for storing the BEM coefficient
%               matrices on the disk.
%               The variables in the file are named Am$ and Bm$,
%               where $ is the index in the m vector.
%               The file also contains the variables 'm' and 'k'. 
%
%  Output parameters:
%
%    -A :       In case no filename is given, A contains the A
%               coefficient matrix, the indexes are:
%               A(M nodes,M nodes,m terms)
%               If a file is used, it returns an empty matrix.
%    -B :       In case no filename is given, B contains the B
%               coefficient matrix, the indexes are:
%               A(M nodes,N elements,nodes in element,m terms)
%               If a file is used, it returns an empty matrix.
%    -C :       real column vector containing the C-constants

M=size(rzb,1);
[nel, ncols]=size(Topology);
nknel=ncols-1;
if nargin<5
   chiefpoints=[];
else
   chiefpoints=varargin{1};
end
nchiefp=size(chiefpoints,1);

if nargin<6
   filename=[];
else
   filename=varargin{2};
end

% Decide whether to save matrices in files or memory
if ~isempty(filename)
   save(filename,'k','m');
end

% Find 'A' and 'B' matrices
for ii=1:M+nchiefp
  if ii <= M
     przb=rzb(ii,:);
  else
     przb=chiefpoints(ii-M,:);
  end
  cii=4*pi*(1+sign(przb(3)))/2; % 4*pi or 0 for exteror/interior domain.(VC)
  Am=zeros(M,length(m));
  Bm=zeros(nel,nknel,length(m));
  for iel=1:nel
    elknrzb(1:nknel,:)=rzb(Topology(iel,1:nknel),:);
    % Singular part
    [g,h,cjj] = intF2(przb,elknrzb,m);
    cii=cii+cjj;
    Am(Topology(iel,1:nknel),:)=Am(Topology(iel,1:nknel),:)+h(1:nknel,:);
    Bm(iel,1:nknel,:)=g(1:nknel,:);
    % Oscillating part
    [g,h] = intF1(przb,elknrzb,k,m);
    Am(Topology(iel,1:nknel),:)=Am(Topology(iel,1:nknel),:)+h(1:nknel,:);
    Bm(iel,1:nknel,:)=Bm(iel,1:nknel,:)+g(1:nknel,:);
  end

  disp(sprintf('Row %g/%g C constant = %1.15e',ii,M+nchiefp,cii/4/pi));
  if ii <= M
     Am(ii,:)=Am(ii,:)-cii;
     C(ii)=cii;
  end
  if ~isempty(filename)  % Store the row on disk for all m's
     for im=1:length(m)
        eval(['Arow' int2str(ii) 'm' int2str(im) '=Am(:,im).'';']);   % Arow$ii$m$im$=Am(:,im).';
        eval(['Brow' int2str(ii) 'm' int2str(im) '=Bm(:,:,im);']);    % Brow$ii$m$im$=Bm(:,:,im);
        save(filename,['Arow' int2str(ii) 'm' int2str(im)]...         % save filename Arow$ii$m$im$ Brow$ii$m$im$ -append
           ,['Brow' int2str(ii) 'm' int2str(im)],'-append');
     end
  else  % Store in memory
     A(ii,:,:)=Am(:,:);
     B(ii,:,:,:)=Bm(:,:,:);
  end
end

% Arrange data in the file as arrays for each m,
% never loading too much at a time
if ~isempty(filename)
%   for im=1:length(m)
%      eval(['Am' int2str(im) '=zeros(M+nchiefp,M);']);            % Am$im$=zeros(M+nchiefp,M);
%      eval(['Bm' int2str(im) '=zeros(M+nchiefp,nel,nkel);']);     % Bm$im$=zeros(M+nchiefp,nel,nkel);
%      for ii=1:M+nchiefp
%         load(filename,['Arow' int2str(ii) 'm' int2str(im)]...       % load filename Arow$ii$m$im$ Brow$ii$m$im$
%            ,['Brow' int2str(ii) 'm' int2str(im)]);
%         eval(['Am' int2str(im) '(ii,:)=Arow' ...                    % Am$im$(ii,:)=Arow$ii$m$im$;
%               int2str(ii) 'm' int2str(im) ';']);
%         eval(['Bm' int2str(im) '(ii,:,:)=Brow' ...                  % Bm$im$(ii,:,:)=Brow$ii$m$im$;
%               int2str(ii) 'm' int2str(im) ';']);
%         eval(['Arow' int2str(ii) 'm' int2str(im) '=[];']);          % Arow$ii$m$im$=[];
%         eval(['Brow' int2str(ii) 'm' int2str(im) '=[];']);          % Brow$ii$m$im$=[];
%         save(filename,['Arow' int2str(ii) 'm' int2str(im)]...       % save filename Arow$ii$m$im$ Brow$ii$m$im$ -append
%            ,['Brow' int2str(ii) 'm' int2str(im)],'-append');
%      end
%      save(filename,['Am' int2str(im)],...                        % save filename Am$im$ Bm$im$ -append
%         ['Bm' int2str(im)],'-append'); 
%      eval(['clear Am' int2str(im) ';']);                         % clear Am$im$;
%      eval(['clear Bm' int2str(im) ';']);                         % clear Bm$im$;
%   end
   A=[];B=[];
end
