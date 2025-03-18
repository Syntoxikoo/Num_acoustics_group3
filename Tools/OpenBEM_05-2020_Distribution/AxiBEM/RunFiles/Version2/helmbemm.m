function [phi,condnums]=helmbemm(rzb,ElemNodeNum,chiefpoints,k,m,beta,dphidn,varargin)

%  [phi,condnums]=flobem(rzb,topology,chiefpoints,k,m,beta,dphidn{,pI,jobname})
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
%    -topology: Topology matrix
%               One row for each element
%               Each row contains the global node number
%               in the column no. corresponding to the local
%               node number. Last column contains the body number.
%               (If left empty it is generated automatically)
%    -chiefpoints: Like 'rzb', but contains CHIEF points instead
%               One row for each chief point
%    -k       : Wavenumber (k=2*pi*f/c) (real scalar)
%    -m       : Circumferential mode numbers (vector of integers>=0)
%    -beta    : beta=j*k*Y,         where Y is the locally reacting
%               admitance (non-dimensionalised wrt 1/(density*c)
%               Complex matrix with one row for each element
%               each column contains the value at the corresponding
%               local node number.
%               A completely empty matrix can be used if the entire
%               body is hard.
%    -dphidn  : Derivative of phi wrt. the normal to the surface
%               Cell array with one cell for each m-value.
%               Each cell contains a complex matrix with the row
%               number corresponding to element number and the
%               column number corresponds to the local node number.
%               A completely empty matrix can be used if the entire
%               body has no excitation.
%    -jobname : Jobname to be used for storing the BEM coefficient
%               Matrices on the disk.
%               The file name in which the A-matrix is stored is
%               Obtained by concatenating jobname,'Am_',m
%               The file name in which the B-matrix is stored is
%               Obtained by concatenating jobname,'Bm_',m
%    -pI      : Incoming field (acoustic velocity potential)
%               Cell array with one cell for each m-value.
%               Each cell contains a complex column vector with one row
%               for each node followed by one row for each chief point
%
%  Output parameters:
%  
%    -phi     : For each m-value a complex column vector containing
%               the solution along the generator is returned. Each row
%               contains the solution at the corresponding node number.
%               If there is more than one m-value, a cell array of
%               column vectors is returned.
%    -condnums: Cell array with one cell for each m-value
%               Each cell contains a real column vector containing the
%               condition numbers of the coefficient matrix for each
%               chief point that is added. First row is without CHIEF.
%               NOTE: Omitting this argument saves a lot of CPU-time

%  Key points about input
[M, ncolrzb] = size(rzb);
if ncolrzb<3
	rzb=[rzb ones(M,1)];
end
NumBodies = rzb(M,3);

nchiefp = size(chiefpoints,1);
if size(chiefpoints,2)<3
	chiefpoints=[chiefpoints ones(nchiefp,1)];
end

if nargin<8
   pI=[];
else
   pI=varargin{1};
end
if isempty(pI)
   pI=zeros(M+nchiefp,1);
end

if nargin<9
   jobname='helmbemm';
else
   jobname=varargin{2};
end

if isempty(ElemNodeNum)
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
        if inode>M
           ElemNodeNum(iel,3)=1;
           break;
        end
     end
   end
end
nel=size(ElemNodeNum,1);
nknel=size(ElemNodeNum,2)-1;
Topology=ElemNodeNum;

% Admitance
if isempty(beta)
   beta=zeros(nel,nknel);
end

% Excitation
if isempty(dphidn)
   dphidn=zeros(nel,nknel);
end

for ik=1:length(k)
   if length(k)>1
      disp(sprintf('Frequency %3d of %4d',ik,length(k)));
   end
   kk=k(ik);

   %  Calculate 'A' and 'B' for all m's (and store on disk)
   [Afiles,Bfiles]=BEMEquat(rzb,ElemNodeNum,kk,m,chiefpoints,jobname);

   % Solve for one 'm' at a time
   for im=1:length(m)
      disp(sprintf('Solving for m=%d',m(im)));

      % load A-matrix from disk
      Afid=fopen(sprintf('%sAm_%d',jobname,im),'r');
      A=complex_fread(Afid,[M M+nchiefp]);
      fclose(Afid);
      A=A.';

      % Calculate right hand side
      if iscell(dphidn)
         dphidnkm=dphidn{ik,im};
      else
         dphidnkm=dphidn;
      end
      if iscell(pI)
         pIkm=pI{ik,im};
      else
         pIkm=pI;
      end
      Bfid=fopen(sprintf('%sBm_%d',jobname,im),'r');
      for ii=1:M+nchiefp
         Bii=complex_fread(Bfid,[nel nknel]);
         rhs(ii,1)=sum(sum(Bii.*dphidnkm));
         for ikn=1:nknel
            A(ii,Topology(:,ikn))=A(ii,Topology(:,ikn)) ...
                                  -(Bii(:,ikn).*beta(:,ikn)).';
         end
      end
      fclose(Bfid);
      rhs=rhs-4*pi*pIkm;

      % Test if A-matrix is singular
      if rank(A,norm(A)*1e-3)<min(size(A))
         warning(sprintf('Coefficient matrix only has practical rank %g', ...
                         rank(A,norm(A)*1e-3)  ));
      end

      %  Solve system of equations:
      phi{ik,im}=A\rhs;

      % Test that the solution works
      if norm(rhs-A*phi{ik,im}) > 1e-2*norm(rhs)
         warning(sprintf('Error (2-norm) of rhs= %g / %g', ...
         norm(rhs-A*phi{ik,im}), norm(rhs)  ));
      end

      if nargout>=2  % Calculate condition numbers if desired
         for ichief=0:nchiefp
            condnums{ik,im}(ichief+1,1)=cond(A(1:M+ichief,:));
         end
      end

      if nargin<9  % Delete files if no jobname was specified
         delete(sprintf('%sAm_%d',jobname,im));
         delete(sprintf('%sBm_%d',jobname,im));
      end
   end
end

% Don't confuse people with cells if they aren't needed
if length(k)==1 & length(m)==1
   phi=phi{1,1};
end




