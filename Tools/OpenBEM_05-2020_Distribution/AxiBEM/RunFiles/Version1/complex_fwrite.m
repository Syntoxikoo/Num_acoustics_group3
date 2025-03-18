function complex_fwrite(fid,A)

% function complex_fwrite(fid,A)
% Write an array of complex numbers to a binary file
% Input:
%   -fid : File id (e.g. from fopen)
%   -A   : The matrix to write to the binary file (2D-matrix)

tmp(1:2:2*size(A,1)-1,1:size(A,2))=real(A);
tmp(2:2:2*size(A,1),1:size(A,2))=imag(A);
fwrite(fid,tmp,'double');
