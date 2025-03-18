function [A]=complex_fread(fid,size_A)

% function [A]=complex_fread(fid,size_A)
% Read a complex matrix from a binary file
% Input:
%   -fid    : File id (e.g. from fopen)
%   -size_A : row vector of integers
%             1st column is number of rows in the output matrix
%             2nd column is number of columns in output matrix
% Output:
%   -A      : The matrix as read columnwise from the binary file

tmp=fread(fid,[2*size_A(1),size_A(2)],'double');
A=tmp(1:2:2*size_A(1)-1,:)+i*tmp(2:2:2*size_A(1),:);
