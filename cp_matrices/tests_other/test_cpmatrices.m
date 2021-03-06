% test_cpmatrices.m
%
% I made this code just to check what the matrices generated by the various
% cp_matrices functions look like.
%
% The band is a small, square domain, inside a larger square domain. So the
% matrices should be arranged as usual, i.e. not rearranged according to
% the indexing of some cp band.
%
% If the cp_matrices code could handle the band being the whole domain,
% this smaller square wouldn't be needed.
%
% input the test matrix at the bottom of the script


dx = 1;

% the big square
% change the grid spacing to get an idea of how the matrix is formed

x1d = (0:dx:100)';
y1d = x1d;

[xx,yy] = meshgrid(x1d,y1d);

% the little square
% four-points-by-four-points in the middle of the bigger square

iMiddle = round(length(x1d)/2);

min = x1d(iMiddle-2);
max = x1d(iMiddle+2);

band = find( (xx<max).*(xx>min) .*(yy<max).*(yy>min) );

%% input test matrix here

L = laplacian_2d_matrix(x1d,y1d,2,band);

% print results:

full(L)
