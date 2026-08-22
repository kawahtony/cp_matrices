function [V,F] = read_off(filename)
%READ_OFF Read an OFF triangular mesh file.
%   [V,F] = READ_OFF(filename) returns vertices V (n-by-3) and faces F
%   (m-by-3) using MATLAB's 1-based indexing.

fid = fopen(filename,'r');
if fid < 0
    error('Could not open file: %s', filename);
end
cleanup = onCleanup(@() fclose(fid));

line = strtrim(fgetl(fid));
while startsWith(line, '#') || isempty(line)
    line = strtrim(fgetl(fid));
end
if ~strcmp(line, 'OFF')
    error('Not a valid OFF file.');
end

line = strtrim(fgetl(fid));
while startsWith(line, '#') || isempty(line)
    line = strtrim(fgetl(fid));
end
counts = sscanf(line, '%d %d %d');
nv = counts(1); nf = counts(2);

V = fscanf(fid, '%f %f %f', [3 nv])';
F = zeros(nf,3);
for k = 1:nf
    vals = fscanf(fid, '%d', 4);
    if vals(1) ~= 3
        error('This reader expects triangular faces only.');
    end
    F(k,:) = vals(2:4)' + 1; % OFF is 0-based; MATLAB is 1-based.
end
end


% function [vertex,face] = read_off(filename)
% 
% % read_off - read data from OFF file.
% %
% %   [vertex,face] = read_off(filename);
% %
% %   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
% %   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
% %
% %   Copyright (c) 2003 Gabriel Peyr
% 
% 
% fid = fopen(filename,'r');
% if( fid==-1 )
%     error('Can''t open the file.');
%     return;
% end
% 
% str = fgets(fid);   % -1 if eof
% if ~strcmp(str(1:3), 'OFF')
%     error('The file is not a valid OFF one.');    
% end
% 
% str = fgets(fid);
% [a,str] = strtok(str); nvert = str2num(a);
% [a,str] = strtok(str); nface = str2num(a);
% 
% 
% 
% [A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
% if cnt~=3*nvert
%     warning('Problem in reading vertices.');
% end
% A = reshape(A, 3, cnt/3);
% vertex = A;
% % read Face 1  1088 480 1022
% [A,cnt] = fscanf(fid,'%d %d %d %d\n', 4*nface);
% if cnt~=4*nface
%     warning('Problem in reading faces.');
% end
% A = reshape(A, 4, cnt/4);
% face = A(2:4,:)+1;
% 
% 
% fclose(fid);
