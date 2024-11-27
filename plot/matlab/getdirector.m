% Get director from frame (in two dimensions) and reshape. Note, dx>0
% always
function [S, dx, dy] = getdirector(fr, varargin)
    
    if (nargin == 1)
        QQxx = fr.QQxx;
        QQyx = fr.QQyx;    
    elseif (varargin{1} == 1)
        QQxx = fr.QQxx1;
        QQyx = fr.QQyx1;
    elseif (varargin{1} == 2)
        QQxx = fr.QQxx2;
        QQyx = fr.QQyx2;
    else
        QQyx = 0;
        QQxx = 0;
    end
    
    S  = sqrt(QQxx.^2 + QQyx.^2);              %Eigenvalues          
    dx = sqrt((ones(size(S)) + QQxx./S)/2);       %Normalized eigenvector
    dy = sign(QQyx).*sqrt((ones(size(S)) - QQxx./S)/2);