% Reshape all fields that need to be reshaped in a frame
% Inputs:
%   fr -- the frame to be reshaped
function answ = reshapeframe(fr)

    % size of domain
    siz = [fr.parameters.LY fr.parameters.LX];
    % reshape all fields
    fr.ff = reshape(fr.ff, fr.parameters.LY, fr.parameters.LX, []);
    fr.QQxx = reshape(fr.QQxx, siz);
    fr.QQyx = reshape(fr.QQyx, siz);
    
    try
        fr.ccA = reshape(fr.ccA, siz);        
    catch ME
        %repackage as a warning        
    end
    
    try
        fr.ccB = reshape(fr.ccB, siz);      
    catch ME
        %repackage as a warning        
    end
    
    try
       fr.phi = reshape(fr.phi, siz); 
    catch ME
        %repackage as a warning
    end
    
    answ = fr;