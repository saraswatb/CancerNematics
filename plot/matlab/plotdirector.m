function qh = plotdirector(S, dx, dy, nb)
  
    
    [x, y] = meshgrid(1:nb:size(S,1), 1:nb:size(S,2));
    %size(x)
    dx = dx(1:nb:end,1:nb:end); dy = dy(1:nb:end,1:nb:end); S = S(1:nb:end,1:nb:end); 
%     size(S)
%     size(dx)
%     size(dy)
   
    %qh = quiver(x'-dx/2, y'-dy/2, (S.*dx), (S.*dy), 'ShowArrowHead', 'Off');
    qh = quiver(y'+dy/2, x'+dx/2, (S.*dx), (S.*dy), 'ShowArrowHead', 'Off');
    %set(qh, 'Color', [37.0/255, 123.0/255, 124.0/255]); 
    set(qh, 'Color', 'k'); 
    
    