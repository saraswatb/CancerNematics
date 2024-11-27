function [Cvals, Cpoly] = cleanContourMatrix(qh)
    %This function reads out a contour matrix qh into easier form
    %Cval: Contour Val, Cpoly: Contour polygon indices array
    Cvals = [];
    Cpoly = {};   
    sizeQH = size(qh);

    i=1; cnt=1;
    while(i<sizeQH(2))
       Cvals = [Cvals; qh(1, i)];
       Ncontour = qh(2, i);
       Cpoly{cnt} = qh(:, i+1:i+Ncontour);
       cnt = cnt+1;
       i = i+Ncontour+1;
    end
end