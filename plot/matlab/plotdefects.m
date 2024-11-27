function [xxpq, yypq, xxnq, yynq, psi, psi_n] = plotdefects(nx, ny, dx, dy,S, varargin)
    qc=calcs(dx, dy);           
    
    if (nargin == 5)
        [xxpq, yypq, xxnq, yynq, psi, psi_n] = chargearray(dx,dy,S,nx ,ny, qc);
    else
        [xxpq, yypq, xxnq, yynq, psi, psi_n] = chargearray(dx,dy,S,nx ,ny, qc, 0);
    end
end

function [xxpq, yypq, xxnq, yynq, psi, psi_n] = chargearray(dx,dy,S,nx ,ny, qc, varargin)
    plotOn = 1;
    
    
    yypq = [];
    xxpq = [];
    yynq = [];
    xxnq = [];
 
    if (nargin >= 7)
        
        plotOn = 0;
    end
    
    for i=2:ny-1; 
        for j=2:nx-1; 
            yy=0; xx=0; n1=0;
            %if(abs(qc(i,j))>0.4&& S(i,j)<0.5)
            if(abs(qc(i,j))>0.4 && S(i,j)>0.1) % changed defect detection thresholds
                ql = sign(qc(i,j)); yy = i; xx = j; n1 = 1; qc(i,j) = 0;
                for ii=-1:1
                    for jj=-1:1
                        if(ql*qc(i+ii,j+jj)>0.4)
                            yy = yy+i+ii; xx = xx+j+jj; n1 = n1+1;
                            qc(i+ii,j+jj) = 0;
                        end
                    end
                end
                if(ql>0)
                    yypq = [yypq; yy/n1]; xxpq = [xxpq; xx/n1];
                else
                    yynq = [yynq; yy/n1]; xxnq = [xxnq; xx/n1];
                end
            end
        end
    end

    if (plotOn == 1)
        plot(xxpq, yypq, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); hold on;
        plot(xxnq, yynq, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    end
    psi=[];

    for P=1:size(xxpq)
        k=1/2;
        Qxx=S.*(dx.*dx-1/2);
        Qyy=S.*(dy.*dy-1/2);
        Qxy=S.*(dx.*dy);
        [dxQxx,dyQxx]=gradient((Qxx));
        [dxQyy,dyQyy]=gradient((Qyy));
        [dxQxy,dyQxy]=gradient((Qxy));
        [dxxQxx,dyxQxx]=gradient((dxQxx));
        [dxyQyy,dyyQyy]=gradient((dyQyy));
        [dxyQxy,dyyQxy]=gradient((dyQxy));    
        x1=floor(xxpq(P));
        x2=floor(xxpq(P));
        x3=ceil(xxpq(P));
        x4=ceil(xxpq(P));
        y1=floor(yypq(P));
        y2=ceil(yypq(P));
        y3=ceil(yypq(P));
        y4=floor(yypq(P));

        [dxQxx,dyQxx]=gradient((Qxx));
        [dxQyy,dyQyy]=gradient((Qyy));
        [dxQxy,dyQxy]=gradient((Qxy));
        dem=(dxQxy(y1,x1)-dyQxx(y1,x1))+(dxQxy(y2,x2)-dyQxx(y2,x2))+(dxQxy(y3,x3)-dyQxx(y3,x3))+(dxQxy(y4,x4)-dyQxx(y4,x4));
        num=(dxQxx(y1,x1)+dyQxy(y1,x1))+(dxQxx(y2,x2)+dyQxy(y2,x2))+(dxQxx(y3,x3)+dyQxy(y3,x3))+(dxQxx(y4,x4)+dyQxy(y4,x4));
        psi=[psi; k/(1-k)*atan2(dem,num)];

    end
    psi_n=[];

    for P=1:size(xxnq)
        k=-1/2;
        Qxx=S.*(dx.*dx-1/2);
        Qyy=S.*(dy.*dy-1/2);
        Qxy=S.*(dx.*dy);
        [dxQxx,dyQxx]=gradient((Qxx));
        [dxQyy,dyQyy]=gradient((Qyy));
        [dxQxy,dyQxy]=gradient((Qxy));
        [dxxQxx,dyxQxx]=gradient((dxQxx));
        [dxyQyy,dyyQyy]=gradient((dyQyy));
        [dxyQxy,dyyQxy]=gradient((dyQxy));    
        x1=floor(xxnq(P));
        x2=floor(xxnq(P));
        x3=ceil(xxnq(P));
        x4=ceil(xxnq(P));
        y1=floor(yynq(P));
        y2=ceil(yynq(P));
        y3=ceil(yynq(P));
        y4=floor(yynq(P));
  
        [dxQxx,dyQxx]=gradient((Qxx));
        [dxQyy,dyQyy]=gradient((Qyy));
        [dxQxy,dyQxy]=gradient((Qxy));

    
        dem=(-dxQxy(y1,x1)-dyQxx(y1,x1))+(-dxQxy(y2,x2)-dyQxx(y2,x2))+(-dxQxy(y3,x3)-dyQxx(y3,x3))+(-dxQxy(y4,x4)-dyQxx(y4,x4));
        num=(dxQxx(y1,x1)-dyQxy(y1,x1))+(dxQxx(y2,x2)-dyQxy(y2,x2))+(dxQxx(y3,x3)-dyQxy(y3,x3))+(dxQxx(y4,x4)-dyQxy(y4,x4));
        psi_n=[psi_n; k/(1-k)*atan2(dem,num)];
 
    end
    if (plotOn == 1)
        quiver_size=5;
        hold on
        dq=quiver(xxpq,yypq,cos(psi)*quiver_size,sin(psi)*quiver_size,0);
        set(dq,'LineWidth',3,'Color', 'r');
        set(dq, 'ShowArrowHead', 'off');
        dq=quiver(xxnq,yynq,cos(psi_n)*quiver_size,sin(psi_n)*quiver_size,0);
        set(dq,'LineWidth',3,'Color','g');
        set(dq, 'ShowArrowHead', 'off');
        dq=quiver(xxnq,yynq,cos(psi_n+pi*2/3)*quiver_size,sin(psi_n+pi*2/3)*quiver_size,0);
        set(dq,'LineWidth',3,'Color','g');
        set(dq, 'ShowArrowHead', 'off');
        dq=quiver(xxnq,yynq,cos(psi_n+pi*4/3)*quiver_size,sin(psi_n+pi*4/3)*quiver_size,0);
        set(dq,'LineWidth',3,'Color','g');
        set(dq, 'ShowArrowHead', 'off');
    end
    
end
  
 
    
    
