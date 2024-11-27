function ViscLyo(subname, varargin)   
    
    autoOn = 0;    
    if(nargin == 2)
        autoOn = 1;
        f1= figure('Renderer', 'painters', 'Position', [10 10 3000 1600]);    
    else
        f1 = figure(1);
    end
        
    
    addpath('./jsonlab/');
    addpath('./archive/');

    filename = strcat('../../', subname)
    moviename = strcat(filename, '/viscLyo.avi');
    figname = strcat(filename, '/LastFrame.jpg');
    dataname = strcat(filename, '/data.mat');
    done_message = "Done";
    
    ndefPQ_list = [];
    ndefNQ_list = [];
    loc_defPQlist = [];
    loc_defNQlist = [];
    loc_defPQ_arr = {};
    loc_defNQ_arr = {};
    
    
    time_list         = [];
    anchor_mean_list  = [];
    anchor_SD_list    = [];
    shape_list        = [];
    shape_list2       = [];
    
    if (autoOn == 1)
       wait_sim_start = 0;
       while (exist(filename, 'dir') ~= 7)           
           wait_sim_start = wait_sim_start+1;
           disp('Waiting')
           pause(60);
           if (wait_sim_start > 120) %wait for 120 minutes
               warning('Frame not found in 2 hours. Code stopped')
               return;
           end
       end       
    end
    ar = loadarchive(filename);


    writerObj = VideoWriter(moviename);
    writerObj.FrameRate = 4; %floor(getnframes(ar)/10);
    open(writerObj);   
    

    startIdx = ar.nstart/ar.ninfo;
    

    try
    for m=startIdx:getnframes(ar)
        m
        if (autoOn == 1)
            frame_name = strcat(filename, '/frame', int2str(ar.ninfo*(m+startIdx)), '.json');
            wait_it = 0;
            while (exist(frame_name, 'file') ~= 2)
               disp('Waiting')
               pause(60);
               wait_it = wait_it+1;
               if (wait_it > 60) %wait for 60 minutes
                   error('Frame not found in an hour.')                   
               end
            end
        end
        
        fr = loadframe(ar, m);
        fr = reshapeframe(fr);
        
        time_list = [time_list; ar.nstart + m*ar.ninfo];
        
        [ S, qx, qy ] = getdirector(fr);
        [ ux, uy ]    = getvelocity(fr);
        phi = fr.phi;
        
        siz = [fr.parameters.LY fr.parameters.LX];
        fr.QQxx2 = reshape(fr.QQxx2, siz);
        fr.QQyx2 = reshape(fr.QQyx2, siz); 
        [ S2, qx2, qy2 ] = getdirector(fr, 2);

        window_size = 3;
        averaging_filter = fspecial('average', [window_size window_size]);
        smoothedS2 = imfilter(S2, averaging_filter);

        clf; sp1=2; sp2=3;                              

        %picture of phi
        subplot(sp1,sp2,1); 
            [~, h] = contourf(phi, 40); colorbar; hold on;
            plotdirector(S, qx, qy, 5); hold on;
            set(h, 'LineStyle', 'none'); title('\phi');
            axis equal tight;  

        %velocity field
        subplot(sp1,sp2,2); 
            plotvelmagnitude(ux,uy); hold on;
            plotvelocity(ux,uy, 5); 
            hold off; title('velocity')
            axis equal tight;  
                  
        %ecm nematic field
        subplot(sp1, sp2, 3);                 
            plotorder(S2); hold on;
            
            plotdirector(S2, qx2, qy2, 4); hold on;
            clim([0, 1.5*mean(S2, "all")])
            [Dxpq2, Dypq2, Dxnq2, Dynq2, ~, ~] = plotdefects(ar.LX, ar.LY, qx2, qy2, (phi<= 0.70)&(smoothedS2>0.01)); hold on;
            %[Dxpq2, Dypq2, Dxnq2, Dynq2, ~, ~] = plotdefects(ar.LX, ar.LY, qx2, qy2, smoothedS2); hold on;
            title('ECM Nematic'); axis equal tight;
            try
                ndefPQ_list = [ndefPQ_list; length(Dxpq2)];
                ndefNQ_list = [ndefNQ_list; length(Dxnq2)];
                locPQ_temp = cat(2, Dxpq2, Dypq2);
                locNQ_temp = cat(2, Dxnq2, Dynq2);
                loc_defPQlist = [loc_defPQlist; locPQ_temp];
                loc_defNQlist = [loc_defNQlist; locNQ_temp];
                loc_defPQ_arr{m+1} = locPQ_temp;
                loc_defNQ_arr{m+1} = locNQ_temp;
            catch
                ndefPQ_list = [ndefPQ_list; 0];
                ndefNQ_list = [ndefNQ_list; 0];   
                loc_defPQ_arr{m+1} = [];
                loc_defNQ_arr{m+1} = [];
            end
            hold off
        
        %defect count
        subplot(sp1, sp2, 4);
        try
            plot(time_list, ndefPQ_list, '--o'); hold on;
            plot(time_list, ndefNQ_list, '--o'); hold off;
            legend('Plus defects', 'Minus defects');
            title('Number of defects'); xlabel('time')
        end
            
        %average anchoring
        subplot(sp1, sp2, 5);
            [dxPhi, dyPhi] = gradient(phi); 
            delPhi = sqrt(dxPhi.^2 + dyPhi.^2);
            dpQdp = (dxPhi.*dxPhi.*fr.QQxx2 - dyPhi.*dyPhi.*fr.QQxx2 + 2*dxPhi.*dyPhi.*fr.QQyx2)./(delPhi.^2.*S2);
            phi_interface_mask = (phi>=0.30) & (phi<=0.70); %delPhi >= 0.05;%
            anchor_ang_deg = 90-acosd(dpQdp(phi_interface_mask))/2;
            anchor_mean_list  = [anchor_mean_list; mean(anchor_ang_deg)];
            anchor_SD_list    = [anchor_SD_list; std(anchor_ang_deg)];
            errorbar(time_list, anchor_mean_list, anchor_SD_list, '--o'); 
            title('Average anchoring'); xlabel('time')
            
        
        subplot(sp1, sp2, 6);
            %shape index 1 = circularity
            total_area = sum(sum(phi));
            total_perimeter = length(anchor_ang_deg);
            shape_index = total_perimeter/sqrt(total_area)/5.65; %5.65 is the empirical shape index for a circle
            shape_list = [shape_list; shape_index];            
            
            yyaxis left
            plot(time_list, shape_list, '--o'); hold on;
            ylabel('Irregularity')
            
            %shape index 2 = eccentricity of bounding box (PBC issue)
            [xmesh, ymesh] = meshgrid(1:ar.LX, 1:ar.LY);
            phi_round = round(phi);
            SumP = sum(sum(phi_round));
            SumPX = sum(sum(phi_round.*xmesh));
            SumPX2 = sum(sum(phi_round.*xmesh.*xmesh));
            SumPY = sum(sum(phi_round.*ymesh));
            SumPXY = sum(sum(phi_round.*xmesh.*ymesh));
            SumPY2 = sum(sum(phi_round.*ymesh.*ymesh));
            Ixx = SumPX2/SumP - SumPX^2/SumP^2;
            Iyy = SumPY2/SumP - SumPY^2/SumP^2;
            Ixy = SumPXY/SumP - SumPY*SumPX/SumP^2;
            I_tensor = [Ixx, Ixy; Ixy, Iyy];
            lambdas = eig(I_tensor);
            shape_index2 = sqrt(max(lambdas)/min(lambdas));
            shape_list2 = [shape_list2; shape_index2];
            
            yyaxis right
            plot(time_list, shape_list2, '--o'); hold on;
            ylabel('Aspect Ratio')
       
        annotation('textbox', [0, 0.5, 0, 0], 'string', strcat('m=', int2str(m)))
        annotation('textbox', [0, 0.47, 0, 0], 'string', strcat('timestep=', int2str(m*ar.ninfo)))
        drawnow;

        frame = getframe(gcf);
        writeVideo(writerObj, frame);       


    end
    
    catch ME
        %close(writerObj);
        %save(dataname, 'ndefPQ_list', 'ndefNQ_list', 'anchor_mean_list', 'anchor_SD_list', 'shape_list', 'shape_list2', 'time_list', 'loc_defPQlist', 'loc_defNQlist', ...
        %'loc_defPQ_arr', 'loc_defNQ_arr');
         %saveas(f1, figname);
        rethrow(ME)
        % switch ME.identifier
        %     case 'MATLAB:audiovideo:VideoWriter:invalidDimensions'
        %         warning('Frame size changed. Video stopped')
        %         done_message= "Frame Size Change";
        %     case ''
        %         if strcmp(ME.message, 'input file does not exist')
        %             warning('Input file missing. Video stopped')
        %             done_message= "Missing Input";
        %         elseif strcmp(ME.message, 'Frame not found.')
        %             warning('Input frame missing. Video stopped')
        %         else
        %             warning(ME.message)
        %         end
        % 
        %     otherwise
        %         warning(ME.identifier)
        % end
    end

           
    annotation('textbox', [0.5, 0.42, 0.1, 0.1], 'string', done_message)
    close(writerObj);
    save(dataname, 'ndefPQ_list', 'ndefNQ_list', 'anchor_mean_list', 'anchor_SD_list', 'shape_list', 'shape_list2', 'time_list', 'loc_defPQlist', 'loc_defNQlist', ...
    'loc_defPQ_arr', 'loc_defNQ_arr');
     saveas(f1, figname);
    %%%load(dataname);   
    
    
end

