function [distP_array, distN_array, drop_size_array, drop_omega_array] = StatsToolbox3_ViscLyo(subname)
    addpath('./jsonlab/');
    addpath('./archive/');
    
    filename = strcat('../../', subname)
    %figname = strcat(filename, '/LastFrame.jpg');
    dataname = strcat(filename, '/ClusterDefectData.mat');
    
    %defect locations
    loidxP = 0; hiidxP = 0; loidxN = 0; hiidxN = 0;
    distP_list = []; distN_list = []; 
    clusterIrreg_P_list = []; clusterIrreg_N_list= []; %specially for one thing
    
    %droplet statistics
    drop_size_list = []; drop_omega_list = [];
    
    %cell array versions of the above
    distP_array={}; distN_array={};
    drop_size_array={}; drop_omega_array={};
    
    %

    try
        load(strcat('../../', subname, '/data.mat'), 'ndefPQ_list', ...
            'ndefNQ_list', 'loc_defPQlist', 'loc_defNQlist', 'loc_defPQ_arr', 'loc_defNQ_arr');
    catch
        try %old analysis
            load(strcat('../../', subname, '/data.mat'), 'ndefPQ_list', ...
            'ndefNQ_list', 'loc_defPQlist', 'loc_defNQlist');
        catch %doesn't work - redo
            disp('Missing data: Running ViscLyo');
            ViscLyo(subname);
            load(strcat('../../', subname, '/data.mat'), 'ndefPQ_list', ...
                'ndefNQ_list', 'loc_defPQlist', 'loc_defNQlist');
        end
    end
    
    ar = loadarchive(filename);
    
    NdefP_offset = getnframes(ar)-length(ndefPQ_list);
    NdefN_offset = getnframes(ar)-length(ndefNQ_list);

    try
    for m=1:getnframes(ar)
        m
        fr = loadframe(ar, m);
        fr = reshapeframe(fr);
        phi = fr.phi;        
                
        [X_grid, Y_grid] = meshgrid(1:ar.LX, 1:ar.LY);
        
        cluster_points = phi>=0.70;
        X_points = X_grid(cluster_points);
        Y_points = Y_grid(cluster_points);       
        data_points = cat(2, X_points, Y_points);       
        
        try
            %figure(2); clf;        
            [domainC, ~] = contour(phi, 2);
            [Cvals, Cpoly] = cleanContourMatrix(domainC);
            %use higher contours for droplets
            droplet_polys = Cpoly(Cvals>0.50);
            num_drops = length(droplet_polys);
            area_drops = zeros(num_drops, 1);
            omega_drops = zeros(num_drops, 1);
            for i = 1:num_drops
                area_drops(i) = polyarea(droplet_polys{i}(1, : ), droplet_polys{i}(2, : ));
                omega_drops(i) = perimeter(polyshape(droplet_polys{i}(1, : ), droplet_polys{i}(2, : )))/sqrt(area_drops(i));                
            end
            drop_size_list = [drop_size_list; area_drops];
            drop_omega_list = [drop_omega_list; omega_drops];
            drop_size_array{m} = area_drops;
            drop_omega_array{m} = omega_drops;
        catch
            disp('Cluster stats error in this timestep. Continuing...');
            drop_size_array{m} = [];
            drop_omega_array{m} = [];
        end

        try
            try
                query_pointsP = loc_defPQ_arr{m};
                query_pointsN = loc_defNQ_arr{m};
            catch
                %extract the defect values at the current timestep
                loidxP = hiidxP+1;
                hiidxP = loidxP + ndefPQ_list(max(m-NdefP_offset, 1))-1;
                loidxN = hiidxN+1;
                hiidxN = loidxN + ndefNQ_list(max(m-NdefN_offset, 1))-1;
                
                query_pointsP = loc_defPQlist(loidxP:hiidxP, :);
                query_pointsN = loc_defNQlist(loidxP:hiidxP, :);
            end
            
            try
                [nearP, distP] = dsearchn(data_points, query_pointsP);
                [nearN, distN] = dsearchn(data_points, query_pointsN);        
                distP_list = [distP_list; distP];
                distN_list = [distN_list; distN];
                distP_array{m} = distP;
                distN_array{m} = distN;
            catch
                distP = [];
                distN = [];
            end
        catch
            disp('Defect stats error. Continuing');
            distP_array{m} = [];
            distN_array{m} = [];
        end
        
            %figure(1); clf;  
            % histogram(distP, 'Normalization', 'probability'); hold on;
            % histogram(distN, 'Normalization', 'probability'); 
            % legend('+ defs', '- defs');         
            % drawnow;
            % %Specifically to calculate defect-irregularity data - disable
            % %comment otherwise for speed
            % for defidx = 1:length(distP)
            %     defloc = data_points(nearP(defidx), :);
            %     write_flag = 0;
            %     for i=1:num_drops
            %         inpoly_flag = inpolygon(defloc(1), defloc(2), droplet_polys{i}(1, : ), droplet_polys{i}(2, : ));
            %         if (inpoly_flag == 1)
            %             clusterIrreg_P_list = [clusterIrreg_P_list; omega_drops(i)];
            %             write_flag = 1;
            %             break;
            %         end
            %     end
            %     if (write_flag == 0) %cluster nor found
            %         clusterIrreg_P_list = [clusterIrreg_P_list; mean(omega_drops)];
            %     end
            % end
            % 
            % for defidx = 1:length(distN)
            %     defloc = data_points(nearN(defidx), :);
            %     write_flag = 0;
            %     for i=1:num_drops
            %         inpoly_flag = inpolygon(defloc(1), defloc(2), droplet_polys{i}(1, : ), droplet_polys{i}(2, : ));
            %         if (inpoly_flag == 1)
            %             clusterIrreg_N_list = [clusterIrreg_N_list; omega_drops(i)];
            %             write_flag = 1;
            %             break;
            %         end
            %     end
            %     if (write_flag == 0) %cluster nor found - shouldn't happen
            %         clusterIrreg_N_list = [clusterIrreg_N_list; mean(omega_drops)];
            %     end
            % end                 

    end
    catch ME
        if (ME.message == "input file does not exist")
            disp('Frame not found. Exiting.')
        else
            disp('Outer Loop error')
            rethrow(ME);
        end
    end

    try
        distP_cut = distP_list(ceil(0.50*length(distP_list)):end);
    catch
        distP_cut = [];
    end
    try
        distN_cut = distN_list(ceil(0.50*length(distN_list)):end);
    catch
        distN_cut = [];
    end
    
    try
        figure(1); clf;
        subplot(1, 2, 1);
            hP = histogram(distP_cut, 'Normalization', 'probability'); hold on;
            hN = histogram(distN_cut, 'Normalization', 'probability'); hold off;
        
        
        binvalsP = movmean(hP.BinEdges, 2);
        binvalsP = binvalsP(2:end);
        pdf_lateP= griddedInterpolant(binvalsP, hP.Values);
        
        binvalsN = movmean(hN.BinEdges, 2);
        binvalsN = binvalsN(2:end);
        pdf_lateN= griddedInterpolant(binvalsN, hN.Values);
        
        clf;
        xrange = 0.6:0.2:60;
        plot(xrange, pdf_lateP(xrange)); hold on;
        plot(xrange, pdf_lateN(xrange)); hold on;
        legend('+ defs', '- defs');        
    catch
        %do nothing
        pdf_lateN = [];
        pdf_lateP = [];

    end

    subplot(1, 2, 2);
        histogram(drop_size_list, 'Normalization', 'probability'); hold on;
        xlabel('size(s)')
        %histogram(drop_omega_list, 'Normalization', 'probability'); hold off;
    
    save(dataname, 'distP_list','distN_list', 'pdf_lateP', 'pdf_lateN', 'drop_size_list', 'drop_omega_list', ...
        'distP_array', 'distN_array', 'drop_size_array', 'drop_omega_array', 'clusterIrreg_P_list' , 'clusterIrreg_N_list');
    %saveas(gcf, figname);
end

