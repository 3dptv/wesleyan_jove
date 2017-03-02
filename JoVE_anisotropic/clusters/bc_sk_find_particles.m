function bc_sk_find_particles(cpvs,cpv_index_good,cls_min,cls_max,tr_min,tr_max,threshold,savepath,filename)
        
    ncams = length(cpvs);
    tr_flag = zeros(1,ncams);
    cls_flag = zeros(1,ncams);
    
    for icam=1:ncams      
            cpv=cpvs(icam);
            frewind(cpv);
            [cpvheader,rcnt] = fread(cpv,20,'*uint8');
            if rcnt ~= 20
                warning('Error reading cpv file header');
            end
            xpix = uint32(cpvheader(5))+uint32(cpvheader(6))*256;
            ypix = uint32(cpvheader(7))+uint32(cpvheader(8))*256;
            
            nframes = size(cpv_index_good(icam).data,1);
            currentFrameNum = zeros(nframes,1);
            
            % create mask & gutter
            xi = [50 50 1230 1230];
            yi = [974 100 100 974];
            mask=poly2mask(xi,yi,double(ypix),double(xpix));
            mask = uint8(mask);

            gutter = find(mask==0);
            
            ntracers = 1;
            nclusters = 1; 
            framecnt = 1;
            frametime = 0;
            
            for iframe=1:nframes
                if mod(iframe,5)==0
                display(sprintf('working on frame: %d of %d camera: %d, file: %s ', iframe, nframes, icam, filename));
                end
                frame = tic;
                fseek(cpv,cpv_index_good(icam).data(iframe,1),'bof');
                currentFrameNum(iframe) = fread(cpv,1,'*uint32');
                [p, rcnt] = fread(cpv,4,'uint8=>uint32');
                pixcnt = bitshift(p(4),14)+bitshift(p(3),6)+bitshift(p(2),-2);
                pixlist = zeros(pixcnt,3,'uint32');
                dat = fread(cpv,[4,pixcnt],'uint8=>uint32');
                last3bitarray=uint32(07)+zeros(1,pixcnt,'uint32');
                pixlist(:,3)= dat(1,:);
                %pixel locations range [0,xpix-1] by [0,ypix -1]
                pixlist(:,1) = dat(2,:)+bitshift(bitand(dat(3,:),last3bitarray),8); % decoding x position
                pixlist(:,2) = bitshift(dat(3,:),-3) + bitshift(dat(4,:),5);% decoding y position

                img = zeros(ypix,xpix,'uint8');
                for i=1:pixcnt
                    if(pixlist(i,1) < xpix && pixlist(i,1) > 0 && pixlist(i,2) < ypix && pixlist(i,2) > 0 && pixlist(i,3) > threshold)
                        img(pixlist(i,2)+1,pixlist(i,1)+1)=pixlist(i,3);
                    end
                end
                
                imgb=im2bw(img,0.01); % convert image to binary
                pix_info=regionprops(imgb,img,'Centroid','WeightedCentroid','Area','PixelIdxList','MaxIntensity','MeanIntensity');
                center=cat(1,pix_info.WeightedCentroid);
                area=cat(1,pix_info.Area);
                maxbright=cat(1,pix_info.MaxIntensity);
                meanbright=cat(1,pix_info.MeanIntensity);

                % find tracers
                ind = find(area>=tr_min(icam) & area<=tr_max(icam));
                if ~isempty(ind)
                    tracers(ntracers:ntracers+size(ind,1)-1,1:2) = center(ind,:); % x and y position
                    tracers(ntracers:ntracers+size(ind,1)-1,3) = maxbright(ind); % brightness
                    tracers(ntracers:ntracers+size(ind,1)-1,4) = cpv_index_good(icam).data(iframe,1); % pos in cpv file
                    tracers(ntracers:ntracers+size(ind,1)-1,5) = size(ind,1); % number of pixels in frame
                    tracers(ntracers:ntracers+size(ind,1)-1,6) = currentFrameNum(iframe); % frame number
                    ntracers=ntracers + size(ind,1);
                end
                
                % find clusters
                center=cat(1,pix_info.Centroid);
                ind = find(area>cls_min(icam) & area<cls_max(icam));
                if ~isempty(ind)
                    ncls=length(ind);
                    scls=ind;
                    for idz=length(ind):-1:1
                        if ~isempty(intersect(pix_info(ind(idz)).PixelIdxList,gutter));
                            ind(idz)=[];
                            scls(idz)=[];
                            ncls=ncls-1;
                        else
                            scls(idz)=size(pix_info(ind(idz)).PixelIdxList,1);
                        end
                    end
                end
                if ~isempty(ind)
                    clusters(nclusters:nclusters+ncls-1,1:2)=center(ind,:); % x,y center position of particle
                    clusters(nclusters:nclusters+ncls-1,3)=meanbright(ind);
%                     clusters(nclusters:nclusters+ncls-1,4)=area(ind); % size of particles
%                     clusters(nclusters:nclusters+ncls-1,4)=scls; % size of particles
                    clusters(nclusters:nclusters+ncls-1,4)=cpv_index_good(icam).data(iframe,1);
                    clusters(nclusters:nclusters+ncls-1,5)=ncls; % number of particles
                    clusters(nclusters:nclusters+ncls-1,6)=currentFrameNum(iframe); % frame number
                    nclusters = nclusters+ncls;
                end
                if mod(iframe,5)==0
                display(sprintf('Time for current frame: %f sec', toc(frame)));
                frametime=frametime+toc(frame);
                display(sprintf('Average time/frame: %f sec\nExpected to finish: %f min\n', frametime/framecnt,((frametime/framecnt)*nframes-framecnt*(frametime/framecnt))/60));
                end
                framecnt=framecnt+1;
            end
            display(sprintf('averaged %f clusters/frame on camera %d',nclusters/nframes,icam))
            display(sprintf('averaged %f tracers/frame on camera %d \n',ntracers/nframes,icam))
            if ntracers~=1;
                write2_gdf(tracers,[savepath filename '_cam' num2str(icam) '_tracers.gdf']);
                tr_flag(icam)=0;
            else
                warning('no tracers found!');
                tr_flag(icam)=1;
            end
            if nclusters~=1;
                write2_gdf(clusters,[savepath filename '_cam' num2str(icam) '_clusters.gdf']);
                cls_flag(icam)=0;
            else
                warning('no clusters found!');
                cls_flag(icam)=1;
            end
            tracers=[];
            clusters=[];
    end
    
    % remove frames that don't have particles on all cameras
    if ~max(tr_flag);
        tr1=read2_gdf([savepath filename '_cam1_tracers.gdf']);
        tr2=read2_gdf([savepath filename '_cam2_tracers.gdf']);
        tr3=read2_gdf([savepath filename '_cam3_tracers.gdf']);
        tr4=read2_gdf([savepath filename '_cam4_tracers.gdf']);

        if (~isempty(tr1) && ~isempty(tr2) && ~isempty(tr3) && ~isempty(tr4) )
            tr1(~ismember(tr1(:,6),tr2(:,6)),:)=[];
            tr1(~ismember(tr1(:,6),tr3(:,6)),:)=[];
            tr1(~ismember(tr1(:,6),tr4(:,6)),:)=[];
            tr2(~ismember(tr2(:,6),tr1(:,6)),:)=[];
            tr3(~ismember(tr3(:,6),tr1(:,6)),:)=[];
            tr4(~ismember(tr4(:,6),tr1(:,6)),:)=[];
            write2_gdf(tr1,[savepath filename '_cam1_tracers_final.gdf']);
            write2_gdf(tr2,[savepath filename '_cam2_tracers_final.gdf']);
            write2_gdf(tr3,[savepath filename '_cam3_tracers_final.gdf']);
            write2_gdf(tr4,[savepath filename '_cam4_tracers_final.gdf']);
        end
    end
    
    if ~max(cls_flag);
        cls1=read2_gdf([savepath filename '_cam1_clusters.gdf']);
        cls2=read2_gdf([savepath filename '_cam2_clusters.gdf']);
        cls3=read2_gdf([savepath filename '_cam3_clusters.gdf']);
        cls4=read2_gdf([savepath filename '_cam4_clusters.gdf']);

        if (~isempty(cls1) && ~isempty(cls2) && ~isempty(cls3) && ~isempty(cls4) )
            cls1(~ismember(cls1(:,6),cls2(:,6)),:)=[];
            cls1(~ismember(cls1(:,6),cls3(:,6)),:)=[];
            cls1(~ismember(cls1(:,6),cls4(:,6)),:)=[];
            cls2(~ismember(cls2(:,6),cls1(:,6)),:)=[];
            cls3(~ismember(cls3(:,6),cls1(:,6)),:)=[];
            cls4(~ismember(cls4(:,6),cls1(:,6)),:)=[];
            write2_gdf(cls1,[savepath filename '_cam1_clusters_final.gdf']);
            write2_gdf(cls2,[savepath filename '_cam2_clusters_final.gdf']);
            write2_gdf(cls3,[savepath filename '_cam3_clusters_final.gdf']);
            write2_gdf(cls4,[savepath filename '_cam4_clusters_final.gdf']);
        end
    end
    
end