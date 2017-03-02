function [cpv_index] = sk_read_cpv(cpvs,nframes)
% Read through the videos and identify places to look for the clusters.

    ncams = length(cpvs);
    cpv_index = struct('data',{zeros(nframes,5),zeros(nframes,5),zeros(nframes,5),zeros(nframes,5)});

    for icam=1:ncams
        cpv=cpvs(icam);
        frewind(cpv); % Start at the beginning of the file

        [cpvheader,rcnt] = fread(cpv,20,'*uint8');
        if rcnt ~= 20
            warning('Error reading cpv file header');
        end
        xpix = uint32(cpvheader(5))+uint32(cpvheader(6))*256;
        ypix = uint32(cpvheader(7))+uint32(cpvheader(8))*256;

        display(sprintf('\t CPV file version: %d',cpvheader(1)));
        display(sprintf('\t Number of Lines in real-time image compression output: %d',cpvheader(2)));
        display(sprintf('\t Image size (x,y): %d by %d pixels\n',xpix, ypix));

        posCurrent=zeros(nframes,1);
        currentFrameNum = zeros(nframes,1);
        pixcnt = zeros(nframes,1);
        fifo_overflow = zeros(nframes,1);
        line_overflow = zeros(nframes,1);

        for iframe=1:nframes
            posCurrent(iframe) = ftell(cpv);
            currentFrameNum(iframe) = fread(cpv,1,'*uint32');
            [p, rcnt] = fread(cpv,4,'uint8=>uint32');
            if feof(cpv)
                break;            
            elseif rcnt ~= 4
                error('Error reading frame #%d',currentFrameNum(iframe));
            end
            fifo_overflow(iframe)=bitget(p(2),1);
            line_overflow(iframe)=bitget(p(2),2);
            pixcnt(iframe) = bitshift(p(4),14)+bitshift(p(3),6)+bitshift(p(2),-2);
            fseek(cpvs(icam),pixcnt(iframe)*4,'cof');

        end

        cpv_index(icam).data(:,1) = posCurrent; 
        cpv_index(icam).data(:,2) = currentFrameNum;
        cpv_index(icam).data(:,3) = pixcnt;
        cpv_index(icam).data(:,4) = fifo_overflow;
        cpv_index(icam).data(:,5) = line_overflow;
    end
end