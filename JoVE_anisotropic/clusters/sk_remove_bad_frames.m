function [cpv_index_good, cpv_info] = sk_remove_bad_frames(cpv_index,startf,endf,min_nbrp)
% Get rid of frames that we don't want to look at or can't look at
% because there is a problem with them.

    cam1=cpv_index(1).data;
    cam2=cpv_index(2).data;
    cam3=cpv_index(3).data;
    cam4=cpv_index(4).data;
    
    tmp = [length(cam1) length(cam2) length(cam3) length(cam4)];
    totlength = tmp;
    

    
    % outside frame range
    cam1(cam1(:,2)<startf,:)=[];
    cam2(cam2(:,2)<startf,:)=[];
    cam3(cam3(:,2)<startf,:)=[];
    cam4(cam4(:,2)<startf,:)=[];
    cam1(cam1(:,2)>endf,:)=[];
    cam2(cam2(:,2)>endf,:)=[];
    cam3(cam3(:,2)>endf,:)=[];
    cam4(cam4(:,2)>endf,:)=[];
    out_of_range = [tmp(1)-length(cam1) tmp(2)-length(cam2) tmp(3)-length(cam3) tmp(4)-length(cam4)];
    display(sprintf('Frames with frame numbers outside %d and %d:',startf, endf));
        display(sprintf('\t camera 1: %5.0f: ',out_of_range(1)));
        display(sprintf('\t camera 2: %5.0f: ',out_of_range(2)));
        display(sprintf('\t camera 3: %5.0f: ',out_of_range(3)));
        display(sprintf('\t camera 4: %5.0f: \n',out_of_range(4)));
    tmp = [length(cam1) length(cam2) length(cam3) length(cam4)];

    % not enough bright pixels
    cam1(cam1(:,3)<min_nbrp(1),:)=[];
    cam2(cam2(:,3)<min_nbrp(2),:)=[];
    cam3(cam3(:,3)<min_nbrp(3),:)=[];
    cam4(cam4(:,3)<min_nbrp(4),:)=[];
    no_bright_pix = [tmp(1)-length(cam1) tmp(2)-length(cam2) tmp(3)-length(cam3) tmp(4)-length(cam4)];
    tmp = [length(cam1) length(cam2) length(cam3) length(cam4)];
    display(sprintf('Number of frames with not enough bright pixels:'));
        display(sprintf('\t camera 1: %5.0f: (required: %d)',no_bright_pix(1),min_nbrp(1)));
        display(sprintf('\t camera 2: %5.0f: (required: %d)',no_bright_pix(2),min_nbrp(2)));
        display(sprintf('\t camera 3: %5.0f: (required: %d)',no_bright_pix(3),min_nbrp(3)));
        display(sprintf('\t camera 4: %5.0f: (required: %d) \n',no_bright_pix(4),min_nbrp(4)));
    display(sprintf('Average number of bright pixels:'));
        display(sprintf('\t camera 1: %5.1f',mean(cam1(:,3))));
        display(sprintf('\t camera 2: %5.1f',mean(cam2(:,3))));
        display(sprintf('\t camera 3: %5.1f',mean(cam3(:,3))));
        display(sprintf('\t camera 4: %5.1f\n',mean(cam4(:,3))));
        
    % fifo_overflow and line_overflow
    cam1(cam1(:,4:5)==1,:)=[];
    cam2(cam2(:,4:5)==1,:)=[];
    cam3(cam3(:,4:5)==1,:)=[];
    cam4(cam4(:,4:5)==1,:)=[];
    fifo = [tmp(1)-length(cam1) tmp(2)-length(cam2) tmp(3)-length(cam3) tmp(4)-length(cam4)];
    tmp = [length(cam1) length(cam2) length(cam3) length(cam4)];
    display(sprintf('Number of frames with too many bright pixels:'));
        display(sprintf('\t camera 1: %5.0f',fifo(1)));
        display(sprintf('\t camera 2: %5.0f',fifo(2)));
        display(sprintf('\t camera 3: %5.0f',fifo(3)));
        display(sprintf('\t camera 4: %5.0f \n',fifo(4)));

    % not on all cameras
    cam1(~ismember(cam1(:,2),cam2(:,2)),:)=[];
    cam1(~ismember(cam1(:,2),cam3(:,2)),:)=[];
    cam1(~ismember(cam1(:,2),cam4(:,2)),:)=[];
    cam2(~ismember(cam2(:,2),cam1(:,2)),:)=[];
    cam3(~ismember(cam3(:,2),cam1(:,2)),:)=[];
    cam4(~ismember(cam4(:,2),cam1(:,2)),:)=[];
    not_all_cams = [tmp(1)-length(cam1) tmp(2)-length(cam2) tmp(3)-length(cam3) tmp(4)-length(cam4)];
    tmp = [length(cam1) length(cam2) length(cam3) length(cam4)];
    display(sprintf('Number of frames with frame numbers that did not appear on other cameras:'));
        display(sprintf('\t camera 1: %5.0f',not_all_cams(1)));
        display(sprintf('\t camera 2: %5.0f',not_all_cams(2)));
        display(sprintf('\t camera 3: %5.0f',not_all_cams(3)));
        display(sprintf('\t camera 4: %5.0f \n',not_all_cams(4)));
        
    % not in order
    for i=(length(cam1)-1):-1:2
        if (cam1(i-1,2)>cam1(i,2) || cam1(i+1,2)<cam1(i,2))
            cam1(i,:)=[];
        end
    end
    for i=(length(cam2)-1):-1:2
        if (cam2(i-1,2)>cam2(i,2) || cam2(i+1,2)<cam2(i,2))
            cam2(i,:)=[];
        end
    end
    for i=(length(cam3)-1):-1:2
        if (cam3(i-1,2)>cam3(i,2) || cam3(i+1,2)<cam3(i,2))
            cam3(i,:)=[];
        end
    end
    for i=(length(cam4)-1):-1:2
        if (cam4(i-1,2)>cam4(i,2) || cam4(i+1,2)<cam4(i,2))
            cam4(i,:)=[];
        end
    end
    not_in_order = [tmp(1)-length(cam1) tmp(2)-length(cam2) tmp(3)-length(cam3) tmp(4)-length(cam4)];
    tmp = [length(cam1) length(cam2) length(cam3) length(cam4)];
    display(sprintf('Number of frames with frame numbers out of order:'));
        display(sprintf('\t camera 1: %5.0f',not_in_order(1)));
        display(sprintf('\t camera 2: %5.0f',not_in_order(2)));
        display(sprintf('\t camera 3: %5.0f',not_in_order(3)));
        display(sprintf('\t camera 4: %5.0f \n',not_in_order(4)));
    display(sprintf('Number of frames left after selection (out of %d):', mean(totlength)));
        display(sprintf('\t camera 1: %5.0f',tmp(1)));
        display(sprintf('\t camera 2: %5.0f',tmp(2)));
        display(sprintf('\t camera 3: %5.0f',tmp(3)));
        display(sprintf('\t camera 4: %5.0f \n',tmp(4)));

    cpv_index_good = struct('data', {cam1,cam2,cam3,cam4});
    cpv_info = [out_of_range, no_bright_pix, fifo, not_all_cams, not_in_order, tmp];
end