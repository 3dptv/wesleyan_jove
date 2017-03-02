function cpv_reader(filename,min_frame,max_frame,k)
% Script to read in a .cpv file and extract pixel lists or images for a
% specified range of frames.

% The following line will need to be changed to fit the file structure on
% your computer.
pathlocation = 'Z:\bcole\JoVE\MATLABcodes';

min_frame_num = min_frame;
max_frame_num = max_frame;

plotter = zeros(max_frame_num-min_frame_num+1,2,4);
for icam=1:4
    
cpv = fopen(strcat(pathlocation,'\video_files\camera',int2str(icam),'\',filename,'.cpv'),'r');

show = 1;
verbose = 1;


%The next matrix is going to be used to plot bright frames against frame
%number for each camera. A row for each frame (+1 to avoid overflows), 
%2 columns (frames, pixels), 4 cameras. 
plotter(:,1,icam) = min_frame_num:max_frame_num;

[cpvheader,rcnt] = fread(cpv,20,'*uint8');
if rcnt ~= 20
    warning('Error reading cpv file header');
end
xpix = uint32(cpvheader(5))+uint32(cpvheader(6))*256;
ypix = uint32(cpvheader(7))+uint32(cpvheader(8))*256;
if verbose
    fprintf('CPV file version: %d\n',cpvheader(1));
    fprintf('Number of Lines in real-time image compression output: %d\n',cpvheader(2));
    fprintf('Image size (x,y): %d by %d pixels \n',xpix, ypix);
end
%other info in file header which is ignored here include exposure, rate, and gain

%allocate an image to fill with pixel brightness data
img = zeros(ypix,xpix,'uint8'); %matlab uses (row,column) indexing
%fseek(cpv,20,'bof');  This command would skip reading the file header

while ~feof(cpv)  %While not at the end of the cpv file
    % Read 4 bytes which contain the frame number.
    frame_num=fread(cpv,1,'*uint32');
    
    % Read the next four bytes for the information about the number of
    % bright pixels and decode it. pixcnt contains the number of bright
    % pixels in the current frame.
    [p,rcnt] = fread(cpv,4,'uint8=>uint32');
    if rcnt ~= 4
        warning('Error in frame header');
    else
        % The data in the second 4 bytes of the frame header is:
        threshold=p(1);  %Should be same in every frame, but is stored here anyway
        fifo_overflow =bitget(p(2), 1);
        line_overflow =bitget(p(2),2);
        pixcnt = bitshift(p(4),14)+bitshift(p(3),6)+bitshift(p(2),-2);
     
        
        if frame_num < min_frame_num
            fseek(cpv,pixcnt*4,'cof');
        elseif frame_num > max_frame_num
            break;
        else
            if verbose
                fprintf('Read frame %d with %d pixels.\t fifo overflow: %d \t line_overflow: %d\n',...
                    frame_num,pixcnt,fifo_overflow,line_overflow);
            end

            pixlist = zeros(pixcnt,3,'uint32');
            dat = fread(cpv,[4,pixcnt],'uint8=>uint32');
            last3bitarray=uint32(07)+zeros(1,pixcnt,'uint32');
            
            pixlist(:,3)= dat(1,:);
            %pixel locations range [0,xpix-1] by [0,ypix -1]
            pixlist(:,1) = dat(2,:)+bitshift(bitand(dat(3,:),last3bitarray),8); % decoding x position
            pixlist(:,2) = bitshift(dat(3,:),-3) + bitshift(dat(4,:),5);% decoding y position
            
            %For brightness vs. frame plot
            plotter(frame_num-min_frame_num+1,2,icam)=pixcnt;
            
        %optionally show the image
            if show
                
            img = zeros(ypix,xpix,'uint8');

            for i=1:pixcnt
                if(pixlist(i,1) < xpix && pixlist(i,2) < ypix && pixlist(i,3) >= threshold)
                    img(pixlist(i,2)+1,pixlist(i,1)+1)=pixlist(i,3);
                else
                    if verbose
                        fprintf('Bad pixel:  x=%d, y=%d, brightness=%d\n',pixlist(i,:));
                    end
                end
            end         
            
                figure(k);imagesc(img,[0,255]);
                pause(0.01)
            end

        end
    end
end
% figure(k+1);plot(plotter(:,1,icam),plotter(:,2,icam),'r:')
end
figure(k+1)
plot(plotter(:,1,1),plotter(:,2,1),'r:',plotter(:,1,2),plotter(:,2,2),'b:',...
    plotter(:,1,3),plotter(:,2,3),'g:',plotter(:,1,4),plotter(:,2,4),'k:')
legend(['Camera 1';'Camera 2';'Camera 3';'Camera 4'])


%fclose(cpv);

end