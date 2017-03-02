function bc_sk_find_frame_colors(data_opt,filename,min_frame,max_frame,imnum)

% This function takes in the output from find_orientation or splice_tracks
% and plots the model on top of the image data, with each arm of the
% model as a line of a different color.
% The other inputs are, in order: filename (as a string: 'st001'), the
% first frame you want to view, the final frame you want to view, the
% figure number you want to use.

% Mac computers address files with the slashes going in the opposite
% direction of the slashes on Windows computers.
mac = 0; % =0 if not on a Mac, =1 if on a Mac.
if mac
    slash = '/';
elseif ~mac
    slash = '\';
end

% The following line will need to be changed to fit the file structure on
% your computer.
pathlocation = strcat('Z:',slash,'bcole',slash,'JoVE',slash,'MATLABcodes');

show = 1; % Optionally show the image.

% Model arm length.
hl = 4;

% Dimensions of the figure when zoomed in on the tetrad. These numbers are
% chosen to match the aspect ratio of the cameras to minimize distortions.
wsize = [93.75,75];


camParaCalib = load([pathlocation slash 'calibration' slash 'dynamic_camParaCalib.mat']);
camParaCalib = camParaCalib.camParaCalib;

% Check to make sure that the frames chosen for viewing are frames that are
% in the file.
frames = data_opt(:,25);
if ~any(frames==min_frame)
    display(sprintf('No data at frame %d, searching for lowest close frame.',min_frame))
    frameposs = min_frame:max_frame;
    for iframe = frameposs
        if any(frames==iframe)
            display(sprintf('Found a match at frame %d. Using that as min_frame.',iframe))
            min_frame = iframe;
            break;
        end
    end
end
minframeline = min(find(frames==min_frame));
if ~any(frames==max_frame)
    display(sprintf('No data at frame %d, searching for highest close frame.',max_frame))
    frameposs = max_frame:-1:min_frame;
    for iframe = frameposs
        if any(frames==iframe)
            display(sprintf('Found a match at frame %d. Using that as max_frame.',iframe))
            max_frame = iframe;
            break;
        end
    end
end
maxframeline = max(find(frames==max_frame));

min_frame = minframeline;
max_frame = maxframeline;

% Load in all the cpv files.
cpv1 = fopen([pathlocation slash 'video_files' slash 'camera1' slash filename '_c1.cpv'],'r');
cpv2 = fopen([pathlocation slash 'video_files' slash 'camera2' slash filename '_c2.cpv'],'r');
cpv3 = fopen([pathlocation slash 'video_files' slash 'camera3' slash filename '_c3.cpv'],'r');
cpv4 = fopen([pathlocation slash 'video_files' slash 'camera4' slash filename '_c4.cpv'],'r');

cpvs = [cpv1, cpv2, cpv3, cpv4];
ncams = length(cpvs);
model = bc_object;

for iframe=min_frame:max_frame
    
posinfile=data_opt(iframe,[11,14,17,20]);
pixcnts = data_opt(iframe,28:31);

% make the model
eul = data_opt(iframe,5:7);
cg = data_opt(iframe,2:4);
ori  = gg_ori(eul);
r_L = zeros(length(model.r_L),3);
r_R = zeros(length(model.r_R),3);
for i=1:length(model.r_L)
    r_L(i,:) = (ori*model.r_L(i,:)')';
    r_R(i,:) = (ori*model.r_R(i,:)')';
end
r_L = bsxfun(@plus, r_L*hl,cg);
r_R = bsxfun(@plus, r_R*hl,cg);

for icam=1:ncams
    
%     Much of this is very similar to read_cpv, but now the intent is to
%     plot the image.
    cpvfile = cpvs(icam);
    frewind(cpvfile);
    [cpvheader,rcnt] = fread(cpvfile,20,'*uint8');
        if rcnt ~= 20
            warning('Error reading cpv file header');
        end
    xpix = uint32(cpvheader(5))+uint32(cpvheader(6))*256;
    ypix = uint32(cpvheader(7))+uint32(cpvheader(8))*256;
    fseek(cpvfile,posinfile(icam),'bof');
    currentFrameNum = fread(cpvfile,1,'*uint32');

    [p, rcnt] = fread(cpvfile,4,'uint8=>uint32');
    threshold=p(1);  %Should be same in every frame, but is stored here anyway
    pixcnt = bitshift(p(4),14)+bitshift(p(3),6)+bitshift(p(2),-2);

    pixlist = zeros(pixcnt,3,'uint32');
    dat = fread(cpvfile,[4,pixcnt],'uint8=>uint32');
    last3bitarray=uint32(07)+zeros(1,pixcnt,'uint32');
    pixlist(:,3)= dat(1,:);
    % pixel locations range [0,xpix-1] by [0,ypix -1]
    pixlist(:,1) = dat(2,:)+bitshift(bitand(dat(3,:),last3bitarray),8); % decoding x position
    pixlist(:,2) = bitshift(dat(3,:),-3) + bitshift(dat(4,:),5);% decoding y position
    img = zeros(ypix,xpix,'uint8');
    for i=1:pixcnt
        if(pixlist(i,1) < xpix && pixlist(i,1) > 0 && pixlist(i,2) < ypix &&...
                pixlist(i,2) > 0 && pixlist(i,3) >= threshold)
            img(pixlist(i,2)+1,pixlist(i,1)+1)=pixlist(i,3);
        end
    end
    
    if show
        center = (calibProj_Tsai(camParaCalib(icam),cg));
        x_R = (calibProj_Tsai(camParaCalib(icam),r_R));
        x_L = (calibProj_Tsai(camParaCalib(icam),r_L));
        xset = [x_L(1,1);x_R(1,1);x_L(2,1);x_R(1,1);x_L(3,1);x_R(1,1);x_L(4,1);x_R(1,1)];
        yset = [x_L(1,2);x_R(1,2);x_L(2,2);x_R(1,2);x_L(3,2);x_R(1,2);x_L(4,2);x_R(1,2)];

        figure(imnum);subplot(2,2,icam),imagesc(img),colormap(gray);
        hold on

        % Many different options for understanding and visualizing how well the
        % model fits the data.
        plot(x_R(1,1),x_R(1,2),'wo','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',5);

    % The following eight lines are in pairs. The first of each pair shows the
    % model data as a line from the center to the end of the arm. The second
    % shows just the endpoint of the model arm as a marker.
        plot(xset(1:2),yset(1:2),'-r','LineWidth',4)
    %     plot(x_L(1,1),x_L(1,2),'wo','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
        plot(xset(3:4),yset(3:4),'-b','LineWidth',4)
    %     plot(x_L(2,1),x_L(2,2),'wo','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',5);
        plot(xset(5:6),yset(5:6),'-g','LineWidth',4)
    %     plot(x_L(3,1),x_L(3,2),'wo','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',5);
        plot(xset(7:8),yset(7:8),'-k','LineWidth',4)
    %     plot(x_L(4,1),x_L(4,2),'wo','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);

    %     title(num2str(frames(iframe))); % Uncomment to have the title of all
    %     of the subplot images be the frame number.

        title(sprintf('Camera %d',icam)); % Uncommment to have the
    %     title of the subplot images be the camera number.

        % Zoom in on the object to get a better view of how well the model is
        % matched.
        axis([center(1)-wsize(1),center(1)+wsize(1),center(2)-wsize(2),center(2)+wsize(2)]);

    %     To get rid of pixel counts on the figure axes, uncomment the following
    %     four lines.
    %     set(gca, 'XTickLabelMode', 'Manual')
    %     set(gca, 'XTick', [])
    %     set(gca, 'YTickLabelMode', 'Manual')
    %     set(gca, 'YTick', [])

        hold off
    end
end

pause(.01)
% Uncomment the following lines to display the number of pixels visible on
% each of the four cameras for the current frame.
fprintf('Read frame %d with\t %d\t %d\t %d\t %d\t pixels.\n',frames(iframe),pixcnts(1),...
    pixcnts(2),pixcnts(3),pixcnts(4));

% Uncomment the following line to display the current frame number.
% fprintf('Read frame %d.\n',frames(iframe));

end

end