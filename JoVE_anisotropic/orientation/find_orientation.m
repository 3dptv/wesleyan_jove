clear
restoredefaultpath

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
pathlocation = strcat('Z:',slash,'bcole',slash,'JoVE',slash,'MATLABcodes',slash,'JoVE_anisotropic');

% Find the correct paths for the programs that will be called.
addpath(strcat(pathlocation, slash, 'clusters', slash))
addpath(strcat(pathlocation, slash, 'orientation', slash))
addpath(strcat(pathlocation, slash, 'extra', slash))
addpath(strcat(pathlocation, slash, 'data', slash));

filename = 'st001';

% Define paths. Now, "filepath" is the path to the output of find_cluster.
filepath = strcat(pathlocation, slash, 'clusters', slash, 'output', slash);
savepath = strcat(pathlocation, slash, 'orientation', slash, 'results', slash);

% Load in all the .cpv files.
cpv1 = fopen([pathlocation slash 'video_files' slash 'camera1' slash filename '_c1.cpv'],'r');
cpv2 = fopen([pathlocation slash 'video_files' slash 'camera2' slash filename '_c2.cpv'],'r');
cpv3 = fopen([pathlocation slash 'video_files' slash 'camera3' slash filename '_c3.cpv'],'r');
cpv4 = fopen([pathlocation slash 'video_files' slash 'camera4' slash filename '_c4.cpv'],'r');

cpvs = [cpv1, cpv2, cpv3, cpv4];
ncams = length(cpvs);

% Calibration file
camParaCalib = load([pathlocation slash 'calibration' slash 'dynamic_camParaCalib.mat']);
camParaCalib = camParaCalib.camParaCalib;

% Define the number of pixels in each image.
xpix = camParaCalib(1).Npixw;
ypix = camParaCalib(1).Npixh;

cpv_index_good = load([filepath filename '_cpv_index_good.mat']);
cpv_index_good = cpv_index_good.cpv_index_good;

% Set parameters as they were in find_cluster.
cls_min = [300, 300, 300, 300];
cls_max = [2000, 2000, 2000, 2000];

cls_sep = 50; % Separation of clusters required
tr_min = [2, 2, 2, 2]; % min area of tracers (in pixels)
tr_max = [25, 25, 25, 25]; % max area of tracers (in pixels)
threshold = 36; % brightness threshold (max. noise brightness)

fps = 450;

% Create a model object.
model = bc_object;

% Load in the file created from the outputs of find_cluster.
dataspot = strcat(filepath, 'final_tracked', slash, filename, '_clusters_final_tracked.gdf');
data = read2_gdf(dataspot);
data(data(:,19)==1,:)=[];

% Separate out a couple of helpful columns of the data file.
times = data(:,5);
cents3d = data(:,2:4);

nframes = length(data);

% Separate out the tracks. Each version of MATLAB seems to do this
% differently, which is why there are a couple of different attempts at it
% in this program. Each seems to work for a different version.
[C ia ic] = unique(data(:,1));
tracklength = [ia;nframes+1];
tracklength = diff(tracklength);
if length(data(data(:,1)==0,1)) ~= tracklength(1)
    [C ia ic] = unique(data(:,1),'first');
    tracklength=ia;
    tracklength(2:end)=diff(tracklength);
    display('First fix');
end
if length(data(data(:,1)==0,1)) ~=tracklength(1)
    [C ia ic] = unique(data(:,1),'first');
    tracklength = [ia;nframes+1];
    tracklength = diff(tracklength);
    display('Second fix')
end

% Create the matrix where the final optimized results will go.
data_opt = zeros(nframes,31);
% Define at which frame and track you choose to start. Track numbering
% starts at 0, track indexing in this program starts at 1.
firstframe = 1;
firsttrack = 1;

% create mask & gutter
xi = [50 50 1230 1230];
yi = [974 100 100 974];
mask=poly2mask(xi,yi,double(ypix),double(xpix));
mask = uint8(mask);
gutter = find(mask==0);

nonlintic=tic;
iframe=firstframe;

% ---------
frametime = 0;
nbpline = [28,29,30,31];
    
for trackid = firsttrack:length(tracklength)
    
    % Tell the program that it is starting a new track.
    new = 1;
    
    for itrack=1:tracklength(trackid);
        
        data_opt(iframe,1:27) = [data(iframe,1) zeros(1,6) data(iframe,5:end) zeros(1,5)];
        
        frame = tic;
        display(sprintf('Current track: %d out of %d tracks',trackid,length(tracklength)));
        display(sprintf('Current frame: %d out of %d in total',iframe,nframes));
        display(sprintf('Current frame: %d out of %d in track',itrack,tracklength(trackid)));

        ncams = size(cpvs,2);
        clusters = struct('ind',{[],[],[],[]},'indmod',{[],[],[],[]});
        
        posinfile=data(iframe,[8,11,14,17]);
        center_2d = struct('center',{data(iframe,6:7),data(iframe,9:10),data(iframe,12:13),data(iframe,15:16)});
    
    % The bulk of this code is about finding the cluster in the file and
    % finding the bounds of that cluster of bright pixels on all four
    % cameras.
        for icam=1:ncams
            nclusters=1;
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
            pixcnt = bitshift(p(4),14)+bitshift(p(3),6)+bitshift(p(2),-2);
            pixlist = zeros(pixcnt,3,'uint32');
            dat = fread(cpvfile,[4,pixcnt],'uint8=>uint32');
            last3bitarray=uint32(07)+zeros(1,pixcnt,'uint32');
            pixlist(:,3)= dat(1,:);
            %pixel locations range [0,xpix-1] by [0,ypix -1]
            pixlist(:,1) = dat(2,:)+bitshift(bitand(dat(3,:),last3bitarray),8); % decoding x position
            pixlist(:,2) = bitshift(dat(3,:),-3) + bitshift(dat(4,:),5);% decoding y position
            img = zeros(ypix,xpix,'uint8');
            for i=1:pixcnt
                if(pixlist(i,1) < xpix && pixlist(i,1) > 0 && pixlist(i,2) < ypix && pixlist(i,2) > 0 ...
                        && pixlist(i,3) > threshold)
                    img(pixlist(i,2)+1,pixlist(i,1)+1)=pixlist(i,3);
                end
            end
            
            imgb=im2bw(img,0.01);
            pix_info=regionprops(imgb,img,'Centroid','Area','PixelValues','PixelList','PixelIdxList');
            center=cat(1,pix_info.Centroid);
            area=cat(1,pix_info.Area);
            ind = find(area>cls_min(icam) & area<cls_max(icam));
            if new
                if ~isempty(ind)
                    close_cls=0;
                    for idm=1:length(ind)
                        if norm(center(ind(idm),1:2)-center_2d(icam).center(1:2))<cls_sep;
                            close_cls=1;
                        end
                    end
                    for idz=length(ind):-1:1
                        if (~isempty(intersect(pix_info(ind(idz)).PixelIdxList,gutter)) && ~close_cls);
                            ind(idz)=[];
                        end
                    end
                end
                if ~isempty(ind)
                    for idz=1:length(ind)
                        clusters(icam).ind(nclusters:nclusters+size(pix_info(ind(idz)).PixelList,1)-1,...
                            1:2)=pix_info(ind(idz)).PixelList; 
                        clusters(icam).ind(nclusters:nclusters+size(pix_info(ind(idz)).PixelList,1)-1,...
                            3)=pix_info(ind(idz)).PixelValues;
                        nclusters = nclusters+size(pix_info(ind(idz)).PixelList,1);
                    end
                end
            else

                % Create the model and orient it according to the Euler
                % angles found for the previous frame.
                % r_R is the center of the particle, which is first guessed
                % to be the geometric center of the bright pixels. Later,
                % the center position is optimized along with the Euler
                % angles to best fit the model.
                r_R = eul_cntr(4:6);
                x_R = single((calibProj_Tsai(camParaCalib(icam),r_R)));
                x_R = repmat(x_R,4,1);
                center_model = x_R;
                
                if ~isempty(ind)
                    for idz=length(ind):-1:1
                        close_cls=0;
                        for idm=1:size(center_model,1)
                             if norm(center(ind(idz),1:2)-center_model(idm,1:2))<cls_sep;
                                 close_cls=1;
                             end
                        end
                        if (~isempty(intersect(pix_info(ind(idz)).PixelIdxList,gutter)) && ~close_cls);
                            ind(idz)=[];
                        end
                    end
                end
                if ~isempty(ind)
                    for idz=1:length(ind)
                        clusters(icam).ind(nclusters:nclusters+size(pix_info(ind(idz)).PixelList,1)-1,...
                            1:2)=pix_info(ind(idz)).PixelList; 
                        clusters(icam).ind(nclusters:nclusters+size(pix_info(ind(idz)).PixelList,1)-1,...
                            3)=pix_info(ind(idz)).PixelValues;
                        nclusters = nclusters+size(pix_info(ind(idz)).PixelList,1);
                    end
                end
            end
            clusters(icam).indmod=clusters(icam).ind;
            data_opt(iframe,nbpline(icam)) = length(clusters(icam).ind);

        end
                
        % Now we do a nonlinear least squares fit to find the orientation.
            nonlin0=tic;
            fval = zeros(1,2);
            if new
                % First guess of Euler angles for a new object is a random 1 x 3 matrix.
                % Otherwise the previous frame's Euler angles are used.
                eul = rand(1,3)*2*pi;
                % First guess of the 3D position is the geometric center of
                % the bright pixels.
                eul_cntr = [eul cents3d(iframe,:)];
            end

            if new
                
                % "fval" is the residual: a way of checking how well the
                % model is matched to the bright pixels on all four
                % cameras. We separate out the residuals from the first 
                % frame of a track (fval(1)) and the rest of a track (fval(2)). 
                [eul_cntr,fval(1),fval_new,out,xflag] = bc_sk_nonlinearopt(eul_cntr,clusters,...
                    camParaCalib,model);
            else 
                if ~mod(itrack,10)
                    new=1;
                end
                [eul_cntr,fval(2),fval_new,out,xflag] = bc_sk_nonlinearopt(eul_cntr,clusters,...
                    camParaCalib,model);
            end
            
            frametime=frametime+toc(frame);
            display(sprintf('\tThis frame: %f \n\tAverage time/frame: %f sec \n\tExpected to finish: %f min',...
                toc(nonlin0), frametime/iframe,((nframes-iframe)*(frametime/iframe))/60));

%             Save the results of the optimization with information about
%             the frame.
             framenum = times(iframe)*450;
             data_opt(iframe,[2:7,23:27]) = [eul_cntr(4:6) eul_cntr(1:3) fval(1:2) framenum fval_new];
            
            iframe=iframe+1;
            new = 0;
    end
    save([savepath filename '_opt_temp.mat'], 'data_opt');
end

nonlint=toc(nonlintic);

% As in find_cluster, display how long the code took to run.
hrs = floor(nonlint/3600);
minutes = floor((nonlint-(hrs*3600))/60);
seco = mod(nonlint-((hrs*3600)+(minutes*60)),60);
display(sprintf('Total time: %5.3f hours, %5.3f minutes, %5.3f sec.', hrs,minutes,seco));
display(sprintf('Nonlinear fit total time: %f seconds',nonlint));

save([savepath filename '_opt.mat'], 'data_opt');
write2_gdf(data_opt,[savepath filename '_opt.gdf']);

% Check the Euler angles to ensure that the algorithm did not find
% different symmetric orientations of the tetrad for consecutive frames.
data_fixed = bc_min_rotation_finder(data_opt);
display('Eulfix complete');

save([savepath filename '_fixed.mat'], 'data_fixed');
write2_gdf(data_fixed,[savepath filename '_fixed.gdf']);
fclose('all');

% Now we splice the tracks together and calculate the rotation rates.
rrates = tic;
data_spliced = splice_tracks(data_fixed);
splicedpath = strcat(savepath, 'spliced', slash);
save([splicedpath filename '_spliced.mat'], 'data_opt');
display('Done splicing tracks')

% Calculating rotation rates goes over many fit lengths, so it makes the
% most sense to save inside of that program.
rratepath = strcat(savepath, 'rotationrates', slash);
calc_rotation(data_spliced, filename, rratepath);

% Display how long these two steps took.
rates = toc(rrates);
hrs = floor(rates/3600);
minutes = floor((rates-(hrs*3600))/60);
seco = mod(rates-((hrs*3600)+(minutes*60)),60);
display(sprintf('Total time: %5.3f hours, %5.3f minutes, %5.3f sec.', hrs,minutes,seco));

% Display the total run time.
ttime = rates + nonlint;
hrs = floor(ttime/3600);
minutes = floor((ttime-(hrs*3600))/60);
seco = mod(ttime-((hrs*3600)+(minutes*60)),60);
display(sprintf('Total program run time: %5.3f hours, %5.3f minutes, %5.3f seconds.', hrs, minutes, seco));

% And display at what time it all finished.
c = clock;
display(sprintf('Finished orientations at %d:%d:%d',c(4),c(5),c(5)));

