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
pathlocation = strcat('Z:',slash,'bcole',slash,'JoVE',slash,'MATLABcodes');

% Find the correct path for the programs that will be called.
addpath(strcat(pathlocation, slash, 'clusters', slash))

% These lines declare, in order: Where to find the calibration
% file, where to find the video file for one camera, and where to put the
% output of this function.
addpath(genpath(strcat(pathlocation, slash, 'calibration', slash)));
filepath = strcat(pathlocation, slash, 'video_files', slash, 'camera1', slash);
savepath = strcat(pathlocation, slash, 'clusters', slash, 'output', slash');

% Because there are usually many video files per data set, the following
% three lines are used to find and identify which file we want to be looking at.
files = dir(filepath);
files(1:3)=[];
files(2:end)=[];

% Displays the file name as a final check that this is the one you want to analyze.
display(files.name)
pause(5.0)

filenames = char.empty(0,4);
filenames{1} = strcat(files.name(1:5),'_c1.cpv');
filenames{2} = strcat(files.name(1:5),'_c2.cpv');
filenames{3} = strcat(files.name(1:5),'_c3.cpv');
filenames{4} = strcat(files.name(1:5),'_c4.cpv');

% Now, where to find all four video files.
path1 = strcat(pathlocation, slash, 'video_files', slash, 'camera1', slash);
path2 = strcat(pathlocation, slash, 'video_files', slash, 'camera2', slash);
path3 = strcat(pathlocation, slash, 'video_files', slash, 'camera3', slash);
path4 = strcat(pathlocation, slash, 'video_files', slash, 'camera4', slash);

% Set parameters to use when searching for particles.
min_nbrp = [300, 300, 300, 300]; % min number of bright pixels required on each camera
cls_min = [300, 300, 300, 300]; % min number of bright pixels for a cluster
cls_max = [2000, 2000, 2000, 2000]; % max number of bright pixels for a cluster
tr_min = [2, 2, 2, 2]; % min area of tracers (in pixels)
tr_max = [25, 25, 25, 25]; % max area of tracers (in pixels)
threshold = 36; % brightness threshold (max. noise brightness)
    
startf = 1; % starting frame number
endf = 4950; % last frame number
nframes = 4950; % total number of frames

% load cpv files
[pathstr,filename,ext] = fileparts([path1, filenames{1}]);
cpv1 = fopen([pathstr,slash,filename,ext],'r');
[pathstr,filename,ext] = fileparts([path2, filenames{2}]);
cpv2 = fopen([pathstr,slash,filename,ext],'r');
[pathstr,filename,ext] = fileparts([path3, filenames{3}]);
cpv3 = fopen([pathstr,slash,filename,ext],'r');
[pathstr,filename,ext] = fileparts([path4, filenames{4}]);
cpv4 = fopen([pathstr,slash,filename,ext],'r');

cpvs = [cpv1, cpv2, cpv3, cpv4];

% read in cpv files
display(['Reading in cpv files: ' filename]);
read_cpv = tic;

cpv_index = sk_read_cpv(cpvs,nframes);

save([savepath filename '_cpv_index.mat'], 'cpv_index');
display(sprintf('Elapsed time: %5.3f sec\n', toc(read_cpv)));

display('Removing bad frames...');
remove_bad = tic;

[cpv_index_good, cpv_info] = sk_remove_bad_frames(cpv_index,startf,endf,min_nbrp);
pause(3.0);

save([savepath filename '_cpv_index_good.mat'], 'cpv_index_good');
save([savepath filename '_cpv_info.txt'], 'cpv_info','-ascii');
display(sprintf('Elapsed time: %5.3f min\n', toc(remove_bad)/60));

display('Now searching for tracers and clusters...')
search = tic;

% Now search for tracers and particles
bc_sk_find_particles(cpvs,cpv_index_good,cls_min,cls_max,tr_min,tr_max,threshold,savepath,filename);

% Now that it's all done, display how long that the program took to
% run.
s = toc(search);
hrs = floor(s/3600);
minutes = floor((s-(hrs*3600))/60);
seco = mod(s-((hrs*3600)+(minutes*60)),60);
display(sprintf('Total time: %5.3f hours, %5.3f minutes, %5.3f sec.', hrs,minutes,seco));

fclose all;

% And display at what time it finished.
c = clock;
display(sprintf('Finished at %d:%d:%d',c(4),c(5),c(5)));