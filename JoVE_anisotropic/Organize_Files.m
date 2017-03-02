% This script should be run first before any of the other programs. It will
% arrange the supplemental files into a file folder structure that the
% subsequent programs will be able to use. It uses basic commands but is
% written in MATLAB because all other programs in this submission are in
% MATLAB.

mac = 1; % If using a Mac, set mac = 1.
if mac
    slash = '';%'\';
elseif ~mac
    slash = '\';
end

% In the following line, enter the location of the folder where the
% supplemental files were saved.
% cd Z:\bcole\JoVE\MATLABcodes\

% Each file will have to be done individually, so there will not be many
% comments throughout. In each section of the script, a folder is created
% and the relevant files are moved into it.

% First, make a folder that everything's going to go into.
mkdir JoVE_anisotropic
parent = fullfile('Jove_anisotropic');
movefile('Use_Instructions.txt',parent)

% Start with the calibration files.
mkdir(fullfile(parent,'calibration'));
target = fullfile(parent,'calibration');
movefile('dynamic_camParaCalib.mat',target)
movefile('README_calib.txt',target)

% Then the cluster-finding codes.
mkdir(fullfile(parent,'clusters'));
target = fullfile(parent,'clusters', slash);
movefile('bc_sk_find_particles.m',target)
movefile('find_clusters.m',target)
movefile('read2_gdf.m',target)
movefile('README_clust.txt',target)
movefile('sk_read_cpv.m',target)
movefile('sk_remove_bad_frames.m',target)
movefile('write2_gdf.m',target)

mkdir(fullfile(parent,'clusters',slash,'output'))
movefile('st001_cpv_index_good.mat',...
    fullfile(parent,'clusters',slash,'output',slash))

mkdir(fullfile(parent,'clusters',slash,'output',slash,'final_matched'))
movefile('st001_clusters_final_matched.gdf',...
    fullfile(parent,'clusters',slash,'output',slash,'final_matched',slash))
mkdir(fullfile(parent,'clusters',slash,'output',slash,'final_tracked'))
movefile('st001_clusters_final_tracked.gdf',...
    fullfile(parent,'clusters',slash,'output',slash,'final_tracked',slash))


% Next, the data analysis codes.
mkdir(fullfile(parent,'data'))
target = fullfile(parent,'data', slash);
movefile('calc_rotation.m',target)
movefile('README_data.txt',target)
movefile('splice_tracks.m',target)

% Next, programs to visualize data and models.
mkdir(fullfile(parent,'extra'))
target = fullfile(parent,'extra', slash);
movefile('bc_cross.m',target)
movefile('bc_degenerateorientations_plots.m',target)
movefile('bc_jack.m',target)
movefile('bc_new_residual_plot.m',target)
movefile('bc_sk_find_frame_colors.m',target)
movefile('bc_sk_plot_model.m',target)
movefile('bc_triad.m',target)
movefile('README_extra.txt',target)


% Next, the codes involved in orientation-finding.
mkdir(fullfile(parent,'orientation'))
target = fullfile(parent,'orientation', slash);
movefile('bc_euler.m',target)
movefile('bc_min_rotation_finder.m',target)
movefile('bc_new_residual.m',target)
movefile('bc_object.m',target)
movefile('bc_sk_gaussian_intensity_multiple_rods.m',target)
movefile('bc_sk_leastSqOriPos.m',target)
movefile('bc_sk_leastSqOriPos1.m',target)
movefile('bc_sk_nonlinearopt.m',target)
movefile('calibProj_Tsai.m',target)
movefile('find_orientation.m',target)
movefile('gg_ori.m',target)
movefile('README_ori.txt',target)

mkdir(fullfile(parent,'orientation',slash,'results'))
target = fullfile(parent,'orientation',slash,'results');
movefile('st001_fixed*',target)
movefile('st001_opt*',target)

mkdir(fullfile(parent,'orientation',slash,'results',slash,'spliced'))
movefile('st001_spliced*',...
    fullfile(parent,'orientation',slash,'results',slash,'spliced', slash))

mkdir(fullfile(parent,'orientation',slash,'results',slash,'rotationrates'))
movefile('README_rot.txt',...
    fullfile(parent,'orientation',slash,'results',slash,'rotationrates',slash))

% Finally, the video files and a script for viewing them.
mkdir(fullfile(parent,'video_files'))
target = fullfile(parent,'video_files', slash);
movefile('cpv_reader.m',target)
movefile('README_videos.txt',target)

for i = 1:4
    mkdir(fullfile(parent,'video_files',slash,'camera',int2str(i)));
    target = fullfile(parent,'video_files',slash,'camera',int2str(i),slash);
    movefile(strcat('README_c',int2str(i),'.txt'),target);
    movefile(strcat('st001_c',int2str(i),'.cpv'),target);
end



