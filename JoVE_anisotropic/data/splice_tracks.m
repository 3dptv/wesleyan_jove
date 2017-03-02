function data_spliced = splice_tracks(data_opt)
% This function's purpose is to make sure that while a single tetrad is in
% view, it is registered as a single continuous track. The function does this by
% checking the end of one track against the beginning of another to see if
% the two tetrads are in similar places with similar orientations at
% approximately the same time. 

% Get rid of any empty rows that might be in the data file.
data_opt(data_opt(:,11)==0,:)=[];

[~,ntracks2check] = mode(data_opt(:,25));
%  While find_orientation does very well for the vast majority (~90%) of
%  frames, there are of course times where the model is an imperfect fit.
%  The following lines are our way of finding these times without having to
%  comb through the data by hand. The "cuts" below are also imperfect, but
%  they do remove of a significant fraction of the bad measurements. 
%  Other experimenters, of course, would find other ways of removing bad
%  measurements.
opterror24 = 2000;
opterror26 = 4000;
nbrp = data_opt(:,28:31);
m = size(nbrp,1);
[mbr, j] = min(nbrp,[],2);
nbrp((1:size(nbrp,1))+(j-1).'*m) = 10000;
[smbrp,j] = min(nbrp,[],2);
nbrp((1:size(smbrp,1))+(j-1).'*m) = 10000;
tmbrp = min(nbrp,[],2);
for k = length(tmbrp):-1:1
    if tmbrp(k) < 1100
        if data_opt(k,24)>opterror24 && data_opt(k,26)>opterror26
            data_opt(k,:) = [];
        elseif data_opt(k,23)>opterror24 && data_opt(k,26)>opterror26
            data_opt(k,:) = [];
        end
    elseif data_opt(k,26) > opterror26 && data_opt(k,27) > opterror26
        data_opt(k,:) = [];
    end
end
    

% Now the tracks can be spliced together.
% This algorithm will not work with the first track one frame long, and we
% cannot calculate rotation rates of tracks one frame long, so we can
% delete the first track if it is one frame.
kk = data_opt(1,1);
if length(data_opt(data_opt(:,1)==kk,1)) == 1
    data_opt(1,:) = [];
end
kk = data_opt(1,1);
if length(data_opt(data_opt(:,1)==kk,1)) == 1
    data_opt(1,:) = [];
end

dt = 1/450;

% These separation parameters are chosen by the experimenter.
max_psep = 1.5;
max_osep = 1;
max_tsep = 30*dt;

% Separate out the tracks. Each version of MATLAB seems to do this
% differently, which is why there are a couple of different attempts at it
% in this program. Each seems to work for a different version.
[C ia ic] = unique(data_opt(:,1));
trackidx = ia(2:end);
tracklength = [ia;length(data_opt)];
tracklength = diff(tracklength);
if length(data_opt(data_opt(:,1)==0,1)) ~= tracklength(1)
    [C ia ic] = unique(data_opt(:,1),'first');
    trackidx = ia;
    tracklength=ia;
    tracklength(2:end)=diff(tracklength);
    display('First trackidx fix');
end
if length(data_opt(data_opt(:,1)==0,1)) ~=tracklength(1)
    [C ia ic] = unique(data_opt(:,1),'first');
    trackidx = ia(2:end);
    tracklength = ia;
    tracklength = diff(tracklength);
    display('Second trackidx fix')
end

if length(data_opt(data_opt(:,1)==0,1)) ~= trackidx(1)
    trackidx = trackidx - 1;
end

for j = 1:ntracks2check+1

    psep = sqrt(sum((data_opt(trackidx(j:end-1)+1,2:4)-data_opt(trackidx(1:end-j),2:4)).^2,2));
    osep = abs(data_opt(trackidx(j:end-1)+1,5:7)-data_opt(trackidx(1:end-j),5:7));
    tsep = data_opt(trackidx(j:end-1)+1,8)-data_opt(trackidx(1:end-j),8);
    k=0;
    for i=1:length(tsep)-1
        if (psep(i)<max_psep && tsep(i)<max_tsep && osep(i,1)<max_osep && osep(i,2)<max_osep...
                && osep(i,3)<max_osep)
            data_opt(trackidx(i+j-1)+1:trackidx(i+j),1) = data_opt(trackidx(i),1);
            k=k+1;
        end
    end

end

% Now, re-sort data_opt so that the first column (track number) is in
% order after the renumbering that happened above.
[~,b] = sort(data_opt(:,1));
data_opt = data_opt(b,:);

data_spliced = data_opt;

end