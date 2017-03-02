function calc_rotation(data_opt, filename, savepath)
% This function calculates an average rotation rate over a range of numbers
% of consecutive frames. See Ref [1] (Marcus, 2014), Figure 4, for the
% fit-length dependence of rotation rate measurements of jacks and crosses.

message = sprintf('Now calculating rotation rates for %s.',filename);
display(message);

% We go out to fit lengths of about 2.5 Kolmogorov times. The videos
% included with this publication have a different Kolmogorov time than is
% reported in the text because they were captured with the grid oscillating
% at a higher frequency and with a slightly different salt solution.
hset = 1:43;

d = data_opt;

for ii = 1:length(hset)

    h = hset(ii);
    
    mintracklength = round(1.25*h);
    data_opt = d;
    
    % Get rid of short tracks.
    [C1 ia1 ic1] = unique(data_opt(:,1),'first');
    trackidx=ia1;
    tracklength=diff(trackidx);
    short=find(tracklength<mintracklength);
    to_del=zeros(size(trackidx));
    for i=1:length(short)
        to_del(trackidx(short(i)):trackidx(short(i)+1)-1)=1;
    end
    data_opt(to_del==1,:)=[];
    [C1 ia1 ic1] = unique(data_opt(:,1),'first');
    trackidx=ia1;
    tracklength1=diff(trackidx);

    tracklength = tracklength1;

    % Now we can start calculating rotation rates.
    ntracks=size(tracklength,1);
    nframes=size(data_opt(:,1),1);

    dt = 1/450;
    fraction_required=.7;

    data_fit = [data_opt(:,1:8) repmat(-1,size(data_opt,1),11) data_opt(:,9:end) zeros(size(data_opt,1),1)];
for trackid=1:length(tracklength);
    display(sprintf('trackID: %d out of %d for halflength %d of %s.\n',...
        trackid, length(tracklength),h,filename));
    eul=data_opt(trackidx(trackid):trackidx(trackid+1)-1,5:7);
    t=data_opt(trackidx(trackid):trackidx(trackid+1)-1,8)./dt; % t should be in frames
    for center=(h+1):tracklength(trackid)-(h)
        % extract the orientations that are within h time steps of the center,
        % call the ones within h time steps a segment.
        beginstep=center-h; %lower index that can be in the fitting window
        endstep=center+h; %upper index
        trange=t(beginstep:endstep);  %first extract the subset that is within h index steps
        eulrange=eul(beginstep:endstep,:);
        tmid=t(center);
        eulmid=eul(center,:);
        f = find(abs(trange-tmid) <= h); %now find the part of the subset that is actually within h time steps.
        frac = length(f)/(2*h+1);

        if frac >= fraction_required
            tseg=trange(f)-tmid;  %tseg contains the segment of time differences of each sample from the center
            eulseg=eulrange(f,:); %eulseg contains the euler angles at the times stored in tseg
            nseg=length(tseg);
            %Initial guess for R0 is center rotation matrix
            R0=gg_ori(eulmid);
            
            % fitseg is a matrix with (3x3x(2h))x3 elements ignore 0 time
            % element
%             p = R0(:,3)';
%             d = -p; % by definition
            fitseg=zeros(3*3*(nseg-1),3);
            datseg=zeros(3*3*(nseg-1),1);
            k=1;
            for i=1:nseg
                Rh=gg_ori(eulseg(i,:));
                if(tseg(i)~=0)
                    for j=1:3
                        % build fitseg such that fitseg*omega=datseg and
                        % find omega via least squares
                        fitseg(k:k+3-1,:)=tseg(i)*[0, -R0(3,j), R0(2,j); R0(3,j) 0 -R0(1,j); -R0(2,j) R0(1,j) 0;];
                        datseg(k:k+3-1)=Rh(:,j)-R0(:,j);
                        k=k+3;
                    end
                end
            end
            % this is where the least squares fit happens: "\"
            w=(fitseg\datseg);
            % now find smoothed euler angles using fitted omega and
            % compute residual
            rot=[0 w(3) -w(2); -w(3) 0 w(1); w(2) -w(1) 0];
            diffsq = 0;
            for i=1:nseg
                if(tseg(i)~=0)
                    datsegR=gg_ori(eulseg(i,:));
                    fitsegR=(tseg(i)*rot+eye(3))*R0;
                    pfit = fitsegR(:,3);
                    dfit = -pfit;
                    eulfit = bc_euler(fitsegR);
                end
                diffsq = diffsq+sum(sum((fitsegR-datsegR).^2));
            end
            % Renormalize pfit and dfit
            pfit = pfit./norm(pfit);
            dfit = dfit./norm(dfit);
            
            % flip w to switch from passive rotation to active rotation
            S = [-1, 0, 0; 0, -1, 0; 0, 0, -1];
            w = S*w;
            wsn = norm(w);
            pdot = cross(w,pfit,1);
            wd = dot(w,dfit,1);
            
            % Here we save the frame number and residual from the vorticity
            % smoothing
            data_fit(trackidx(trackid)+center-1,5:7)=eulfit;
            data_fit(trackidx(trackid)+center-1,9:11)=w./dt;
            data_fit(trackidx(trackid)+center-1,12)=wsn/dt;
            data_fit(trackidx(trackid)+center-1,13)=wd/dt;
            data_fit(trackidx(trackid)+center-1,14:16)=pfit;
            data_fit(trackidx(trackid)+center-1,17:19)=pdot./dt;
            data_fit(trackidx(trackid)+center-1,33)=1;
            data_fit(trackidx(trackid)+center-1,37)=diffsq;

        end
    end
end

% Getting rid of any empty lines that resulted from not being able to fit
% around that frame.
data_fit(data_fit(:,12)==-1,:) = [];

save([savepath filename '_f' num2str(h) '.mat'], 'data_fit');
end

end
