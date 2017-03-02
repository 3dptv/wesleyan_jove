function int_out=bc_sk_gaussian_intensity_multiple_rods(pix_xy,x_L,x_R,I0,armwidth)
% This function projects the model tetrad onto the camera image and
% computes what the brightness of the model would be if seen from that
% camera. This corresponds to step 5.1.3.1 of the Protocol.

% Inputs:  pix_xy:  (N_pixels by 2) matrix containing 2D (x,y) coordinates of bright pixels from data images
%             X_L:  (N_arms by 2) matrix containing 2D (x,y) coordinates of Left ends of the arms
%             X_R:  (N_arms by 2) matrix containing 2D (x,y) coordinates of Right ends of the arms
%                         (center of the tetrad)
%              I0:   Maximum intensity that should be used in returned intensity
% Outputs:int_out:  (N_pixels by 1) matrix with the intensity at each pixel in pix_xy

narms = single(size(x_L,1));
npix=single(size(pix_xy,1));
intensity=zeros(npix,narms+1,'single');

pmid = 0.5*(x_L+x_R); % center point of each arm
armlength = sqrt(sum((x_L-x_R).^2,2))/2;  % length of each arm
theta=atan2(x_R(:,2)-x_L(:,2),x_R(:,1)-x_L(:,1));  % angle of each arm

for i=1:narms
    arm_mid = pmid(i,:);
    rot2D=[cos(theta(i)),sin(theta(i));-sin(theta(i)),cos(theta(i))];
    arm_coords=(rot2D * ([pix_xy(:,1)-arm_mid(1),pix_xy(:,2)-arm_mid(2)])')';
%     Gaussian times Fermi:
    intensity(:,i)=I0*exp(-arm_coords(:,2).^2/(1*armwidth^2)).*(1./(1+exp((abs(arm_coords(:,1))-armlength(i))./(.25*armwidth))));
end

% Putting a circle at the center of the tetrad so that there isn't an
% intensity dropoff in the center, where often the image is brightest.

centered_xy = pix_xy - repmat(x_R(1,:),npix,1);
xcent = centered_xy(:,1)';
ycent = centered_xy(:,2)';

% Fermi intensity across the circle:
intensity(:,narms+1) = I0*1./(1+exp(((xcent.^2+ycent.^2)-armwidth^2)/(1.75*armwidth)));

int_out = max(intensity,[],2);
end
