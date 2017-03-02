function [diffsq1,diffsq2] = bc_new_residual(eul,data2d,camParaCalib,param,model,cg)
% This program calculates a final residual value -- a way of measuring the
% difference between our final model and the actual image data. A nonlinear
% fit could be done on one or both of the residuals calculated here, but we
% determined that this was too computationally intensive for the small
% change in final orientations that would result to be worth it.

ncams = length(camParaCalib);
d = param.diamet;
hl = param.armlength;

%make and orient the model.
ori  = gg_ori(eul);

r_L = zeros(length(model.r_L),3);
r_R = zeros(length(model.r_R),3);

for i=1:length(model.r_L)
    r_L(i,:) = (ori*model.r_L(i,:)')';
    r_R(i,:) = (ori*model.r_R(i,:)')';
end
r_L = bsxfun(@plus, r_L*hl,cg);
r_R = bsxfun(@plus, r_R*hl,cg);

diffsq1 = 0;
diffsq2 = 0;
thresh = .5;

for icam=1:ncams

    % Determine the end points of the objects projected onto each camera
    x_R = single((calibProj_Tsai(camParaCalib(icam),r_R))); % FLOAT
    x_L = single((calibProj_Tsai(camParaCalib(icam),r_L)));
    center = (calibProj_Tsai(camParaCalib(icam),cg));
    
    % First we organize the experimental data
    pixd = single(data2d(icam).ind(:,1:2));
    vd = ones(size(data2d(icam).ind(:,3)));
    
    imgd=zeros(1024,1280,'single');
	imgm=zeros(1024,1280,'single');
    
%     Create a large box of pixels around the center.
    xcuts = single([1280:-1:round(center(1)+196),round(center(1)-196):-1:1]);
    ycuts = single([1024:-1:round(center(2)+153),round(center(2)-153):-1:1]);

    xcuts(xcuts<0) = [];
    xcuts(xcuts==0) = [];
    xcuts(xcuts>1280) = [];
    ycuts(ycuts<0) = [];
    ycuts(ycuts==0) = [];
    ycuts(ycuts>1024) = [];
    
    xpix = camParaCalib(icam).Npixw;
    ypix = camParaCalib(icam).Npixh;    
    pos=zeros(size(vd));
    
    pos(:)=sub2ind([ypix,xpix],pixd(:,2),pixd(:,1));
    imgd(pos)=vd;
    
    imgd(:,xcuts) = [];
    imgd(ycuts,:) = [];

%   Now create the list of the pixels within that box to be put into
%   gaussian_intensity. For this, we need the Cartesian product.
    xvals = single(round(center(1)-196):round(center(1)+196));
    xvals(xvals<=0) = [];
    xvals(xvals>1280) = [];
      
    yvals = single(round(center(2)-153):round(center(2)+153));
    yvals(yvals<=0) = [];
    yvals(yvals>1024) = [];
      
    [pixx,pixy] = meshgrid(xvals,yvals);
    pixd_new = [pixx(:),pixy(:)];

    I0 = single(1);
    % Compute the intensity of the model everywhere in the box.
    vm = bc_sk_gaussian_intensity_multiple_rods(pixd_new,x_L,x_R,I0,d);     
    
    for pix = 1:length(vm)
        imgm(pixd_new(pix,2),pixd_new(pix,1)) = vm(pix);
    end
    imgm(:,xcuts) = [];
    imgm(ycuts,:) = [];
%     
    % Now we compute the sum of the difference between the intesities
    % squared at the points computed.
    diffsq1 = diffsq1 + sum(sum((imgm-imgd).^2));
    
    % Then, we do the same thing but after putting a lower bound on the
    % brightness of the model pixels that we are interested in.
    imgm(imgm<thresh) = 0;
    diffsq2 = diffsq2 + sum(sum((imgm-imgd).^2));

end

% To see an output of this function, run bc_new_residual_plot
end
