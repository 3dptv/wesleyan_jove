function bc_new_residual_plot(eul,data2d,camParaCalib,param,model,cg)

% This function shows the way the data is interpreted in bc_new_residual.

ncams = length(camParaCalib);
d = param.diamet;
hl = param.armlength;

%make the model
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

    % Determine the end points of the objects projected onto each camera
    x_R = single((calibProj_Tsai(camParaCalib(icam),r_R))); % FLOAT
    x_L = single((calibProj_Tsai(camParaCalib(icam),r_L)));
    center = (calibProj_Tsai(camParaCalib(icam),cg));
    
    % First we organize the experimental data
    
    pixd = single(data2d(icam).ind(:,1:2));
    vd = ones(size(data2d(icam).ind(:,3)));
    
    imgd=zeros(1024,1280,'single');
	imgm=zeros(1024,1280,'single');
    imgm_t=zeros(1024,1280,'single');
    
%     Create the box of pixels around the center.
    xcuts = [1280:-1:(round(center(1)+93.75)),round(center(1)-93.75):-1:1];
    ycuts = [1024:-1:round(center(2)+75),round(center(2)-75):-1:1];
    xcuts(xcuts<0) = [];
    xcuts(xcuts==0) = [];
    xcuts(xcuts>1280) = [];
    ycuts(ycuts<0) = [];
    ycuts(ycuts>1024) = [];
    ycuts(ycuts==0) = [];
    
    xpix = camParaCalib(icam).Npixw;
    ypix = camParaCalib(icam).Npixh;    
    pos=zeros(size(vd));
    
    pos(:)=sub2ind([ypix,xpix],pixd(:,2),pixd(:,1));
    imgd(pos)=vd;
    
    for x = xcuts
            imgd(:,x) = [];
    end
    for y = ycuts
            imgd(y,:) = [];
    end


%     Now create the list of the pixels within that box to be put into
%     gaussian_intensity. For this, we need the Cartesian product.
%     
      xvals = round(center(1)-93.75):round(center(1)+93.75);
      xvals(xvals<0) = [];
      xvals(xvals==0) = [];
      xvals(xvals>1280) = [];
      
      yvals = round(center(2)-75):round(center(2)+75);
      yvals(yvals<0) = [];
      yvals(yvals==0) = [];
      yvals(yvals>1024) = [];
      
      [pixx,pixy] = meshgrid(xvals,yvals);
      pixd_new = [pixx(:),pixy(:)];

    I0 = single(1);
    % Compute the intensity of the model everywhere in the box.
    vm = bc_sk_gaussian_intensity_multiple_rods(pixd_new,x_L,x_R,I0,d);     

    for pix = 1:length(vm)
        imgm(pixd_new(pix,2),pixd_new(pix,1)) = vm(pix);
    end
    for x = xcuts
            imgm(:,x) = [];
    end
    for y = ycuts
            imgm(y,:) = [];
    end
    
    vm_t = vm;
    vm_t(vm_t<.6) = 0;
    vm_t(vm_t>=.6) = 1;
       
    for pix = 1:length(vm_t)
        imgm_t(pixd_new(pix,2),pixd_new(pix,1)) = vm_t(pix);
    end
    for x = xcuts
            imgm_t(:,x) = [];
    end
    for y = ycuts
            imgm_t(y,:) = [];
    end
%     vm = sk_model_intensity(pixd(:,1),pixd(:,2),x_L(:,1),x_L(:,2),x_R(:,1),x_R(:,2),I0,d);     


% ---------- Plotting-----------
%  There are many plotting options, each of which displays the data in a
%  slightly different way.
    
    figure(23);subplot(2,2,icam),imagesc(imgm(:,:))

%     The following lines are various other ways of viewing what this
%     program is doing. The one chosen above is chosen to display the model
%     that is created.
%     figure(10);subplot(2,2,icam),imagesc(imgd(:,:))
%     figure(24);subplot(2,2,icam),imagesc((imgm(:,:)-imgd(:,:)))
%     figure(25);subplot(2,2,icam),imagesc((imgm_t(:,:)-imgd(:,:)))
%     figure(26);subplot(2,2,icam),imagesc((imgm_t(:,:)))
    hold on
%     Ways of viewing the model particle on top of the camera images.
%     new_center = [93.75,75];
%     x_Rp = repmat(new_center,4,1);
%     x_Lp = x_L - x_R + x_Rp;
%     xset = [x_Lp(1,1);x_Rp(1,1);x_Lp(2,1);x_Rp(1,1);x_Lp(3,1);x_Rp(1,1);x_Lp(4,1)];
%     yset = [x_Lp(1,2);x_Rp(1,2);x_Lp(2,2);x_Rp(1,2);x_Lp(3,2);x_Rp(1,2);x_Lp(4,2)];
%     plot([x_L(:,1);x_R(end,1)],[x_L(:,2);x_R(end,2)],'-k','LineWidth',4);
%     plot([x_L(:,1);x_R(end,1)],[x_L(:,2);x_R(end,2)],'-wo','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',7);
%     plot(xset,yset,'-k','LineWidth',4);
%     plot(x_Rp(1,1),x_Rp(1,2),'wo','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',7);
%     plot(x_Lp(:,1),x_Lp(:,2),'wo','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',7);
    hold off
    title(sprintf('Camera %d',icam));
    set(gca, 'XTickLabelMode', 'Manual')
    set(gca, 'XTick', [])
    set(gca, 'YTickLabelMode', 'Manual')
    set(gca, 'YTick', [])
%     
    pause(0.01);
    
end

end
