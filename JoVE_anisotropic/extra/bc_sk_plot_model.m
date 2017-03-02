function bc_sk_plot_model(eul_cntr,data2d,camParaCalib,param,model)
% This function is organized just like the leastSqOriPos functions, but it
% has the plotting portions uncommented.

ncams = length(camParaCalib);
cg = eul_cntr(4:6);
d = param.diamet;
hl = param.armlength;
imgd=zeros(ncams,1024,1280,'single');
imgm=zeros(ncams,1024,1280,'single');

%make the model
eul = eul_cntr(1:3);
ori = gg_ori(eul);

r_L = zeros(length(model.r_L),3);
r_R = zeros(length(model.r_R),3);

for i=1:length(model.r_L)
    r_L(i,:) = (ori*model.r_L(i,:)')';
    r_R(i,:) = (ori*model.r_R(i,:)')';
end
r_L = bsxfun(@plus, r_L*hl,cg);
r_R = bsxfun(@plus, r_R*hl,cg);

for icam=1:ncams
    xpix = camParaCalib(icam).Npixw;
    ypix = camParaCalib(icam).Npixh;
    % Determine the end points of the objects projected onto each camera
    x_R = single((calibProj_Tsai(camParaCalib(icam),r_R))); % FLOAT
    x_L = single((calibProj_Tsai(camParaCalib(icam),r_L)));

    xset = [x_L(1,1);x_R(1,1);x_L(2,1);x_R(1,1);x_L(3,1);x_R(1,1);x_L(4,1)];
    yset = [x_L(1,2);x_R(1,2);x_L(2,2);x_R(1,2);x_L(3,2);x_R(1,2);x_L(4,2)];
    
    % First we organize the experimental data
    
    pixd = single(data2d(icam).ind(:,1:2));
    vd = data2d(icam).ind(:,3);
    
    % We use the average intensity on each image to set the maximum value
    % in our model.
    I0 = single(mean(vd)); % FLOAT
    % Compute the intensity of the model at the bright points in the data
    vm = bc_sk_gaussian_intensity_multiple_rods(pixd,x_L,x_R,I0,d);
    pos=zeros(size(vm));

    pos(:)=sub2ind([ypix,xpix],pixd(:,2),pixd(:,1));
    imgd(icam,pos)=vd;
    imgm(icam,pos)=vm;
    
    center = (calibProj_Tsai(camParaCalib(icam),cg));
    figure(5)
    subplot(2,2,icam),imagesc(squeeze(imgd(icam,:,:))); % change imgm -> imgd to see image data.
    hold on
%     The following lines are all different ways to see the data and the
%     model.
%     colormap(gray)
%     scatter(pixd(:,1),pixd(:,2),1,vm);
%     plot([x_L(:,1);x_R(end,1)],[x_L(:,2);x_R(end,2)],'-k','LineWidth',4);
    plot(xset,yset,'-k','LineWidth',4);
    plot(x_R(1,1),x_R(1,2),'wo','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',7);
    plot(x_L(:,1),x_L(:,2),'wo','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',7);

%     axis equal
    axis([center(1)-93.75,center(1)+93.75,center(2)-75,center(2)+75]);
    hold off;
    
end
% A way of getting the MATLAB figure to update if this is run over multiple
% frames.
pause(0.001);

end
