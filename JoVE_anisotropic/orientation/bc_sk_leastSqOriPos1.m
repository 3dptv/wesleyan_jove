function diffsq = bc_sk_leastSqOriPos1(eul,data2d,camParaCalib,param,model,cg)
% The first step in the nonlinear search. This function calculates
% the difference in intensity between the model and images on all four
% cameras.

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

diffsq = 0;
% To see a plot of some of the steps of the nonlinear search,
% uncomment the following five lines. There are hundreds of steps in the
% search; using a random number shows a small sample of those steps.
% a = randn(1);
% if round(abs(a)) == 2
%     imgd = zeros(ncams,1024,1280);
%     imgm = zeros(ncams,1024,1280);
% end

for icam=1:ncams
    % Determine the end points of the objects projected onto each camera
    x_R = single((calibProj_Tsai(camParaCalib(icam),r_R))); % FLOAT
    x_L = single((calibProj_Tsai(camParaCalib(icam),r_L)));

    % Experimental data is flattened to be all ones to reduce the effects
    % of uneven illumination of the particles.
    pixd = single(data2d(icam).ind(:,1:2));
    vd = ones(size(data2d(icam).ind(:,3)));
    
    I0 = single(1);
    % Compute the intensity of the model at the bright points in the data
    vm = bc_sk_gaussian_intensity_multiple_rods(pixd,x_L,x_R,I0,d);
    data2d(icam).indmod(:,3) = vm;
    
    % To view the progress of the nonlinear search, uncomment
    % the following lines.
%     if round(abs(a)) == 2
%         pos=zeros(size(vm));
%         xpix = camParaCalib(icam).Npixw;
%         ypix = camParaCalib(icam).Npixh;
%         xset = [x_L(1,1);x_R(1,1);x_L(2,1);x_R(1,1);x_L(3,1);x_R(1,1);x_L(4,1);x_R(1,1)];
%         yset = [x_L(1,2);x_R(1,2);x_L(2,2);x_R(1,2);x_L(3,2);x_R(1,2);x_L(4,2);x_R(1,2)];
%         pos(:)=sub2ind([ypix,xpix],pixd(:,2),pixd(:,1));
%         vd = data2d(icam).ind(:,3);
%         imgd(icam,pos)=vd;
%         imgm(icam,pos)=vm;
%         center = (calibProj_Tsai(camParaCalib(icam),cg));
%         figure(1)
%         subplot(2,2,icam),imagesc(squeeze(imgd(icam,:,:))); % ,colormap(gray)
%         hold on;
%         scatter(pixd(:,1),pixd(:,2),1,vd);

% %   Each arm of the model presented as a different color line.
%         plot(xset(1:2),yset(1:2),'-r','LineWidth',4)
%         plot(xset(3:4),yset(3:4),'-b','LineWidth',4)
%         plot(xset(5:6),yset(5:6),'-g','LineWidth',4)
%         plot(xset(7:8),yset(7:8),'-k','LineWidth',4)
        
% %   Circles at the ends of the arms to make them stand out.
%         plot(x_R(1,1),x_R(1,2),'wo','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',7);
%         plot(x_L(:,1),x_L(:,2),'wo','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',7);

%         axis([center(1)-93.75,center(1)+93.75,center(2)-75,center(2)+75]);
%         title(sprintf('Camera %d',icam));
%         hold off;
%     end

    
    % Now we compute the sum of the difference between the intesities
    % squared at only the points just computed.
    
    diffsq = diffsq + sum((vm-vd).^2);
    
end

end
