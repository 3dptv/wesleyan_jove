function [eul_cntr,fval,fnew,out,xflag] = bc_sk_nonlinearopt(eul_cntr,data2d,camParaCalib,model)

eul = eul_cntr(1:3);
cg = eul_cntr(4:6);

% Parameters for our model.
param.diamet = 6;
param.armlength = 4;

param.TolX = 1e-6;
param.TolFun = 1e-6;

%     First, we optimize only the Euler angles, using the first guess of
%     the particle's initial position.
nonlin1 = tic;
options = optimset('TolX',param.TolX,'TolFun',param.TolFun,'Display','final','LargeScale','off'); 
[eul,fval,xflag,out] = fminsearch(@(e) bc_sk_leastSqOriPos1(e,data2d,camParaCalib,param,model,cg)...
    , eul,options);
display(sprintf('\tElapsed time for first fminsearch: %f sec', toc(nonlin1)));

%     Checking to make sure no Euler angles are outside of [-pi,pi]
for i = 1:3
       if eul(i) > pi
           eul(i) = eul(i) - 2*pi;
       elseif eul(i) < -pi
           eul(i) = eul(i) + 2*pi;
       end
end

eul_cntr(1:3) = eul;

%     Second, we optimize the 3D center and the Euler angles together.
nonlin2 = tic;
[eul_cntr,fval,xflag,out] = fminsearch(@(e) bc_sk_leastSqOriPos(e,data2d,camParaCalib,param,model)...
    , eul_cntr,options);
display(sprintf('\tElapsed time for second fminsearch: %f sec', toc(nonlin2)));

%     Checking again to make sure no Euler angles are outside of [-pi,pi]
eul = eul_cntr(1:3);
for i = 1:3
       if eul(i) > pi
           eul(i) = eul(i) - 2*pi;
       elseif eul(i) < -pi
           eul(i) = eul(i) + 2*pi;
       end
end
eul_cntr(1:3) = eul;

% Calculating two other measures of how well our model fits the data.
[fnew1,fnew2] = bc_new_residual(eul,data2d,camParaCalib,param,model,cg);
fnew = [fnew1,fnew2];

% Uncomment one of the following two lines to see the final matched model.
% bc_sk_plot_model(eul_cntr,data2d,camParaCalib,param,model);
% bc_new_residual_plot(eul,data2d,camParaCalib,param,model,cg);
end
