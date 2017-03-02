function model = bc_cross

arm1 = [0,.5,0];
arm2 = [.5,0,0];
arm3 = [0,-.5,0];
arm4 = [-.5,0,0];

arms = [arm1;arm2;arm3;arm4];
model.r_L = arms;
model.r_R = zeros(size(arms,1),3);

end