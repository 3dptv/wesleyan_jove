function model = bc_jack

arm1 = [0,.5,0];
arm2 = [0,-.5,0];
arm3 = [.5,0,0];
arm4 = [-.5,0,0];
arm5 = [0,0,.5];
arm6 = [0,0,-.5];

arms = [arm1;arm2;arm3;arm4;arm5;arm6];
model.r_L = arms;
model.r_R = zeros(size(arms,1),3);

end