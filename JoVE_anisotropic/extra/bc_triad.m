function model = bc_triad

arm1 = [0,1,0];
arm2 = [cos(7*pi/6),sin(7*pi/6),0];
arm3 = [cos(-pi/6),sin(-pi/6),0];
arms = [arm1;arm2;arm3];
model.r_L = arms;
model.r_R = zeros(size(arms,1),3);

end