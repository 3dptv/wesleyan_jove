function model = bc_object
% This function makes the tetrad that can be taken into other functions
% and rotated.
 
       % unit length arms of model tetrad       
       arm1 = [2*sqrt(2)/3,0,1/3];
       arm2 = [-sqrt(2)/3,-sqrt(2/3),1/3];
       arm3 = [-sqrt(2)/3,sqrt(2/3),1/3];
       arm4 = [0,0,-1];  
       arms = [arm1;arm2;arm3;arm4];
       model.r_L = arms;
       model.r_R = zeros(size(arms,1),3);

end