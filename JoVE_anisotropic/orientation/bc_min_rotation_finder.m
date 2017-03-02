function data_out = bc_min_rotation_finder(data_opt)

% This program is step 5.1.6 of the protocol: checking to make sure that
% equivalent symmetric particle orientations are used in successive frames.

data_out = data_opt;

% Tetrads --------------------------------------------------------------

% Building the rotation matrix that keeps arm3 constant while rotating the
% rest of the tetrad around it by 120 degrees. So R3*arm4 = arm1, R3*arm1 =
% arm2, R3*arm2 = arm4, R3*arm3 = arm3.

phi3 = -pi/2;
theta3 = pi/2 + asin(1/3);

% Solving analytically for psi leads to a quadratic equation with these
% solutions:
psi3 = pi/6;
psi3a = 5*pi/6;

R3 = gg_ori([phi3, theta3, psi3]);

% Rotation about arm4 is obvious, since phi is rotation about the z axis
% and arm4 = [0;0;-1]. The phi that accomplishes this is phi = -2*pi/3. phi
% = 2*pi/3 would also work, but the arms are defined clockwise from one
% another and phi = -2*pi/3 brings us from arm1->arm2->arm3.
phi4 = -2*pi/3;
theta4 = 0;
psi4 = 0;

R4 = gg_ori([phi4, theta4, psi4]);

% From these two insights, we can construct the other two rotation
% matrices. Because phi is rotation about z and determines which arm will
% be kept stationary, adding multiples of 2*pi/3 to phi3 gives us the other
% two rotation matrices. So, to keep arm1 steady, we need to go 2*pi/3 away.
% -pi/2 + 2*pi/3 = pi/6.
phi1 = pi/6;
theta1 = theta3;
psi1 = psi3;

R1 = gg_ori([phi1, theta1, psi1]);

% Rotation around arm2 needs to be done analytically. The only 
% change from arm3 here is that pi/2 is added to theta (asin(1/3)+pi/2=acos(-1/3))

phi2 = -pi/2;
theta2 = acos(-1/3);
% psi3a works here where psi3 doesn't.
psi2 = psi3a;

R2 = gg_ori([phi2, theta2, psi2]);

% Now we create one matrix with all of the possible equivalent
% orientations. Checking bc_degenerateorientations_plots confirms that this
% is the correct set of rotation matrices.

fit_mat = zeros(3,3,12);

fit_mat(:,:,1) = R1;
fit_mat(:,:,2) = R1*R1*R4;
fit_mat(:,:,3) = R4*R1;

fit_mat(:,:,4) = R1*R4;
fit_mat(:,:,5) = R1*R1*R4*R4;
fit_mat(:,:,6) = R4*R3;

fit_mat(:,:,7) = R1*R4*R4;
fit_mat(:,:,8) = R4*R3*R3;
fit_mat(:,:,9) = R4*R4*R1*R1;

fit_mat(:,:,10) = R4;
fit_mat(:,:,11) = R4*R4;
fit_mat(:,:,12) = [1,0,0;0,1,0;0,0,1];

% The rest of this code tracks through the Euler angles found in the
% orientation analysis and checks for times where different degenerate
% orientations were found for sequential frames. Such frames are smoothed
% back to the degenerate orientation used for the previous frame.

frames = data_opt(:,25);
N = length(frames);

eul = data_opt(:,5:7);
ori = gg_ori(eul(1,:));

match = zeros(size(fit_mat,3),1);

for i=0:N-1

    orin = gg_ori(eul(i+1,:));
    for j = 1:size(fit_mat,3)
        orij = orin*fit_mat(:,:,j)';
        match(j) = norm(ori-orij);
    end

%  Find the rotation matrix with the smallest norm and make that the new
%  matrix.

    [~,k] = min(match);

    orif = orin*fit_mat(:,:,k)';
    ori = orif;

    data_out(i+1,5:7) = bc_euler(orif);

end

end