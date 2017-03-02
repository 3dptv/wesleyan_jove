% This produces figures showing all of the equivalent orientations for each
% type of particle discussed in this paper.

% Mac computers address files with the slashes going in the opposite
% direction of the slashes on Windows computers.
mac = 0; % =0 if not on a Mac, =1 if on a Mac.
if mac
    slash = '/';
elseif ~mac
    slash = '\';
end

% The following line will need to be changed to fit the file structure on
% your computer.
pathlocation = strcat('Z:',slash,'bcole',slash,'JoVE',slash,'MATLABcodes');

% Find the correct paths for the programs that will be called.
addpath(strcat(pathlocation, slash, 'orientation', slash))
addpath(strcat(pathlocation, slash, 'extra', slash))


% Begin with all of the rotation matrices from bc_min_rotation_finder that
% produce equivalent tetrad orientations.


% -------------------Tetrad Degenerate Orientations----------
model = bc_object;


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
% orientations.

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


% ------End Tetrad Matrices---------

% Plot tetrads
eul = zeros(1,3);
eul(2) = eul(2) + pi;

r_L = zeros(length(model.r_L),3);
r_R = zeros(length(model.r_R),3);

figure(1);hold on
subplot(4,3,1:12)

for i=1:12

    ori = gg_ori(eul);
    ori = ori*fit_mat(:,:,i);
    for j=1:length(model.r_L)
        r_L(j,:) = (ori*model.r_L(j,:)')';
        r_R(j,:) = (ori*model.r_R(j,:)')';
    end
    
    xset = [r_L(1,1);0;r_L(2,1);0;r_L(3,1);0;r_L(4,1);0];
    yset = [r_L(1,2);0;r_L(2,2);0;r_L(3,2);0;r_L(4,2);0];
    zset = [r_L(1,3);0;r_L(2,3);0;r_L(3,3);0;r_L(4,3);0];

    subplot(4,3,i);hold on
    
    plot3(xset(1:2),yset(1:2),zset(1:2),'-b','LineWidth',4);
    plot3(xset(3:4),yset(3:4),zset(3:4),'-r','LineWidth',4);
    plot3(xset(5:6),yset(5:6),zset(5:6),'-k','LineWidth',4);
    plot3(xset(7:8),yset(7:8),zset(7:8),'-g','LineWidth',4);
    axis off
    view(100,15);
    axis equal
    hold off
    
end


% -----End Tetrads--------

model = bc_triad;

% -----Triad Degenerate Orientations------

% Triads are easier than tetrads because they are planar. The first three
% degenerate cases are formed by changing just phi by multiples of 2*pi/3.
% Then the triad is flipped over by setting theta = pi and reset by setting
% psi = pi. Then phi = multiples of 2*pi/3 again.

fit_mat = zeros(3,3,6);
fit_mat(:,:,1) = [1,0,0;0,1,0;0,0,1];
fit_mat(:,:,2) = gg_ori([2*pi/3,0,0]);
fit_mat(:,:,3) = gg_ori([4*pi/3,0,0]);
fit_mat(:,:,4) = gg_ori([0,pi,pi]);
fit_mat(:,:,5) = gg_ori([2*pi/3,pi,pi]);
fit_mat(:,:,6) = gg_ori([4*pi/3,pi,pi]);


% Plot triads
eul = zeros(1,3);

r_L = zeros(length(model.r_L),3);
r_R = zeros(length(model.r_R),3);

idx = reshape(1:6,2,[])';
idx = idx(:);

for i=1:6

    ori = gg_ori(eul);
    ori = ori*fit_mat(:,:,i);
    for j=1:length(model.r_L)
        r_L(j,:) = (ori*model.r_L(j,:)')';
        r_R(j,:) = (ori*model.r_R(j,:)')';
    end
    
    xset = [r_L(1,1);0;r_L(2,1);0;r_L(3,1);0];
    yset = [r_L(1,2);0;r_L(2,2);0;r_L(3,2);0];
    
    figure(2)
    subplot(3,2,idx(i));hold on
    
    plot(xset(1:2),yset(1:2),'-b','LineWidth',4);
    plot(xset(3:4),yset(3:4),'-r','LineWidth',4);
    plot(xset(5:6),yset(5:6),'-k','LineWidth',4);
    axis off
    axis equal
    hold off
    
end

% ------End Triads-------


% ------Cross Degenerate Orienations--------

% Finding the cross degenerate orientations is much the same process
% as for the triads, except that phi changes by multiples of pi/2 each
% time, not 2*pi/3.

fit_mat = zeros(3,3,8);

for i = 1:4
    fit_mat(:,:,i) = gg_ori([i*pi/2,0,0]);
end

% Because of the symmetry of crosses, setting psi = pi here isn't strictly
% necessary to get an equivalent object. But it is visually useful to start
% both rows with the same arm on top of the image.

for i = 1:4
    fit_mat(:,:,i+4) = gg_ori([i*pi/2,pi,pi]);
end

model = bc_cross;

% Plot crosses
eul = zeros(1,3);
eul(2) = eul(2) + pi;

r_L = zeros(length(model.r_L),3);
r_R = zeros(length(model.r_R),3);

idx = reshape(1:8,4,[])';
idx = idx(:);

figure(3);hold on
subplot(4,2,1:8)

for i=1:8

    ori = gg_ori(eul);
    ori = ori*fit_mat(:,:,i);
    for j=1:length(model.r_L)
        r_L(j,:) = (ori*model.r_L(j,:)')';
        r_R(j,:) = (ori*model.r_R(j,:)')';
    end
    
    xset = [r_L(1,1);0;r_L(2,1);0;r_L(3,1);0;r_L(4,1);0];
    yset = [r_L(1,2);0;r_L(2,2);0;r_L(3,2);0;r_L(4,2);0];

    figure(3)
    subplot(4,2,idx(i));hold on
    
    plot(xset(1:2),yset(1:2),'-b','LineWidth',4);
    plot(xset(3:4),yset(3:4),'-r','LineWidth',4);
    plot(xset(5:6),yset(5:6),'-k','LineWidth',4);
    plot(xset(7:8),yset(7:8),'-g','LineWidth',4);
    axis off
    axis equal tight
    hold off
    
end

% --------End Crosses--------

% -------Jack Degenerate Orientations-------

% Jacks have more than double the number of degenerate orientations of any
% of the other objects, so they will be split into two figures.
% Nevertheless, because the arms of a jack are along the coordinate axes,
% their equivalent rotations are straightforward.

fit_mat=zeros(24,3,3);
Ea=zeros(24,3);
cnt=1;
phi=0;
theta=0;
for psi=0:pi/2:3*pi/2
    Ea(cnt,:)=[phi,theta,psi];
D = [cos(phi), sin(phi), 0; -sin(phi), cos(phi), 0; 0, 0, 1];
C = [1, 0, 0; 0, cos(theta), sin(theta); 0, -sin(theta), cos(theta)];
B = [cos(psi), sin(psi), 0; -sin(psi), cos(psi), 0; 0, 0, 1];

fit_mat(cnt,:,:) = B * C * D;
cnt=cnt+1;
end

phi=0;
theta=pi/2;
for psi=0:pi/2:3*pi/2
        Ea(cnt,:)=[phi,theta,psi];
D = [cos(phi), sin(phi), 0; -sin(phi), cos(phi), 0; 0, 0, 1];
C = [1, 0, 0; 0, cos(theta), sin(theta); 0, -sin(theta), cos(theta)];
B = [cos(psi), sin(psi), 0; -sin(psi), cos(psi), 0; 0, 0, 1];

fit_mat(cnt,:,:) = B * C * D;
cnt=cnt+1;
end

phi=0;
theta=pi;
for psi=0:pi/2:3*pi/2
        Ea(cnt,:)=[phi,theta,psi];
D = [cos(phi), sin(phi), 0; -sin(phi), cos(phi), 0; 0, 0, 1];
C = [1, 0, 0; 0, cos(theta), sin(theta); 0, -sin(theta), cos(theta)];
fit_mat(cnt+24,:,:) = B * C * D;
fit_mat(cnt+24,:,1) = -fit_mat(cnt+24,:,1);
B = [cos(psi), sin(psi), 0; -sin(psi), cos(psi), 0; 0, 0, 1];

fit_mat(cnt,:,:) = B * C * D;
cnt=cnt+1;
end

phi=pi/2;
theta=pi/2;
for psi=0:pi/2:3*pi/2
        Ea(cnt,:)=[phi,theta,psi];
D = [cos(phi), sin(phi), 0; -sin(phi), cos(phi), 0; 0, 0, 1];
C = [1, 0, 0; 0, cos(theta), sin(theta); 0, -sin(theta), cos(theta)];
B = [cos(psi), sin(psi), 0; -sin(psi), cos(psi), 0; 0, 0, 1];

fit_mat(cnt,:,:) = B * C * D;
cnt=cnt+1;
end

phi=pi;
theta=pi/2;
for psi=0:pi/2:3*pi/2
        Ea(cnt,:)=[phi,theta,psi];
D = [cos(phi), sin(phi), 0; -sin(phi), cos(phi), 0; 0, 0, 1];
C = [1, 0, 0; 0, cos(theta), sin(theta); 0, -sin(theta), cos(theta)];
B = [cos(psi), sin(psi), 0; -sin(psi), cos(psi), 0; 0, 0, 1];

fit_mat(cnt,:,:) = B * C * D;
cnt=cnt+1;
end

phi=3*pi/2;
theta=pi/2;
for psi=0:pi/2:3*pi/2
        Ea(cnt,:)=[phi,theta,psi];
D = [cos(phi), sin(phi), 0; -sin(phi), cos(phi), 0; 0, 0, 1];
C = [1, 0, 0; 0, cos(theta), sin(theta); 0, -sin(theta), cos(theta)];
B = [cos(psi), sin(psi), 0; -sin(psi), cos(psi), 0; 0, 0, 1];

fit_mat(cnt,:,:) = B * C * D;
cnt=cnt+1;
end

% ------End Jack Matrices---------

% Plot jacks
model = bc_jack;

eul = zeros(1,3);

r_L = zeros(length(model.r_L),3);
r_R = zeros(length(model.r_R),3);

idx = reshape(1:12,3,[])';
idx = idx(:);

for i=1:12

    ori = gg_ori(eul);
    ori = ori*squeeze(fit_mat(i,:,:));
    for j=1:length(model.r_L)
        r_L(j,:) = (ori*model.r_L(j,:)')';
        r_R(j,:) = (ori*model.r_R(j,:)')';
    end
    
    xset = [r_L(1,1);0;r_L(2,1);0;r_L(3,1);0;r_L(4,1);0;r_L(5,1);0;r_L(6,1);0];
    yset = [r_L(1,2);0;r_L(2,2);0;r_L(3,2);0;r_L(4,2);0;r_L(5,2);0;r_L(6,2);0];
    zset = [r_L(1,3);0;r_L(2,3);0;r_L(3,3);0;r_L(4,3);0;r_L(5,3);0;r_L(6,3);0];

    figure(4)
    subplot(4,3,idx(i));hold on
    
    plot3(xset(1:2),yset(1:2),zset(1:2),'-b','LineWidth',4);
    plot3(xset(3:4),yset(3:4),zset(3:4),'-r','LineWidth',4);
    plot3(xset(5:6),yset(5:6),zset(5:6),'-k','LineWidth',4);
    plot3(xset(7:8),yset(7:8),zset(7:8),'-g','LineWidth',4);
    plot3(xset(9:10),yset(9:10),zset(9:10),'-m','LineWidth',4);
    plot3(xset(11:12),yset(11:12),zset(11:12),'-c','LineWidth',4);
    axis off
    view(46,15);
    axis equal
    hold off
    
end

for i=13:24

    ori = gg_ori(eul);
    ori = ori*squeeze(fit_mat(i,:,:));
    for j=1:length(model.r_L)
        r_L(j,:) = (ori*model.r_L(j,:)')';
        r_R(j,:) = (ori*model.r_R(j,:)')';
    end
    
    xset = [r_L(1,1);0;r_L(2,1);0;r_L(3,1);0;r_L(4,1);0;r_L(5,1);0;r_L(6,1);0];
    yset = [r_L(1,2);0;r_L(2,2);0;r_L(3,2);0;r_L(4,2);0;r_L(5,2);0;r_L(6,2);0];
    zset = [r_L(1,3);0;r_L(2,3);0;r_L(3,3);0;r_L(4,3);0;r_L(5,3);0;r_L(6,3);0];

    figure(5)
    subplot(4,3,idx(i-12));hold on
    
    plot3(xset(1:2),yset(1:2),zset(1:2),'-b','LineWidth',4);
    plot3(xset(3:4),yset(3:4),zset(3:4),'-r','LineWidth',4);
    plot3(xset(5:6),yset(5:6),zset(5:6),'-k','LineWidth',4);
    plot3(xset(7:8),yset(7:8),zset(7:8),'-g','LineWidth',4);
    plot3(xset(9:10),yset(9:10),zset(9:10),'-m','LineWidth',4);
    plot3(xset(11:12),yset(11:12),zset(11:12),'-c','LineWidth',4);
    axis off
    view(46,15);
    axis equal
    hold off
    
end
