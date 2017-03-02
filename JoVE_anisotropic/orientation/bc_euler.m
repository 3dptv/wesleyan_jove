function [eul1,eul2] = bc_euler(o)
% Find Euler angles from a rotation matrix using A from the Goldstein
% convention. Accounted for the degeneracy in the Euler angles and
% found a proper triplet.


% amount residual at Euler angle pole
rp=abs(o(3,3)-1);
rm=abs(o(3,3)+1);

% Boolean saying if the Euler angle is too close to pole
boolp = rp<1e-6;
boolm = rm<1e-6;

if ( boolp || boolm )
    ps1 = 0; % Gimbal lock, psi is undefined
    if boolp
        th1 = rp;
        ph1 = atan2(o(1,2),o(1,1)); % If we are rotating the coordinate system
%         ph1 = atan2(o(1,2),-o(1,1)); % If we rotate the object.
    else
        th1 = pi-rm; % same as -pi in our case
        ph1 = atan2(o(1,2),o(1,1)); % If we are rotating the coordinate 
                                    % system or rotating the object. In
                                    % this case, they have the same
                                    % condition.
    end
    
    % There is only one set of angles in this case
    ph2=ph1;
    th2=th1;
    ps2=ps1;
    
else
    % There are two possible Euler representations except in the case
    % handled above.
    
    % first set of angles
%     th1 = acos(o(3,3));
    th1 = real(acos(o(3,3)));
    ph1 = atan2(o(3,1)/sin(th1),-o(3,2)/sin(th1));
    ps1 = atan2(o(1,3)/sin(th1),o(2,3)/sin(th1));
    
    % second set of angles
    th2 = -th1;
    ph2 = atan2(o(3,1)/sin(th2),-o(3,2)/sin(th2));
    ps2 = atan2(o(1,3)/sin(th2),o(2,3)/sin(th2));
    
end

    % Return the angles
    eul1 = [ph1,th1,ps1];
    eul2 = [ph2,th2,ps2];

end