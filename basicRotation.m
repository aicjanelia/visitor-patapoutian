function R = basicRotation(theta,axis)

if strcmp(axis,'x') || strcmp(axis,'X')
    R = [1 0 0 ; 0 cos(theta) -sin(theta) ; 0 sin(theta) cos(theta)];
elseif strcmp(axis,'y') || strcmp(axis,'Y')
    R = [cos(theta) 0 sin(theta) ; 0 1 0 ; -sin(theta) 0 cos(theta)];
elseif strcmp(axis,'z') || strcmp(axis,'Z')
    R = [cos(theta) -sin(theta) 0 ; sin(theta) cos(theta) 0 ; 0 0 1];
else
    error('Specify axis as x, y, or z')
end