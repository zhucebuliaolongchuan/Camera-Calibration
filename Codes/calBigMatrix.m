P = zeros(12,12);
world = [67,0,56,1;179,0,112,1;95,0,196,1;0,122,224,1;0,150,140,1;0,66,112,1];
% world = world_coords;
camera = [1849.6,1975.9;2552.6,1663.7;2024.0,964.6;833.8,727.2;639.6,1380.8;1188.0,1566.2];
% camera = camera_coords;

disp('P=');
disp(P);

for p_row = 1:1:12
    if mod(p_row,2) == 1
        world_row = int32((p_row + 1) / 2);
        for col = 1:1:4
            P(p_row,col) = world(world_row,col);
            P(p_row,col+8) = world(world_row,col) * camera(world_row,1) * (-1);
        end
    else
        world_row = int32(p_row / 2);
        for col = 5:1:8
            P(p_row,col) = world(world_row,col-4);
            P(p_row,col+4) = world(world_row,col-4) * camera(world_row,2) * (-1);
        end
    end
end

% Perform SVD of P
[U,S,V] = svd(P);
[min_val, min_index] = min(diag(S(1:12,1:12)));
% m is given by right singular vector of minimum singular value
m = V(1:12, min_index);
disp(m);
M = [m(1),m(2),m(3),m(4);m(5),m(6),m(7),m(8);m(9),m(10),m(11),m(12)];
disp('M =');
disp(M);

% The origin of the world coordinate is in front of the camera
p = 1 / norm(M(3,1:3));
disp('p =');
disp(p);

% Cal U0  
U0 = dot(M(1,1:3), M(3,1:3)) * (p ^ 2);
disp('U0 =')
disp(U0)

% Cal V0
V0 = dot(M(2,1:3), M(3,1:3)) * (p ^ 2);
disp('V0 =')
disp(V0)

% Cal Theta
a = cross(M(1,1:3), M(3,1:3));
b = cross(M(2,1:3), M(3,1:3));
cos_theta = - dot(a,b) / norm(a) * norm(b);
disp('cos_theta='); 
disp(cos_theta)
theta_in_radian = acos(cos_theta);
theta_in_degree = acosd(cos_theta);
disp('theta in radians=');
disp(theta_in_radian);
disp('theta in degree=');
disp(theta_in_degree);
sin_theta = sin(theta_in_radian);
disp('sin_theta=');
disp(sin_theta);

% Cal alpha
alpha = (p ^ 2) * norm(cross(M(1,1:3),M(3,1:3))) * sin(theta_in_radian);
disp('alpha=');
disp(alpha);

% Cal beta
beta = (p ^ 2) * norm(cross(M(2,1:3),M(3,1:3))) * sin(theta_in_radian);
disp('beta=');
disp(beta);

% Cal R1
R1 = cross(M(2,1:3),M(3,1:3)) / norm(cross(M(2,1:3),M(3,1:3)));
disp('R1=');
disp(R1);

% Cal R3
R3 = p * M(3,1:3);
disp('R3=');
disp(R3);

% Cal R2
R2 = cross(R3, R1);
disp('R2=');
disp(R2);

% Cal Rotation in Euler
rotate_matrix = [R1;R2;R3];
rotate_euler = rotm2eul(rotate_matrix);
disp('rotate_euler=');
disp(rotate_euler);
rotate_degree = radtodeg(rotate_euler);
disp('rotate_degree=');
disp(rotate_degree);

% Cal Translation Vector
K = [alpha, -alpha * cot(theta_in_radian), U0; 0,beta / sin(theta_in_radian), V0; 0, 0, 1];
disp('K=');
disp(K);
b = [M(1,4); M(2,4); M(3,4)];
t = p * inv(K) * b;
disp('t=');
disp(t);

% Validate the decomposition  p M = K(R t)
M_k = [R1, t(1); R2, t(2); R3, t(3)];
disp(p * M);
disp(K * M_k);
disp(p * M == K * M_k)

% Reconstruct: Use the calibrate matrix to reconstruct the world coordinate to pixel coordinate
reconstruct_cam_coord = zeros(160,2);
for point = 1:1:160
    world_coord = transpose(world_coords(point,:));
    scale = M(3,:) * world_coord;
    cam_coord = M * world_coord / scale;
    reconstruct_cam_coord(point,1) = cam_coord(1);
    reconstruct_cam_coord(point,2) = cam_coord(2);
end

figure(1);
img = imread('Image1.jpeg');
% img = imread('bonus_1.jpeg');
image(img);
hold on;

% show the selected pointed by being marked by green circle
scatter(camera(:,1), camera(:,2), 80, 'g');
% show the reconstruct points at the matrix via the estimate calibration matrix
scatter(reconstruct_cam_coord(:,1), reconstruct_cam_coord(:,2), 50, 'r.');