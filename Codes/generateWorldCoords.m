% Automatically generate the world coordinates
world_coords = zeros(160, 4);
count = 1;

% Generate coordinates of the points at the left side of the checkboard
init_left_point = [0, 10, 0, 1];
for i = 1:1:8
    % bottom up method
    point_coord = init_left_point;
    point_coord(3) = init_left_point(3) - 28;
    for j = 1:1:10
        % add z
        point_coord = [init_left_point(1), init_left_point(2), point_coord(3) + 28, 1];
        world_coords(count, :) = point_coord;
        count = count + 1;
    end
    % add y
    init_left_point = [init_left_point(1), init_left_point(2) + 28, init_left_point(3), 1];
end

% Generate coordinates of the points at the right side of the checkboard
init_right_point = [11, 0, 0, 1];
for p = 1:1:8
    % bottom up method
    point_coord = init_right_point;
    point_coord(3) = init_right_point(3) - 28;
    for q = 1:1:10
        % add z
        point_coord = [init_right_point(1), init_right_point(2), point_coord(3) + 28, 1];
        world_coords(count, :) = point_coord;
        count = count + 1;
    end
    % add x
    init_right_point = [init_right_point(1) + 28, init_right_point(2), init_right_point(3), 1];
end

disp(world_coords);