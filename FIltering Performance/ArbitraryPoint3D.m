function [meas_dist_vec,points_3d] = ArbitraryPoint3D(noise_Var)
    x_points = 630*rand(100, 1);
    y_points = 540*rand(100, 1);
    t_points = linspace(1,100);
    
    x_points(1,1) = 0;
    y_points(1,1) = 0;
    
    x_interp = interp1(t_points, x_points, linspace(1, 5, 100), 'spline');
    y_interp = interp1(t_points, y_points, linspace(1, 5, 100), 'spline');
    
    points_3d = [x_interp', y_interp'];
    
    %disp(points_3d);
    
    ANK_1 = [0,0];
    ANK_2 = [630,0];
    ANK_3 = [0,540];
    
    true_dist_vec = zeros(100,3);
    meas_dist_vec = zeros(100,3);
    noise_Std = sqrt(noise_Var);
    
    for ind = 1:1:100
        true_dist_vec(ind,1) = sqrt(...
            (points_3d(ind,1) - ANK_1(1,1))^2+...
            (points_3d(ind,2) - ANK_1(1,2))^2);
        meas_dist_vec(ind,1) = true_dist_vec(ind,1) + noise_Std*randn;
    
        true_dist_vec(ind,2) = sqrt(...
            (points_3d(ind,1) - ANK_2(1,1))^2+...
            (points_3d(ind,2) - ANK_2(1,2))^2);
        meas_dist_vec(ind,2) = true_dist_vec(ind,2) + noise_Std*randn;
    
        true_dist_vec(ind,3) = sqrt(...
            (points_3d(ind,1) - ANK_3(1,1))^2+...
            (points_3d(ind,2) - ANK_3(1,2))^2);
        meas_dist_vec(ind,3) = true_dist_vec(ind,3) + noise_Std*randn;
    end
    %{
    % 3D plot
    plot3(points_3d(:,1), points_3d(:,2), points_3d(:,3), 'o-');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Interpolated 3D Points');
    grid on;
    axis equal;
    %}

    
end