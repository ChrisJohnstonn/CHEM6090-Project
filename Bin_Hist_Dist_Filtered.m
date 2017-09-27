bins = 600;

c_xyz = vertcat(result_c1_xyz,result_c2_xyz);
c_xz_vel = zeros(length(c_xyz),1);
c_xyz_vel = zeros(length(c_xyz),1);
c_angle_dev = zeros(length(c_xyz),1);

result_filtered_simulations = zeros(length(c_xyz),2);

vel_min = 11000;
vel_max = 11500;

for i = 1:length(c_xyz)
    %determine angles
    c_angle_dev(i) = abs(atan2(c_xyz(i,1),c_xyz(i,3))*180/pi)-90;
    %determine velocity in x-z/xyz plane
    c_xz_vel(i) = sqrt(c_xyz(i,1)^2 + c_xyz(i,2)^2);
    c_xyz_vel(i) = sqrt(c_xyz(i,1)^2 + c_xyz(i,2)^2 + c_xyz(i,3)^2);
    %c_xz_vel(i) = abs(c_xyz(i,1));
    
    if c_xz_vel(i) > vel_min
        if  c_xz_vel(i) < vel_max
            if i <= length(c_xyz)/2
                result_filtered_simulations(i,1) = ceil(i/Rotation_Steps);
                result_filtered_simulations(i,2) = mod(i,300) - 1;
                
                if result_filtered_simulations(i,2) == -1
                    result_filtered_simulations(i,2) = 299;
                end
                
          
            else
                result_filtered_simulations(i,1) = ceil((i/Rotation_Steps)) - Simulations_Amount;
                result_filtered_simulations(i,2) = mod(i,300)-1;
                
                if result_filtered_simulations(i,2) == -1
                    result_filtered_simulations(i,2) = 299;
                end
            end
        end
    end
end

c_vel_max = max(c_xz_vel);
c_vel_min = min(c_xz_vel);
c_vel_diff = 1.2*c_vel_max - 0.8*c_vel_min;
c_vel_steps = c_vel_diff / bins;

x_vel_axis = (0.8*c_vel_min:c_vel_steps:1.2*c_vel_max);

result_dist_matrix = zeros(length(x_vel_axis),2);
result_dist_matrix(:,1) = x_vel_axis;

for i = 1:length(c_xz_vel)
    vel_x_bin =  ( c_xz_vel(i) - x_vel_axis(1) ) / c_vel_steps ;
    
    l = floor(vel_x_bin + 1);
    u = ceil(vel_x_bin + 1);
    b1 = x_vel_axis(l);
    b2 = x_vel_axis(u);
    l_p = ((b2)-(c_xz_vel(i)))/(b2-b1);
    u_p = 1 - l_p;
    
    if i <= length(c_xyz)/2
        result_dist_matrix(l,2) = result_dist_matrix(l,2) + l_p*result_intensity(ceil(i/Rotation_Steps),1);
        result_dist_matrix(u,2) = result_dist_matrix(u,2) + u_p*result_intensity(ceil(i/Rotation_Steps),1);
    else
        result_dist_matrix(l,2) = result_dist_matrix(l,2) + l_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
        result_dist_matrix(u,2) = result_dist_matrix(u,2) + u_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
    end

end

%create result matrix
result_histogram_matrix = zeros(bins,bins);

%find min/max xyz values for each carbon
c_max = max(c_xyz);
c_min = min(c_xyz);
c_diff = 1.2*c_max - 1.2*c_min; %add empty space on either side
c_steps = c_diff / bins;

%define axis limits
x_axis = (1.2*c_min(1):c_steps(1):1.2*c_max(1));
z_axis = (1.2*c_min(2):c_steps(2):1.2*c_max(2));

for i = 1:length(c_xyz)
    %determine which bin simulation should go into
    xbin = (c_xyz(i,1) - x_axis(1)) / c_steps(1);
    zbin = (c_xyz(i,2) - z_axis(1)) / c_steps(2);
    
    lx = floor(xbin + 1);
    ux = ceil(xbin + 1);
    b1x = x_axis(lx);
    b2x = x_axis(ux);
    lx_p = ((b2x)-(c_xyz(i,1)))/(b2x-b1x);
    ux_p = 1 - lx_p;
    
    lz = floor(zbin + 1);
    uz = ceil(zbin + 1);
    b1z = z_axis(lz);
    b2z = z_axis(uz);
    lz_p = ((b2z)-(c_xyz(i,2)))/(b2z-b1z);
    uz_p = 1 - lz_p;
    
    %add intensity to result matrix
    if i <= length(c_xyz)/2
        result_histogram_matrix(lz,lx) = result_histogram_matrix(lz,lx) + lz_p*lx_p*result_intensity(ceil(i/Rotation_Steps),1);
        result_histogram_matrix(uz,lx) = result_histogram_matrix(uz,lx) + uz_p*lx_p*result_intensity(ceil(i/Rotation_Steps),1);
        result_histogram_matrix(lz,ux) = result_histogram_matrix(lz,ux) + lz_p*ux_p*result_intensity(ceil(i/Rotation_Steps),1);
        result_histogram_matrix(uz,ux) = result_histogram_matrix(uz,ux) + uz_p*ux_p*result_intensity(ceil(i/Rotation_Steps),1);
    else
        result_histogram_matrix(lz,lx) = result_histogram_matrix(lz,lx) + lz_p*lx_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
        result_histogram_matrix(uz,lx) = result_histogram_matrix(uz,lx) + uz_p*lx_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
        result_histogram_matrix(lz,ux) = result_histogram_matrix(lz,ux) + lz_p*ux_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
        result_histogram_matrix(uz,ux) = result_histogram_matrix(uz,ux) + uz_p*ux_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
    end
end

%create result matrix
result_histogram_ccangles_matrix = zeros(bins,bins);

%find min/max xyz values for each carbon
ccangles1_max = -50;
ccangles1_min = 50;

ccangles2_max = -50;
ccangles2_min = 50;

ccangles1_diff = ccangles1_max - ccangles1_min; 
ccangles1_steps = ccangles1_diff / bins;

ccangles2_diff = ccangles2_max - ccangles2_min; 
ccangles2_steps = ccangles2_diff / bins;


%define axis limits
ccangles_x_axis = (ccangles1_min(1):ccangles1_steps(1):ccangles1_max(1));
ccangles_z_axis = (ccangles2_min(1):ccangles2_steps(1):ccangles2_max(1));

for i = 1:(length(c_angle_dev)/2)
    %determine which bin simulation should go into
  
    xbin = (-c_angle_dev(i,1) - ccangles_x_axis(1)) / ccangles1_steps(1);
    zbin = (c_angle_dev(i+(length(c_angle_dev)/2),1) - ccangles_z_axis(1)) / ccangles2_steps(1);
    
    lx = floor(xbin + 1);
    ux = ceil(xbin + 1);
    b1x = ccangles_x_axis(lx);
    b2x = ccangles_x_axis(ux);
    lx_p = ((b2x)-(-c_angle_dev(i,1)))/(b2x-b1x);
    ux_p = 1 - lx_p;

    lz = floor(zbin + 1);
    uz = ceil(zbin + 1);
    b1z = ccangles_z_axis(lz);
    b2z = ccangles_z_axis(uz);
    lz_p = ((b2z)-(c_angle_dev(i+(length(c_angle_dev)/2),1)))/(b2z-b1z);
    uz_p = 1 - lz_p;
    
    %add intensity to result matrix
    if i <= length(c_angle_dev)/2
        result_histogram_ccangles_matrix(lz,lx) = result_histogram_ccangles_matrix(lz,lx) + lz_p*lx_p*result_intensity(ceil(i/Rotation_Steps),1);
        result_histogram_ccangles_matrix(uz,lx) = result_histogram_ccangles_matrix(uz,lx) + uz_p*lx_p*result_intensity(ceil(i/Rotation_Steps),1);
        result_histogram_ccangles_matrix(lz,ux) = result_histogram_ccangles_matrix(lz,ux) + lz_p*ux_p*result_intensity(ceil(i/Rotation_Steps),1);
        result_histogram_ccangles_matrix(uz,ux) = result_histogram_ccangles_matrix(uz,ux) + uz_p*ux_p*result_intensity(ceil(i/Rotation_Steps),1);
    else
        result_histogram_ccangles_matrix(lz,lx) = result_histogram_ccangles_matrix(lz,lx) + lz_p*lx_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
        result_histogram_ccangles_matrix(uz,lx) = result_histogram_ccangles_matrix(uz,lx) + uz_p*lx_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
        result_histogram_ccangles_matrix(lz,ux) = result_histogram_ccangles_matrix(lz,ux) + lz_p*ux_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
        result_histogram_ccangles_matrix(uz,ux) = result_histogram_ccangles_matrix(uz,ux) + uz_p*ux_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
    end
    
end

%create result matrix
result_histogram_ccangvel_matrix = zeros(bins,bins);

%find min/max xyz values for each carbon
ccangvel1_max = -50;
ccangvel1_min = 50;

ccangvel2_max = max(c_xyz_vel);
ccangvel2_min = min(c_xyz_vel);

ccangvel1_diff = ccangvel1_max - ccangvel1_min; 
ccangvel1_steps = ccangvel1_diff / bins;

ccangvel2_diff = 1.2*ccangvel2_max - 0.8*ccangvel2_min; 
ccangvel2_steps = ccangvel2_diff / bins;


%define axis limits
ccangvel_x_axis = (ccangvel1_min(1):ccangvel1_steps(1):ccangvel1_max(1));
ccangvel_z_axis = (0.8*ccangvel2_min(1):ccangvel2_steps(1):1.2*ccangvel2_max(1));

for i = 1:(length(c_angle_dev)/2)
    %determine which bin simulation should go into
  
    xbin = (c_angle_dev(i,1) - ccangvel_x_axis(1)) / ccangvel1_steps(1);
    zbin = (c_xyz_vel(i,1) - ccangvel_z_axis(1)) / ccangvel2_steps(1);
    
    lx = floor(xbin + 1);
    ux = ceil(xbin + 1);
    b1x = ccangvel_x_axis(lx);
    b2x = ccangvel_x_axis(ux);
    lx_p = ((b2x)-(c_angle_dev(i,1)))/(b2x-b1x);
    ux_p = 1 - lx_p;

    lz = floor(zbin + 1);
    uz = ceil(zbin + 1);
    b1z = ccangvel_z_axis(lz);
    b2z = ccangvel_z_axis(uz);
    lz_p = ((b2z)-(c_xyz_vel(i,1)))/(b2z-b1z);
    uz_p = 1 - lz_p;
    
    %add intensity to result matrix
    if i <= length(c_angle_dev)/2
        result_histogram_ccangvel_matrix(lz,lx) = result_histogram_ccangvel_matrix(lz,lx) + lz_p*lx_p*result_intensity(ceil(i/Rotation_Steps),1);
        result_histogram_ccangvel_matrix(uz,lx) = result_histogram_ccangvel_matrix(uz,lx) + uz_p*lx_p*result_intensity(ceil(i/Rotation_Steps),1);
        result_histogram_ccangvel_matrix(lz,ux) = result_histogram_ccangvel_matrix(lz,ux) + lz_p*ux_p*result_intensity(ceil(i/Rotation_Steps),1);
        result_histogram_ccangvel_matrix(uz,ux) = result_histogram_ccangvel_matrix(uz,ux) + uz_p*ux_p*result_intensity(ceil(i/Rotation_Steps),1);
    else
        result_histogram_ccangvel_matrix(lz,lx) = result_histogram_ccangvel_matrix(lz,lx) + lz_p*lx_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
        result_histogram_ccangvel_matrix(uz,lx) = result_histogram_ccangvel_matrix(uz,lx) + uz_p*lx_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
        result_histogram_ccangvel_matrix(lz,ux) = result_histogram_ccangvel_matrix(lz,ux) + lz_p*ux_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
        result_histogram_ccangvel_matrix(uz,ux) = result_histogram_ccangvel_matrix(uz,ux) + uz_p*ux_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
    end
    
end


%create result matrix
result_histogram_ccvel_matrix = zeros(bins,bins);

%find min/max xyz values for each carbon
ccvel1_max = max(c_xz_vel(1:(length(c_xz_vel)/2)));
ccvel1_min = min(c_xz_vel(1:(length(c_xz_vel)/2)));

ccvel2_max = max(c_xz_vel(((length(c_xz_vel)/2)+1):length(c_xz_vel)));
ccvel2_min = min(c_xz_vel(((length(c_xz_vel)/2)+1):length(c_xz_vel)));

ccvel1_diff = 1.2*ccvel1_max - 0.8*ccvel1_min; %add empty space on either side
ccvel1_steps = ccvel1_diff / bins;

ccvel2_diff = 1.2*ccvel2_max - 0.8*ccvel2_min; %add empty space on either side
ccvel2_steps = ccvel2_diff / bins;

%define axis limits
ccvel_x_axis = (0.8*ccvel1_min(1):ccvel1_steps(1):1.2*ccvel1_max(1));
ccvel_z_axis = (0.8*ccvel2_min(1):ccvel2_steps(1):1.2*ccvel2_max(1));

for i = 1:(length(c_xz_vel)/2)
    %determine which bin simulation should go into
  
    xbin = (c_xz_vel(i,1) - ccvel_x_axis(1)) / ccvel1_steps(1);
    zbin = (c_xz_vel(i+(length(c_xz_vel)/2),1) - ccvel_z_axis(1)) / ccvel2_steps(1);

    lx = floor(xbin + 1);
    ux = ceil(xbin + 1);
    b1x = ccvel_x_axis(lx);
    b2x = ccvel_x_axis(ux);
    lx_p = ((b2x)-(c_xz_vel(i,1)))/(b2x-b1x);
    ux_p = 1 - lx_p;

    lz = floor(zbin + 1);
    uz = ceil(zbin + 1);
    b1z = ccvel_z_axis(lz);
    b2z = ccvel_z_axis(uz);
    lz_p = ((b2z)-(c_xz_vel(i+(length(c_xz_vel)/2),1)))/(b2z-b1z);
    uz_p = 1 - lz_p;

    %add intensity to result matrix
    if i <= length(c_xz_vel)/2
        result_histogram_ccvel_matrix(lz,lx) = result_histogram_ccvel_matrix(lz,lx) + lz_p*lx_p*result_intensity(ceil(i/Rotation_Steps),1);
        result_histogram_ccvel_matrix(uz,lx) = result_histogram_ccvel_matrix(uz,lx) + uz_p*lx_p*result_intensity(ceil(i/Rotation_Steps),1);
        result_histogram_ccvel_matrix(lz,ux) = result_histogram_ccvel_matrix(lz,ux) + lz_p*ux_p*result_intensity(ceil(i/Rotation_Steps),1);
        result_histogram_ccvel_matrix(uz,ux) = result_histogram_ccvel_matrix(uz,ux) + uz_p*ux_p*result_intensity(ceil(i/Rotation_Steps),1);
    else
        result_histogram_ccvel_matrix(lz,lx) = result_histogram_ccvel_matrix(lz,lx) + lz_p*lx_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
        result_histogram_ccvel_matrix(uz,lx) = result_histogram_ccvel_matrix(uz,lx) + uz_p*lx_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
        result_histogram_ccvel_matrix(lz,ux) = result_histogram_ccvel_matrix(lz,ux) + lz_p*ux_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
        result_histogram_ccvel_matrix(uz,ux) = result_histogram_ccvel_matrix(uz,ux) + uz_p*ux_p*result_intensity(ceil((i-(length(c_xyz)/2))/Rotation_Steps),1);
    end
    
end

result_filtered_simulations( ~any(result_filtered_simulations,2), : ) = [];

bar(result_dist_matrix(:,1),result_dist_matrix(:,2),'BarWidth',1);
xlabel('xy velocity');
ylabel('intensity');
title('4He - Velocity in xy direction');

figure

scatter3(result_c1_xyz(:,1),result_c1_xyz(:,2),result_c1_xyz(:,3),1)
hold all
scatter3(result_c2_xyz(:,1),result_c2_xyz(:,2),result_c2_xyz(:,3),1);
xlabel('x vel');
ylabel('y vel');
zlabel('z vel');
title('4He - Final Carbon Velocity Vectors Resolved');

figure

imagesc(x_axis,z_axis,result_histogram_matrix);
xlabel('x velocity');
ylabel('z velocity');
set(gca,'YDir','normal')

figure

imagesc(ccvel_x_axis,ccvel_z_axis,result_histogram_ccvel_matrix);
title('4He - Comparison of Carbon Velocities')
xlabel('carbon 1 xz velocity');
ylabel('carbon 2 xz velocity');
set(gca,'YDir','normal')

figure

imagesc(ccangles_x_axis,ccangles_z_axis,result_histogram_ccangles_matrix);
title('4He - Comparison of Carbon Angle Deviations')
xlabel('carbon 1 angle deviation');
ylabel('carbon 2 angle deviation');
set(gca,'YDir','normal')

figure

imagesc(ccangvel_x_axis,ccangvel_z_axis,result_histogram_ccangvel_matrix);
title('4He - Comparison of Carbon Velocity / Angle')
xlabel('carbon angle');
ylabel('carbon velocity');
set(gca,'YDir','normal')