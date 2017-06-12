%Seperate result_inital_positions out for each atom
h1_start_all = result_initial_positions(:,1:3);
c1_start_all = result_initial_positions(:,4:6);
c2_start_all = result_initial_positions(:,7:9);
h2_start_all = result_initial_positions(:,10:12);
he_start_all = result_initial_positions(:,13:12+3*He_Atoms); 

Filtered_Simulations = length(result_filtered_simulations);

result_mid_he_xyz = zeros(Filtered_Simulations*He_Atoms, 3);
result_mid_he_xy = zeros(Filtered_Simulations*He_Atoms, 2);
result_mid_he_z_angle = zeros(Filtered_Simulations*(He_Atoms-1),1);
result_mid_he_plane_angle = zeros(Filtered_Simulations*(He_Atoms-1),1);

result_he_he_plane_angle = zeros(Filtered_Simulations*(0.5*He_Atoms*(He_Atoms-1)),1);



groups = zeros(1,Filtered_Simulations*He_Atoms);

%calculate components of midpoint-helium vector in relevant directions

for i = 1:Filtered_Simulations
    j = result_filtered_simulations(i,1);
  midpoint = [(result_initial_positions(j,4)+result_initial_positions(j,7))/2 ...
              (result_initial_positions(j,5)+result_initial_positions(j,8))/2 ...
              (result_initial_positions(j,6)+result_initial_positions(j,9))/2];
  
  init_direction_vector = [(result_initial_positions(j,4)-result_initial_positions(j,7)) ...
                           (result_initial_positions(j,5)-result_initial_positions(j,8)) ...
                           (result_initial_positions(j,6)-result_initial_positions(j,9))];
                  
  norm_direction_vector = init_direction_vector/norm(init_direction_vector);
  
  mid_he_direction_vector = [(midpoint(1)-result_initial_positions(j,13)) ...
                             (midpoint(2)-result_initial_positions(j,14)) ...
                             (midpoint(3)-result_initial_positions(j,15))];
  
  c = dot(mid_he_direction_vector,norm_direction_vector);
  c2 = init_direction_vector * (c/norm(init_direction_vector));
 
  perp_vector =  mid_he_direction_vector - c2; 
  z_axis = cross(init_direction_vector,perp_vector);
   
   %v_rot = v*cos(theta) + (k x v)*sin(theta) + k(k.v)(1-cos(theta)
   theta = result_filtered_simulations(i,2) * (Rotation_Amount)/(Rotation_Steps);
   mid_he_rot = mid_he_direction_vector*cos(theta) + cross(norm_direction_vector,mid_he_direction_vector)*sin(theta) + ...
                norm_direction_vector*dot(norm_direction_vector,mid_he_direction_vector)*(1-cos(theta));
   
   mid_he_xyz = [(dot(mid_he_rot, norm_direction_vector)) ...
                 (dot(mid_he_rot, perp_vector/norm(perp_vector))) ...
                 (dot(mid_he_rot, z_axis/norm(z_axis)))];
  
             
   result_mid_he_xyz(i * He_Atoms - (He_Atoms - 1),:) = mid_he_xyz;
   groups(1,i * He_Atoms - (He_Atoms - 1)) = 0;

   if He_Atoms > 1 
        for k = 1:(He_Atoms - 1)
            new_mid_he_direction_vector = [(midpoint(1)-result_initial_positions(j,13+3*k)) ...
                                       (midpoint(2)-result_initial_positions(j,14+3*k)) ...
                                       (midpoint(3)-result_initial_positions(j,15+3*k))];
            mid_he_rot = new_mid_he_direction_vector*cos(theta) + cross(norm_direction_vector,new_mid_he_direction_vector)*sin(theta) + ...
                         norm_direction_vector*dot(norm_direction_vector,new_mid_he_direction_vector)*(1-cos(theta));
            mid_he_xyz = [(dot(mid_he_rot, norm_direction_vector)) ...
                          (dot(mid_he_rot, perp_vector/norm(perp_vector))) ...
                          (dot(mid_he_rot, z_axis/norm(z_axis)))];
            result_mid_he_xyz(i * He_Atoms - (He_Atoms - 1) + k,:) = mid_he_xyz;
            groups(1,i * He_Atoms - (He_Atoms - 1) + k) = 1 + (k-1);
            
            result_mid_he_z_angle(i*(He_Atoms - 1) - ((He_Atoms - 1) - k),1) = 90 - acosd( dot(new_mid_he_direction_vector,z_axis) / (norm(new_mid_he_direction_vector) * norm(z_axis)));
            
            mid_he_plane_proj = new_mid_he_direction_vector - ((dot(new_mid_he_direction_vector,z_axis/norm(z_axis))) * z_axis/norm(z_axis));
            result_mid_he_plane_angle(i*(He_Atoms - 1) - ((He_Atoms - 1) - k),1) = acosd( dot(perp_vector, mid_he_plane_proj) / (norm(perp_vector) * norm(mid_he_plane_proj)));
            
            result_he_he_plane_angle(i*(He_Atoms - 1) - ((He_Atoms - 1) - k),1) = acosd( dot(mid_he_direction_vector, mid_he_plane_proj) / (norm(mid_he_direction_vector) * norm(mid_he_plane_proj)));
            
        end
        
        
    end
            
end

%plot results

bins = 600;

result_initial_histogram_matrix = zeros(bins, bins);

mid_he_max = max(result_mid_he_xyz);
mid_he_min = min(result_mid_he_xyz);

mid_he_diff = 1.2 * mid_he_max - [1.2,1.2,1.2].*[mid_he_min];
mid_he_steps = mid_he_diff / bins;

%define axis
mid_he_x_axis = (1.2*mid_he_min(1):mid_he_steps(1):1.2*mid_he_max(1));
mid_he_y_axis = (1.2*mid_he_min(2):mid_he_steps(2):1.2*mid_he_max(2));
mid_he_z_axis = (1.2*mid_he_min(3):mid_he_steps(3):1.2*mid_he_max(3));


for i = 1:Filtered_Simulations
    %determine which bin simulation should go into
    
    %in plane
%     xbin = ceil( (result_mid_he_xyz(i * He_Atoms - (He_Atoms - 1),1) - mid_he_x_axis(1)) / mid_he_steps(1));
%     ybin = ceil( (result_mid_he_xyz(i * He_Atoms - (He_Atoms - 1),2) - mid_he_y_axis(1)) / mid_he_steps(2));
    
    %out of plane
    if He_Atoms > 1
        for j = 1:(He_Atoms - 1)
            xbin = ceil( (result_mid_he_xyz(i * He_Atoms - (He_Atoms - (j+1)),1) - mid_he_x_axis(1)) / mid_he_steps(1));
            ybin = ceil( (result_mid_he_xyz(i * He_Atoms - (He_Atoms - (j+1)),2) - mid_he_y_axis(1)) / mid_he_steps(2));
            
            result_initial_histogram_matrix(ybin,xbin) = result_initial_histogram_matrix(ybin,xbin) + result_intensity(result_filtered_simulations(i,1));
        end
    end
    
    %add intensity to result matrix
    %result_initial_histogram_matrix(ybin,xbin) = result_initial_histogram_matrix(ybin,xbin) + result_intensity(i);

end


mid_he_z_angle_max = 90;
mid_he_z_angle_min = -90;
mid_he_z_angle_diff = mid_he_z_angle_max - mid_he_z_angle_min;
mid_he_z_angle_steps = mid_he_z_angle_diff / bins;

mid_he_z_angle_axis = (mid_he_z_angle_min:mid_he_z_angle_steps:mid_he_z_angle_max);

result_angle_matrix_mid_he_z = zeros(length(mid_he_z_angle_axis),2);
result_angle_matrix_mid_he_z(:,1) = mid_he_z_angle_axis;

for i = 1:length(result_mid_he_z_angle)
    
    mid_he_z_angle_x_bin = (result_mid_he_z_angle(i) - mid_he_z_angle_axis(1) ) / mid_he_z_angle_steps ;
    l = floor(mid_he_z_angle_x_bin + 1);
    u = ceil(mid_he_z_angle_x_bin + 1);
    b1 = mid_he_z_angle_axis(l);
    b2 = mid_he_z_angle_axis(u);
    l_p = ((b2)-(result_mid_he_z_angle(i)))/(b2-b1);
    u_p = 1 - l_p;
    
    result_angle_matrix_mid_he_z(l, 2) = result_angle_matrix_mid_he_z(l, 2) + l_p*result_intensity(ceil(result_filtered_simulations(ceil(i/(He_Atoms-1)),1)/(He_Atoms-1)));
    result_angle_matrix_mid_he_z(u, 2) = result_angle_matrix_mid_he_z(u, 2) + u_p*result_intensity(ceil(result_filtered_simulations(ceil(i/(He_Atoms-1)),1)/(He_Atoms-1)));
    
end

mid_he_plane_angle_max = 180;
mid_he_plane_angle_min = 0;
mid_he_plane_angle_diff = mid_he_plane_angle_max - mid_he_plane_angle_min;
mid_he_plane_angle_steps = mid_he_plane_angle_diff / bins;

mid_he_plane_angle_axis = (mid_he_plane_angle_min:mid_he_plane_angle_steps:mid_he_plane_angle_max);

result_angle_matrix_mid_he_plane = zeros(length(mid_he_plane_angle_axis),2);
result_angle_matrix_mid_he_plane(:,1) = mid_he_plane_angle_axis;

for i = 1:length(result_mid_he_plane_angle)
    
    %mid_he_plane_angle_x_bin = ceil( (result_mid_he_plane_angle(i) - mid_he_plane_angle_axis(1) ) / mid_he_plane_angle_steps );
    mid_he_plane_angle_x_bin = ( result_mid_he_plane_angle(i) - mid_he_plane_angle_axis(1) ) / mid_he_plane_angle_steps;
    l = floor(mid_he_plane_angle_x_bin + 1);
    u = ceil(mid_he_plane_angle_x_bin + 1);
    b1 = mid_he_plane_angle_axis(l);
    b2 = mid_he_plane_angle_axis(u);
    l_p = ((b2)-(result_mid_he_plane_angle(i)))/(b2-b1);
    u_p = 1 - l_p;
    
    %result_angle_matrix_mid_he_plane(mid_he_plane_angle_x_bin, 2) = result_angle_matrix_mid_he_plane(mid_he_plane_angle_x_bin, 2) + result_intensity(ceil(i/(He_Atoms-1)));
    result_angle_matrix_mid_he_plane(l, 2) = result_angle_matrix_mid_he_plane(l, 2) + l_p*result_intensity(ceil(result_filtered_simulations(ceil(i/(He_Atoms-1)),1)/(He_Atoms-1)));
    result_angle_matrix_mid_he_plane(u, 2) = result_angle_matrix_mid_he_plane(u, 2) + u_p*result_intensity(ceil(result_filtered_simulations(ceil(i/(He_Atoms-1)),1)/(He_Atoms-1)));
end

he_he_plane_angle_max = 180;
he_he_plane_angle_min = 0;
he_he_plane_angle_diff = he_he_plane_angle_max - he_he_plane_angle_min;
he_he_plane_angle_steps = he_he_plane_angle_diff / bins;

he_he_plane_angle_axis = (he_he_plane_angle_min:he_he_plane_angle_steps:he_he_plane_angle_max);

result_angle_matrix_he_he_plane = zeros(length(he_he_plane_angle_axis),2);
result_angle_matrix_he_he_plane(:,1) = he_he_plane_angle_axis;

for i = 1:length(result_he_he_plane_angle)
    
    
    he_he_plane_angle_x_bin = ( result_he_he_plane_angle(i) - he_he_plane_angle_axis(1) ) / he_he_plane_angle_steps;
    l = floor(he_he_plane_angle_x_bin + 1);
    u = ceil(he_he_plane_angle_x_bin + 1);
    b1 = he_he_plane_angle_axis(l);
    b2 = he_he_plane_angle_axis(u);
    l_p = ((b2)-(result_he_he_plane_angle(i)))/(b2-b1);
    u_p = 1 - l_p;
    
    %result_angle_matrix_mid_he_plane(mid_he_plane_angle_x_bin, 2) = result_angle_matrix_mid_he_plane(mid_he_plane_angle_x_bin, 2) + result_intensity(ceil(i/(He_Atoms-1)));
    result_angle_matrix_he_he_plane(l, 2) = result_angle_matrix_he_he_plane(l, 2) + l_p*result_intensity(ceil(result_filtered_simulations(ceil(i/(He_Atoms-1)),1)/(He_Atoms-1)));
    result_angle_matrix_he_he_plane(u, 2) = result_angle_matrix_he_he_plane(u, 2) + u_p*result_intensity(ceil(result_filtered_simulations(ceil(i/(He_Atoms-1)),1)/(He_Atoms-1)));
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

for i = 1:Filtered_Simulations
    %determine which bin simulation should go into
    Simulation_Number = result_filtered_simulations(i,1)-1;
    Rotation = result_filtered_simulations(i,2);
    xbin = (c_xz_vel(Simulation_Number*Rotation_Steps + Rotation,1) - ccvel_x_axis(1)) / ccvel1_steps(1);
    zbin = (c_xz_vel((Simulation_Number*Rotation_Steps + Rotation)+0.5*length(c_xz_vel),1) - ccvel_z_axis(1)) / ccvel2_steps(1);

    lx = floor(xbin + 1);
    ux = ceil(xbin + 1);
    b1x = ccvel_x_axis(lx);
    b2x = ccvel_x_axis(ux);
    lx_p = ((b2x)-(c_xz_vel(Simulation_Number*Rotation_Steps + Rotation,1)))/(b2x-b1x);
    ux_p = 1 - lx_p;

    lz = floor(zbin + 1);
    uz = ceil(zbin + 1);
    b1z = ccvel_z_axis(lz);
    b2z = ccvel_z_axis(uz);
    lz_p = ((b2z)-(c_xz_vel((Simulation_Number*Rotation_Steps + Rotation)+0.5*length(c_xz_vel),1)))/(b2z-b1z);
    uz_p = 1 - lz_p;

    %add intensity to result matrix
    if i <= length(c_xz_vel)/2
        result_histogram_ccvel_matrix(lz,lx) = result_histogram_ccvel_matrix(lz,lx) + lz_p*lx_p*result_intensity(Simulation_Number+1,1);
        result_histogram_ccvel_matrix(uz,lx) = result_histogram_ccvel_matrix(uz,lx) + uz_p*lx_p*result_intensity(Simulation_Number+1,1);
        result_histogram_ccvel_matrix(lz,ux) = result_histogram_ccvel_matrix(lz,ux) + lz_p*ux_p*result_intensity(Simulation_Number+1,1);
        result_histogram_ccvel_matrix(uz,ux) = result_histogram_ccvel_matrix(uz,ux) + uz_p*ux_p*result_intensity(Simulation_Number+1,1);
    else
        result_histogram_ccvel_matrix(lz,lx) = result_histogram_ccvel_matrix(lz,lx) + lz_p*lx_p*result_intensity(Simulation_Number+1,1);
        result_histogram_ccvel_matrix(uz,lx) = result_histogram_ccvel_matrix(uz,lx) + uz_p*lx_p*result_intensity(Simulation_Number+1,1);
        result_histogram_ccvel_matrix(lz,ux) = result_histogram_ccvel_matrix(lz,ux) + lz_p*ux_p*result_intensity(Simulation_Number+1,1);
        result_histogram_ccvel_matrix(uz,ux) = result_histogram_ccvel_matrix(uz,ux) + uz_p*ux_p*result_intensity(Simulation_Number+1,1);
    end
    
end

imagesc(mid_he_x_axis,mid_he_y_axis,result_initial_histogram_matrix);
xlabel('x position');
ylabel('y position');
figure

% scatter3(result_mid_he_xyz(:,1),result_mid_he_xyz(:,2),result_mid_he_xyz(:,3),1)
% xlabel('x position');
% ylabel('y position');
% zlabel('z position');
% figure
% gscatter3(result_mid_he_xyz(:,1),result_mid_he_xyz(:,2),result_mid_he_xyz(:,3),groups);
% xlabel('x position');
% ylabel('y position');
% zlabel('z position');
% title('2He - He Starting positions (xz / 2.14 - 2.19e4)');
% figure

bar(result_angle_matrix_mid_he_z(:,1),result_angle_matrix_mid_he_z(:,2),'BarWidth',1)
xlabel('Angle between Plane and Helium');
ylabel('Intensity');
xlim([-90 90]);
title('2He - Angle between Plane and Helium (N/S) (xz / 2.14 - 2.19e4)')
figure

bar(result_angle_matrix_mid_he_plane(:,1),result_angle_matrix_mid_he_plane(:,2),'BarWidth',1)
xlabel('Angle on Plane');
ylabel('Intensity');
title('2He - Angle on Plane (xz / 2.14 - 2.19e4)')
xlim([0 180]);
figure

bar(result_angle_matrix_he_he_plane(:,1),result_angle_matrix_he_he_plane(:,2),'BarWidth',1)
xlabel('Angle on Plane between Heliums');
ylabel('Intensity');
title('2He - Angle on Plane between Heliums (xz / 2.14 - 2.19e4)')
xlim([0 180]);

figure

imagesc(ccvel_x_axis,ccvel_z_axis,result_histogram_ccvel_matrix);
title('2He - Comparison of Carbon Velocities')
xlabel('carbon 1 xz velocity');
ylabel('carbon 2 xz velocity');
set(gca,'YDir','normal')