%Seperate result_inital_positions out for each atom
h1_start_all = result_initial_positions(:,1:3);
c1_start_all = result_initial_positions(:,4:6);
c2_start_all = result_initial_positions(:,7:9);
h2_start_all = result_initial_positions(:,10:12);
he_start_all = result_initial_positions(:,13:12+3*He_Atoms); 


result_mid_he_xyz = zeros(Simulations_Amount*He_Atoms, 3);
result_mid_he_xy = zeros(Simulations_Amount*He_Atoms, 2);
result_mid_he_z_angle = zeros(Simulations_Amount*(He_Atoms-1),1);
result_mid_he_plane_angle = zeros(Simulations_Amount*(He_Atoms-1),1);
result_test = zeros(Simulations_Amount*(He_Atoms-1),4);


groups = zeros(1,Simulations_Amount*He_Atoms);

%calculate components of midpoint-helium vector in relevant directions

for i = 1:Simulations_Amount
  midpoint = [(result_initial_positions(i,4)+result_initial_positions(i,7))/2 ...
              (result_initial_positions(i,5)+result_initial_positions(i,8))/2 ...
              (result_initial_positions(i,6)+result_initial_positions(i,9))/2];
  
  init_direction_vector = [(result_initial_positions(i,4)-result_initial_positions(i,7)) ...
                           (result_initial_positions(i,5)-result_initial_positions(i,8)) ...
                           (result_initial_positions(i,6)-result_initial_positions(i,9))];
                  
  norm_direction_vector = init_direction_vector/norm(init_direction_vector);
  
  mid_he_direction_vector = [(midpoint(1)-result_initial_positions(i,13)) ...
                             (midpoint(2)-result_initial_positions(i,14)) ...
                             (midpoint(3)-result_initial_positions(i,15))];
  
  c = dot(mid_he_direction_vector,norm_direction_vector);
  c2 = init_direction_vector * (c/norm(init_direction_vector));
 
  perp_vector =  mid_he_direction_vector - c2; 
  z_axis = cross(init_direction_vector,perp_vector);
 
   mid_he_xyz = [(dot(mid_he_direction_vector, norm_direction_vector)) ...
                 (dot(mid_he_direction_vector, perp_vector/norm(perp_vector))) ...
                 (dot(mid_he_direction_vector, z_axis/norm(z_axis)))];
  
   result_mid_he_xyz(i * He_Atoms - (He_Atoms - 1),:) = mid_he_xyz;
   groups(1,i * He_Atoms - (He_Atoms - 1)) = 0;

   if He_Atoms > 1 
        for j = 1:(He_Atoms - 1)
            mid_he_direction_vector = [(midpoint(1)-result_initial_positions(i,13+3*j)) ...
                                       (midpoint(2)-result_initial_positions(i,14+3*j)) ...
                                       (midpoint(3)-result_initial_positions(i,15+3*j))];
            mid_he_xyz = [(dot(mid_he_direction_vector, norm_direction_vector)) ...
                          (dot(mid_he_direction_vector, perp_vector/norm(perp_vector))) ...
                          (dot(mid_he_direction_vector, z_axis/norm(z_axis)))];
            result_mid_he_xyz(i * He_Atoms - (He_Atoms - 1) + j,:) = mid_he_xyz;
            groups(1,i * He_Atoms - (He_Atoms - 1) + j) = 1 + (j-1);
            
            result_mid_he_z_angle(i*(He_Atoms - 1) - ((He_Atoms - 1) - j),1) = 90 - acosd( dot(mid_he_direction_vector,z_axis) / (norm(mid_he_direction_vector) * norm(z_axis)));
            
            mid_he_plane_proj = mid_he_direction_vector - ((dot(mid_he_direction_vector,z_axis/norm(z_axis))) * z_axis/norm(z_axis));
            result_mid_he_plane_angle(i*(He_Atoms - 1) - ((He_Atoms - 1) - j),1) = acosd( dot(perp_vector, mid_he_plane_proj) / (norm(perp_vector) * norm(mid_he_plane_proj)));
            
            if result_intensity(i) > 2.7
                
                result_test(i*(He_Atoms - 1) - ((He_Atoms - 1) - j),1:3) = mid_he_xyz;
                result_test(i*(He_Atoms - 1) - ((He_Atoms - 1) - j),4) = result_intensity(i);
            
            end
        end
        
        
    end
            
end

%plot results

bins = 500;

result_initial_histogram_matrix = zeros(bins, bins);

mid_he_max = max(result_mid_he_xyz);
mid_he_min = min(result_mid_he_xyz);

mid_he_diff = 1.2 * mid_he_max - [1.2,1.2,1.2].*[mid_he_min];
mid_he_steps = mid_he_diff / bins;

%define axis
mid_he_x_axis = (1.2*mid_he_min(1):mid_he_steps(1):1.2*mid_he_max(1));
mid_he_y_axis = (1.2*mid_he_min(2):mid_he_steps(2):1.2*mid_he_max(2));
mid_he_z_axis = (1.2*mid_he_min(3):mid_he_steps(3):1.2*mid_he_max(3));


for i = 1:Simulations_Amount
    %determine which bin simulation should go into
    
    %in plane
%     xbin = ceil( (result_mid_he_xyz(i * He_Atoms - (He_Atoms - 1),1) - mid_he_x_axis(1)) / mid_he_steps(1));
%     ybin = ceil( (result_mid_he_xyz(i * He_Atoms - (He_Atoms - 1),2) - mid_he_y_axis(1)) / mid_he_steps(2));
    
    %out of plane
    if He_Atoms > 1
        for j = 1:(He_Atoms - 1)
            xbin = ceil( (result_mid_he_xyz(i * He_Atoms - (He_Atoms - (j+1)),1) - mid_he_x_axis(1)) / mid_he_steps(1));
            ybin = ceil( (result_mid_he_xyz(i * He_Atoms - (He_Atoms - (j+1)),2) - mid_he_y_axis(1)) / mid_he_steps(2));
            
            result_initial_histogram_matrix(ybin,xbin) = result_initial_histogram_matrix(ybin,xbin) + result_intensity(i);
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
    
    result_angle_matrix_mid_he_z(l, 2) = result_angle_matrix_mid_he_z(l, 2) + l_p*result_intensity(ceil(i/(He_Atoms-1)));
    result_angle_matrix_mid_he_z(u, 2) = result_angle_matrix_mid_he_z(u, 2) + u_p*result_intensity(ceil(i/(He_Atoms-1)));
    
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
    result_angle_matrix_mid_he_plane(l, 2) = result_angle_matrix_mid_he_plane(l, 2) + l_p*result_intensity(ceil(i/(He_Atoms-1)));
    result_angle_matrix_mid_he_plane(u, 2) = result_angle_matrix_mid_he_plane(u, 2) + u_p*result_intensity(ceil(i/(He_Atoms-1)));
end

for i = 1:length(result_test)
   if result_test(Simulations_Amount*(He_Atoms-1) - (i - 1),4) == 0
       result_test(Simulations_Amount*(He_Atoms-1) - (i - 1),:) = [];
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
gscatter3(result_mid_he_xyz(:,1),result_mid_he_xyz(:,2),result_mid_he_xyz(:,3),groups);
xlabel('x position');
ylabel('y position');
zlabel('z position');
figure

bar(result_angle_matrix_mid_he_z(:,1),result_angle_matrix_mid_he_z(:,2),'BarWidth',1)
xlabel('Angle between Plane and Helium');
ylabel('Intensity');
xlim([-90 90]);
title('4He - Angle between Plane and Helium (N/S)')
figure

bar(result_angle_matrix_mid_he_plane(:,1),result_angle_matrix_mid_he_plane(:,2),'BarWidth',1)
xlabel('Angle on Plane');
ylabel('Intensity');
title('4He - Angle on Plane')
xlim([0 180]);