%Seperate result_inital_positions out for each atom
h1_start_all = result_initial_positions(:,1:3);
c1_start_all = result_initial_positions(:,4:6);
c2_start_all = result_initial_positions(:,7:9);
h2_start_all = result_initial_positions(:,10:12);
he_start_all = result_initial_positions(:,13:12+3*He_Atoms); 


result_mid_he_xyz = zeros(Simulations_Amount*He_Atoms, 3);
result_mid_he_xy = zeros(Simulations_Amount*He_Atoms, 2);

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
  
  %perp_vector_1 = cross(init_direction_vector,result_final_velocity_vector(i,1:3));
%   perp_vector_1 = cross(init_direction_vector,[1,0,0]);
%   norm_perp_vector_1 = perp_vector_1 / norm(perp_vector_1);
%   
%   perp_vector_2 = cross(init_direction_vector, perp_vector_1);
%   norm_perp_vector_2 = perp_vector_2 / norm(perp_vector_2);
  
  c = dot(mid_he_direction_vector,norm_direction_vector);
  c2 = init_direction_vector * (c/norm(init_direction_vector));
 
  perp_vector =  mid_he_direction_vector - c2; 
  z_axis = cross(init_direction_vector,perp_vector);
 
   mid_he_xyz = [(dot(mid_he_direction_vector, norm_direction_vector)) ...
                 (dot(mid_he_direction_vector, perp_vector/norm(perp_vector))) ...
                 (dot(mid_he_direction_vector, z_axis/norm(z_axis)))];
  
%     mid_he_xy = [dot(mid_he_direction_vector, norm_direction_vector) ...
%                  dot(mid_he_direction_vector, perp_vector/norm(perp_vector))];
   result_mid_he_xyz(i * He_Atoms - (He_Atoms - 1),:) = mid_he_xyz;
   groups(1,i * He_Atoms - (He_Atoms - 1)) = 0;
%    result_mid_he_xy(i,:) = mid_he_xy;
   
    if He_Atoms > 1 
        for j = 1:(He_Atoms - 1)
            mid_he_direction_vector = [(midpoint(1)-result_initial_positions(i,13+3*j)) ...
                                       (midpoint(2)-result_initial_positions(i,14+3*j)) ...
                                       (midpoint(3)-result_initial_positions(i,15+3*j))];
            mid_he_xyz = [(dot(mid_he_direction_vector, norm_direction_vector)) ...
                          (dot(mid_he_direction_vector, perp_vector/norm(perp_vector))) ...
                          (dot(mid_he_direction_vector, z_axis/norm(z_axis)))];
            result_mid_he_xyz(i * He_Atoms - (He_Atoms - 1) + j,:) = mid_he_xyz;
            groups(1,i * He_Atoms - (He_Atoms - 1) + j) = 1;
        end
    end
            
end

%plot results

bins = 300;

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