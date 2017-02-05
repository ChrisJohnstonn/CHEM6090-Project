%Seperate result_inital_positions out for each atom
h1_start_all = result_initial_positions(:,1:3);
c1_start_all = result_initial_positions(:,4:6);
c2_start_all = result_initial_positions(:,7:9);
h2_start_all = result_initial_positions(:,10:12);
he_start_all = result_initial_positions(:,13:12+3*He_Atoms); 


result_mid_he_xyz = zeros(Simulations_Amount, 3);

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
  perp_vector_1 = cross(init_direction_vector,[1,0,0]);
  norm_perp_vector_1 = perp_vector_1 / norm(perp_vector_1);
  
  perp_vector_2 = cross(init_direction_vector, perp_vector_1);
  norm_perp_vector_2 = perp_vector_2 / norm(perp_vector_2);
  
  %calculate components of mid_he vector along initial direction vector and
  %2 perpendicular vector
  mid_he_xyz = [(dot(mid_he_direction_vector, norm_direction_vector)) ...
                (dot(mid_he_direction_vector, norm_perp_vector_1)) ...
                (dot(mid_he_direction_vector, norm_perp_vector_2))];
  
  result_mid_he_xyz(i,:) = mid_he_xyz;
            
end

%plot results

bins = 400;

result_initial_histogram_matrix = zeros(bins, bins);

mid_he_max = max(result_mid_he_xyz);
mid_he_min = min(result_mid_he_xyz);

mid_he_diff = 1.2*mid_he_max - 1.2*mid_he_min;
mid_he_steps = mid_he_diff / bins;

%define axis
mid_he_x_axis = (1.2*mid_he_min(1):mid_he_steps(1):1.2*mid_he_max(1));
mid_he_z_axis = (1.2*mid_he_min(3):mid_he_steps(3):1.2*mid_he_max(3));

for i = 1:length(result_mid_he_xyz)
    %determine which bin simulation should go into
    xbin = ceil( (result_mid_he_xyz(i,1) - mid_he_x_axis(1)) / mid_he_steps(1));
    zbin = ceil( (result_mid_he_xyz(i,3) - mid_he_z_axis(1)) / mid_he_steps(3));
    
    %add intensity to result matrix
    result_initial_histogram_matrix(zbin,xbin) = result_initial_histogram_matrix(zbin,xbin) + result_intensity(i);

end

imagesc(mid_he_x_axis,mid_he_z_axis,result_initial_histogram_matrix);
xlabel('x position');
ylabel('z position');