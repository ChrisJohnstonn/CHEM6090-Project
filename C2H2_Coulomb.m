
%Ask user for data file and read it
rawdata = uigetfile('.psi_asc','Please select a data file for C2H2');
data = dlmread(rawdata);

%Retrieve several properties from data file.
Data_Length = size(data,1);
Total_Atoms = data(1,1);
He_Atoms = Total_Atoms - 2;
Simulations_Amount = (Data_Length - 1)/(Total_Atoms + 1); %Each simulation has Total_Atoms + 1 lines.
%Simulations_Amount = 1;

%Time Step
dt = 1e-15; %In seconds
Total_Steps = 1000;

%Constants
%Atomic Mass Unit
amu = 1.6605e-27;

%Meters in 1 Bohr radius
bohr_to_meter = 5.2918e-11;

%Masses

m_c = mass('C') * amu;
m_h = mass('H') * amu;
m_he = mass('He') * amu;

%C2H2 bond lengths
r_CC = 1.20830e-10;
r_CH = 1.05756e-10;

%Coulomb's Constant
k = 8.9876e09;
%Elementary Charge
q = 1.602e-19;
%Charges on atoms
q_c = 2*q;
q_h = 1*q;
q_he = 1*q;

%Create Results Matrices
result_intensity = zeros(Simulations_Amount,1);
result_initial_vector = zeros(Simulations_Amount,3);
result_final_velocity_vector = zeros(Simulations_Amount,6);

Total_Kinetic = 0;


for i = 1:Simulations_Amount
    
    if mod(i,100) == 0
        disp(fprintf(string('Starting simulation ') + i + string(' of ') + Simulations_Amount));
    end
   
    %Setup up simulation
    
    %Populate intensity results
    result_intensity(i) = data((i*(Total_Atoms+1)) - (He_Atoms + 1),1);
    
    %Calculate starting x,y,z positions.
    %Results in data given by 2 pseudo-atoms.
    %Calculate actual positions by finding midpoint and
    %using bond lengths to determine coordinates.
    %
    %Pull pseudoatom coordinates from data and convert to meters from atomic units.
    p1_coordinates = [data((i*(Total_Atoms+1)) - (He_Atoms),1)*bohr_to_meter ... 
                      data((i*(Total_Atoms+1)) - (He_Atoms),2)*bohr_to_meter ...
                      data((i*(Total_Atoms+1)) - (He_Atoms),3)*bohr_to_meter];
    p2_coordinates = [data((i*(Total_Atoms+1)) - (He_Atoms - 1),1)*bohr_to_meter ... 
                      data((i*(Total_Atoms+1)) - (He_Atoms - 1),2)*bohr_to_meter ...
                      data((i*(Total_Atoms+1)) - (He_Atoms - 1),3)*bohr_to_meter];
    %Find midpoint of 2 pseudo atoms.
    midpoint_coordinates = [((p1_coordinates(1,1) + p2_coordinates(1,1))/2) ...
                            ((p1_coordinates(1,2) + p2_coordinates(1,2))/2) ...
                            ((p1_coordinates(1,3) + p2_coordinates(1,3))/2)];
    
    %Create and normalise vector between 2 pseudoatoms.
    v = [((p1_coordinates(1,1) - p2_coordinates(1,1))) ...
         ((p1_coordinates(1,2) - p2_coordinates(1,2))) ...
         ((p1_coordinates(1,3) - p2_coordinates(1,3)))]; 
    
    u = v/norm(v);
    
    %Using normalised vector, midpoint and bond lengths, determine C2H2
    %atom starting positions.
    %A point distance d from point (x,y,z) has coordinates (x,y,z)+ d*u where u
    %is the normalised vector along the line.
    
    h1_start = midpoint_coordinates + (0.5*r_CC + r_CH)*(u);
    c1_start = midpoint_coordinates + (0.5*r_CC)*(u);
    c2_start = midpoint_coordinates - (0.5*r_CC)*(u);
    h2_start = midpoint_coordinates - (0.5*r_CC + r_CH)*(u);
    
    he_start = zeros(He_Atoms,3);
    for j = 1:He_Atoms %Create matrix with each row being xyz coordinates of a He atom in angstroms.
        he_start(j,1) = data((i*(Total_Atoms+1)) - (He_Atoms - (1+j)),1)*bohr_to_meter;
        he_start(j,2) = data((i*(Total_Atoms+1)) - (He_Atoms - (1+j)),2)*bohr_to_meter;
        he_start(j,3) = data((i*(Total_Atoms+1)) - (He_Atoms - (1+j)),3)*bohr_to_meter;
    end   
    
    
    %Create array to keep track of current properties through simulation
    %mass|charge|r_x|r_y|r_z|v|v_x|v_y|v_z|f|f_x|f_y|f_z|f_x_o|f_y_o|f_z_o|ke
    h1_properties = [m_h q_h h1_start(1,1) h1_start(1,2) h1_start(1,3) 0 0 0 0 0 0 0 0 0 0 0 0];
    h2_properties = [m_h q_h h2_start(1,1) h2_start(1,2) h2_start(1,3) 0 0 0 0 0 0 0 0 0 0 0 0];
    c1_properties = [m_c q_c c1_start(1,1) c1_start(1,2) c1_start(1,3) 0 0 0 0 0 0 0 0 0 0 0 0];
    c2_properties = [m_c q_c c2_start(1,1) c2_start(1,2) c2_start(1,3) 0 0 0 0 0 0 0 0 0 0 0 0];
    
    he_properties = zeros(He_Atoms, 17);
    for c = 1:He_Atoms
        he_properties(c,:) = [m_he q_he he_start(c,1) he_start(c,2) he_start(c,3) 0 0 0 0 0 0 0 0 0 0 0 0];
    end
    
    
    for m = 1:Total_Steps
       %Velocity Verlet Algorithm
       %1. Given current position and velocity, compute forces.
       %2. Update position
       %3. Repeat
       
       %Distance Between Atoms
       %Pythagoras 
       
       %Ignoring Helium initially.
       h1_distance = zeros(1,Total_Atoms);
       h1_distance(1,1) = 0; 
       h1_distance(1,2) = (sqrt((h1_properties(3) - c1_properties(3))^2 + (h1_properties(4) - c1_properties(4))^2 + (h1_properties(5) - c1_properties(5))^2)); 
       h1_distance(1,3) = (sqrt((h1_properties(3) - c2_properties(3))^2 + (h1_properties(4) - c2_properties(4))^2 + (h1_properties(5) - c2_properties(5))^2)); 
       h1_distance(1,4) = (sqrt((h1_properties(3) - h2_properties(3))^2 + (h1_properties(4) - h2_properties(4))^2 + (h1_properties(5) - h2_properties(5))^2));
       
       h2_distance = zeros(1,Total_Atoms);
       h2_distance(1,1) = (sqrt((h2_properties(3) - h1_properties(3))^2 + (h2_properties(4) - h1_properties(4))^2 + (h2_properties(5) - h1_properties(5))^2));
       h2_distance(1,2) = (sqrt((h2_properties(3) - c1_properties(3))^2 + (h2_properties(4) - c1_properties(4))^2 + (h2_properties(5) - c1_properties(5))^2));
       h2_distance(1,3) = (sqrt((h2_properties(3) - c2_properties(3))^2 + (h2_properties(4) - c2_properties(4))^2 + (h2_properties(5) - c2_properties(5))^2));
       h2_distance(1,4) = 0;
       
       c1_distance = zeros(1,Total_Atoms);
       c1_distance(1,1) = (sqrt((c1_properties(3) - h1_properties(3))^2 + (c1_properties(4) - h1_properties(4))^2 + (c1_properties(5) - h1_properties(5))^2));
       c1_distance(1,2) = 0;
       c1_distance(1,3) = (sqrt((c1_properties(3) - c2_properties(3))^2 + (c1_properties(4) - c2_properties(4))^2 + (c1_properties(5) - c2_properties(5))^2));
       c1_distance(1,4) = (sqrt((c1_properties(3) - h2_properties(3))^2 + (c1_properties(4) - h2_properties(4))^2 + (c1_properties(5) - h2_properties(5))^2));
       
       c2_distance = zeros(1,Total_Atoms);
       c2_distance(1,1) = (sqrt((c2_properties(3) - h1_properties(3))^2 + (c2_properties(4) - h1_properties(4))^2 + (c2_properties(5) - h1_properties(5))^2));
       c2_distance(1,2) = (sqrt((c2_properties(3) - c1_properties(3))^2 + (c2_properties(4) - c1_properties(4))^2 + (c2_properties(5) - c1_properties(5))^2));
       c2_distance(1,3) = 0; 
       c2_distance(1,4) = (sqrt((c2_properties(3) - h2_properties(3))^2 + (c2_properties(4) - h2_properties(4))^2 + (c2_properties(5) - h2_properties(5))^2));
       
       he_distance = zeros(He_Atoms,Total_Atoms); 
       
       for n = 1:He_Atoms
           h1_distance(1,n+4) = (sqrt((h1_properties(3) - he_properties(n,3))^2 + (h1_properties(4) - he_properties(n,4))^2 + (h1_properties(5) - he_properties(n,5))^2));
           h2_distance(1,n+4) = (sqrt((h2_properties(3) - he_properties(n,3))^2 + (h2_properties(4) - he_properties(n,4))^2 + (h2_properties(5) - he_properties(n,5))^2));
           c1_distance(1,n+4) = (sqrt((c1_properties(3) - he_properties(n,3))^2 + (c1_properties(4) - he_properties(n,4))^2 + (c1_properties(5) - he_properties(n,5))^2));
           c2_distance(1,n+4) = (sqrt((c2_properties(3) - he_properties(n,3))^2 + (c2_properties(4) - he_properties(n,4))^2 + (c2_properties(5) - he_properties(n,5))^2)); 
           
           he_distance(n,1) = (sqrt((he_properties(n,3) - h1_properties(3))^2 + (he_properties(n,4) - h1_properties(4))^2 + (he_properties(n,5) - h1_properties(5))^2));
           he_distance(n,2) = (sqrt((he_properties(n,3) - c1_properties(3))^2 + (he_properties(n,4) - c1_properties(4))^2 + (he_properties(n,5) - c1_properties(5))^2)); 
           he_distance(n,3) = (sqrt((he_properties(n,3) - c2_properties(3))^2 + (he_properties(n,4) - c2_properties(4))^2 + (he_properties(n,5) - c2_properties(5))^2)); 
           he_distance(n,4) = (sqrt((he_properties(n,3) - h2_properties(3))^2 + (he_properties(n,4) - h2_properties(4))^2 + (he_properties(n,5) - h2_properties(5))^2));
           
           for o = 1:He_Atoms    
                he_distance(o,n+4) = (sqrt((he_properties(o,3) - he_properties(n,3))^2 + (he_properties(o,4) - he_properties(n,4))^2 + (he_properties(o,5) - he_properties(n,5))^2));
           end
           
       end
        
       %Force Between Atoms
       %Coulombic Forces -> F = (k*q1*q2)/(r^2)
       h1_forces = zeros(1,Total_Atoms);
       h1_forces(1,1) = 0;
       h1_forces(1,2) = (k * h1_properties(2) * c1_properties(2))/(h1_distance(1,2)^2);
       h1_forces(1,3) = (k * h1_properties(2) * c2_properties(2))/(h1_distance(1,3)^2);
       h1_forces(1,4) = (k * h1_properties(2) * h2_properties(2))/(h1_distance(1,4)^2);
       
       h2_forces = zeros(1,Total_Atoms);
       h2_forces(1,1) = (k * h2_properties(2) * h1_properties(2))/(h2_distance(1,1)^2);
       h2_forces(1,2) = (k * h2_properties(2) * c1_properties(2))/(h2_distance(1,2)^2);
       h2_forces(1,3) = (k * h2_properties(2) * c2_properties(2))/(h2_distance(1,3)^2);
       h2_forces(1,4) = 0;
       
       c1_forces = zeros(1,Total_Atoms);
       c1_forces(1,1) = (k * c1_properties(2) * h1_properties(2))/(c1_distance(1,1)^2);
       c1_forces(1,2) = 0;
       c1_forces(1,3) = (k * c1_properties(2) * c2_properties(2))/(c1_distance(1,3)^2);
       c1_forces(1,4) = (k * c1_properties(2) * h2_properties(2))/(c1_distance(1,4)^2);
       
       c2_forces = zeros(1,Total_Atoms);
       c2_forces(1,1) = (k * c2_properties(2) * h1_properties(2))/(c2_distance(1,1)^2);
       c2_forces(1,2) = (k * c2_properties(2) * c1_properties(2))/(c2_distance(1,2)^2);
       c2_forces(1,3) = 0;
       c2_forces(1,4) = (k * c2_properties(2) * h2_properties(2))/(c2_distance(1,4)^2);
       
       he_forces = zeros(He_Atoms,Total_Atoms);
       for l = 1:He_Atoms
           h1_forces(1,l+4) = (k * h1_properties(2) * he_properties(l,2))/(h1_distance(1,l+4)^2);
           h2_forces(1,l+4) = (k * h2_properties(2) * he_properties(l,2))/(h2_distance(1,l+4)^2);
           c1_forces(1,l+4) = (k * c1_properties(2) * he_properties(l,2))/(c1_distance(1,l+4)^2);
           c2_forces(1,l+4) = (k * c2_properties(2) * he_properties(l,2))/(c2_distance(1,l+4)^2);
           
           he_forces(l,1) = (k * he_properties(l,2) * h1_properties(2))/(he_distance(l,1)^2);
           he_forces(l,2) = (k * he_properties(l,2) * c1_properties(2))/(he_distance(l,2)^2);
           he_forces(l,3) = (k * he_properties(l,2) * c2_properties(2))/(he_distance(l,3)^2);
           he_forces(l,4) = (k * he_properties(l,2) * h2_properties(2))/(he_distance(l,4)^2);
           
           for p = 1:He_Atoms
               if he_distance(p,l+4) > 0
                    he_forces(p,l+4) = (k * he_properties(p,2) * he_properties(l,2))/(he_distance(p,l+4)^2);
               end
           end
       
       end
       
       %Split Forces into x,y,z and total.
       %x,y,z
       for p = 1:3 %switches through x, y and z.
           h1_properties(10+p) = (h1_forces(1,2) * ((h1_properties(2+p) - c1_properties(2+p))/(h1_distance(1,2)))) + ...
                               (h1_forces(1,3) * ((h1_properties(2+p) - c2_properties(2+p))/(h1_distance(1,3)))) + ...
                               (h1_forces(1,4) * ((h1_properties(2+p) - h2_properties(2+p))/(h1_distance(1,4))));

           h2_properties(10+p) = (h2_forces(1,1) * ((h2_properties(2+p) - h1_properties(2+p))/(h2_distance(1,1)))) + ...
                               (h2_forces(1,2) * ((h2_properties(2+p) - c1_properties(2+p))/(h2_distance(1,2)))) + ...
                               (h2_forces(1,3) * ((h2_properties(2+p) - c2_properties(2+p))/(h2_distance(1,3))));  

           c1_properties(10+p) = (c1_forces(1,1) * ((c1_properties(2+p) - h1_properties(2+p))/(c1_distance(1,1)))) + ...
                               (c1_forces(1,3) * ((c1_properties(2+p) - c2_properties(2+p))/(c1_distance(1,3)))) + ...
                               (c1_forces(1,4) * ((c1_properties(2+p) - h2_properties(2+p))/(c1_distance(1,4))));

           c2_properties(10+p) = (c2_forces(1,1) * ((c2_properties(2+p) - h1_properties(2+p))/(c2_distance(1,1)))) + ...
                               (c2_forces(1,2) * ((c2_properties(2+p) - c1_properties(2+p))/(c2_distance(1,2)))) + ...
                               (c2_forces(1,4) * ((c2_properties(2+p) - h2_properties(2+p))/(c2_distance(1,4))));                       

           for a = 1:He_Atoms %C2H2 - He interactions
                h1_properties(10+p) = h1_properties(10+p) + ...
                                    (h1_forces(1,a+4) * ((h1_properties(2+p) - he_properties(a,2+p))/(h1_distance(1,a+4))));
                h2_properties(10+p) = h2_properties(10+p) + ...
                                    (h2_forces(1,a+4) * ((h2_properties(2+p) - he_properties(a,2+p))/(h2_distance(1,a+4))));
                c1_properties(10+p) = c1_properties(10+p) + ...
                                    (c1_forces(1,a+4) * ((c1_properties(2+p) - he_properties(a,2+p))/(c1_distance(1,a+4))));
                c2_properties(10+p) = c2_properties(10+p) + ...
                                    (c2_forces(1,a+4) * ((c2_properties(2+p) - he_properties(a,2+p))/(c2_distance(1,a+4))));       

                he_properties(a,10+p) =  (he_forces(a,1) * ((he_properties(a,2+p) - h1_properties(2+p))/(he_distance(a,1)))) + ...
                                       (he_forces(a,2) * ((he_properties(a,2+p) - c1_properties(2+p))/(he_distance(a,2)))) + ...
                                       (he_forces(a,3) * ((he_properties(a,2+p) - c2_properties(2+p))/(he_distance(a,3)))) + ...
                                       (he_forces(a,4) * ((he_properties(a,2+p) - h2_properties(2+p))/(he_distance(a,4))));

                for b = 1:He_Atoms %He - He interactions.
                    if he_distance(a,b+4) > 0
                        he_properties(a,10+p) = he_properties(a,10+p) + ...
                                          (he_forces(a,b+4) * ((he_properties(a,2+p) - he_properties(b,2+p))/(he_distance(a,b+4))));
                    end

                end

           end
           
       end
       
       %Calculating x,y,z velocities
       % v(new) = v(old) + ((f(old) + f(curr))/2m)*dt
       
       for p = 1:3 %switches through x,y and z.
           
       h1_properties(6+p) = h1_properties(6+p) + ((h1_properties(13+p) + h1_properties(10+p))/(2*m_h))*dt;
       c1_properties(6+p) = c1_properties(6+p) + ((c1_properties(13+p) + c1_properties(10+p))/(2*m_c))*dt;
       c2_properties(6+p) = c2_properties(6+p) + ((c2_properties(13+p) + c2_properties(10+p))/(2*m_c))*dt;
       h2_properties(6+p) = h2_properties(6+p) + ((h2_properties(13+p) + h2_properties(10+p))/(2*m_h))*dt;
       
           for j = 1:He_Atoms
                he_properties(j,6+p) = he_properties(j,6+p) + ((he_properties(j,13+p) + he_properties(j,10+p))/(2*m_he))*dt;
           end
           
       end
       
       %Update total velocity.
       h1_properties(6) = sqrt(h1_properties(7)^2 + h1_properties(8)^2 + h1_properties(9)^2);
       c1_properties(6) = sqrt(c1_properties(7)^2 + c1_properties(8)^2 + c1_properties(9)^2);
       c2_properties(6) = sqrt(c2_properties(7)^2 + c2_properties(8)^2 + c2_properties(9)^2);
       h2_properties(6) = sqrt(h2_properties(7)^2 + h2_properties(8)^2 + h2_properties(9)^2);
       
       for b = 1:He_Atoms
            he_properties(b,6) = sqrt((he_properties(b,7))^2 + (he_properties(b,8))^2 + (he_properties(b,9))^2);
       end
       
       %Positions
       %New Positions
       % r(new) = r(old) + v(t)*dt + (f(t)/2m) * dt * dt
       
       for p = 1:3 %switch through x,y and z.
           
           h1_properties(2+p) = h1_properties(2+p) + h1_properties(6+p) * dt + ((h1_properties(10+p)) / (2*m_h))*dt*dt;
           c1_properties(2+p) = c1_properties(2+p) + c1_properties(6+p) * dt + ((c1_properties(10+p)) / (2*m_c))*dt*dt;
           c2_properties(2+p) = c2_properties(2+p) + c2_properties(6+p) * dt + ((c2_properties(10+p)) / (2*m_c))*dt*dt;
           h2_properties(2+p) = h2_properties(2+p) + h2_properties(6+p) * dt + ((h2_properties(10+p)) / (2*m_h))*dt*dt;
 
           for b = 1:He_Atoms
                he_properties(b,2+p) = he_properties(b,2+p) + he_properties(b,6+p) * dt + ((he_properties(b,10+p)) / (2*m_he))*dt*dt;
           end
       end
       
       %Update old forces for next loop
       
       for p = 1:3 %switch through x,y and z.
           
           h1_properties(13+p) = h1_properties(10+p);
           c1_properties(13+p) = c1_properties(10+p);
           c2_properties(13+p) = c2_properties(10+p);
           h2_properties(13+p) = h2_properties(10+p);

           for b = 1:He_Atoms
                he_properties(b,13+p) = he_properties(b,10+p);
           end
       end
       
       %calculate kinetic energy and convert from J to eV
       h1_properties(17) = 0.5 * h1_properties(1) * h1_properties(6) * h1_properties(6) * 6.2415e18;
       c1_properties(17) = 0.5 * c1_properties(1) * c1_properties(6) * c1_properties(6) * 6.2415e18;
       c2_properties(17) = 0.5 * c2_properties(1) * c2_properties(6) * c2_properties(6) * 6.2415e18;
       h2_properties(17) = 0.5 * h2_properties(1) * h2_properties(6) * h2_properties(6) * 6.2415e18;
       
       for b = 1:He_Atoms
          he_properties(b,17) = 0.5 * he_properties(b,1) * he_properties(b,6) * he_properties(b,6) * 6.2415e18; 
       end
       
       %result_test1(m) = h1_properties(6);
       %result_test2(m) = c1_properties(6);
       %result_test3(m) = he_properties(1,6);
    end
    
    %Calculate KE and convert to eV
    Total_Kinetic = h1_properties(17) + c1_properties(17) + c2_properties(17) + h2_properties(17);
    for b = 1:He_Atoms
       Total_Kinetic = Total_Kinetic + he_properties(b,17); 
    end
    
    %populate results
    result_initial_vector(i,:) = v;
    
    c1_vel = [c1_properties(7) ...
              c1_properties(8) ...
              c1_properties(9)];
    c2_vel = [c2_properties(7) ...
              c2_properties(8) ...
              c2_properties(9)];
    result_final_velocity_vector(i,1:3) = c1_vel;
    result_final_velocity_vector(i,4:6) = c2_vel;
    
end




