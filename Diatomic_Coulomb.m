function Diatomic_Coulomb(atom1, atom1_x, atom1_y, atom1_z, atom1_charge, atom2, atom2_x, atom2_y, atom2_z, atom2_charge)

%Time Step
dt = 1e-15; %In seconds
total_steps = 1000;
%Constants
%Avogadro's Constant
avc = 6.0221e23; %#ok<NASGU>
%Atomic Mass Unit
amu = 1.6605e-27;
%Atomic Mass
am_1 = mass(atom1);
am_2 = mass(atom2);
%Mass in kg

%m_1 = (am_1 / avc) * 0.001;
%m_2 = (am_2 / avc) * 0.001;

m_1 = am_1 * amu;
m_2 = am_2 * amu;

%Coulomb's Constant
k = 8.9876e09;
%Elementary Charge
q = 1.602e-19;
%Charges on atoms
q_1 = atom1_charge * q;
q_2 = atom2_charge * q;

%Initial Conditions
%Position - Bond Length of CO = 112.8pm
%Assume parallel to x plane and centre of molecule is origin.
%Distance in Meters
r_1_x_curr = atom1_x;
r_1_y_curr = atom1_y;
r_1_z_curr = atom1_z;
r_2_x_curr = atom2_x;
r_2_y_curr = atom2_y;
r_2_z_curr = atom2_z;
%Velocity
%Assume stationary
v_1_x_curr = 0;
v_1_y_curr = 0;
v_1_z_curr = 0;
v_2_x_curr = 0;
v_2_y_curr = 0;
v_2_z_curr = 0;
%Force
%Assume 0 at t = 0
f_1_x_curr = 0;
f_1_y_curr = 0;
f_1_z_curr = 0;
f_2_x_curr = 0;
f_2_y_curr = 0;
f_2_z_curr = 0;

%Create Result Arrays
result_velocity_1_x = zeros(total_steps,1);
result_velocity_1_y = zeros(total_steps,1);
result_velocity_1_z = zeros(total_steps,1);
result_velocity_2_x = zeros(total_steps,1);
result_velocity_2_y = zeros(total_steps,1);
result_velocity_2_z = zeros(total_steps,1);
result_force_1_x = zeros(total_steps,1);
result_force_1_y = zeros(total_steps,1);
result_force_1_z = zeros(total_steps,1);
result_force_2_x = zeros(total_steps,1);
result_force_2_y = zeros(total_steps,1);
result_force_2_z = zeros(total_steps,1);
result_position_1_x = zeros(total_steps,1);
result_position_1_y = zeros(total_steps,1);
result_position_1_z = zeros(total_steps,1);
result_position_2_x = zeros(total_steps,1);
result_position_2_y = zeros(total_steps,1);
result_position_2_z = zeros(total_steps,1);

result_energy_1 = zeros(total_steps,1);
result_energy_2 = zeros(total_steps,1);

for i = 1:total_steps
   %Velocity Verlet Algorithm
   %1. Given current position and velocity, compute forces.
   %2. Update position
   %3. Repeat
   %Distance Between Atoms
   %Pythagoras
   r_1_2 = sqrt((r_1_x_curr - r_2_x_curr)^2 + (r_1_y_curr - r_2_y_curr)^2 + (r_1_z_curr - r_2_z_curr)^2);
   r_2_1 = r_1_2;
   
   %Force Between Atoms
   %Coulombic Forces -> F = (k*q1*q2)/(r^2)
   f_1_2 = (k * q_1 * q_2) / ((r_1_2)^2);
   f_2_1 = f_1_2;
   
   %Resolve Forces into x and y components
   %f_x = f * cos(theta)
   %f_y = f * sin(theta)
   %theta is angle between force and x direction
   f_1_2_x = f_1_2 * (r_1_x_curr - r_2_x_curr) / r_1_2;
   f_1_2_y = f_1_2 * (r_1_y_curr - r_2_y_curr) / r_1_2;
   f_1_2_z = f_1_2 * (r_1_z_curr - r_2_z_curr) / r_1_2;
   
   f_2_1_x = f_2_1 * (r_2_x_curr - r_1_x_curr) / r_2_1;
   f_2_1_y = f_2_1 * (r_2_y_curr - r_1_y_curr) / r_2_1;
   f_2_1_z = f_2_1 * (r_2_z_curr - r_2_z_curr) / r_2_1;
   
   
   
   f_1_x_curr = f_1_2_x;
   f_1_y_curr = f_1_2_y;
   f_1_z_curr = f_1_2_z;
   
   f_2_x_curr = f_2_1_x;
   f_2_y_curr = f_2_1_y;
   f_2_z_curr = f_2_1_z;
   
   %Velocities
   % v(new) = v(old) + (f/2m)*dt
   v_1_x_curr = v_1_x_curr + (f_1_x_curr / (2*m_1))*dt;
   v_1_y_curr = v_1_y_curr + (f_1_y_curr / (2*m_1))*dt;
   v_1_z_curr = v_1_z_curr + (f_1_z_curr / (2*m_1))*dt;
   %display(v_1_x_curr);
   
   v_2_x_curr = v_2_x_curr + (f_2_x_curr / (2*m_2))*dt;
   v_2_y_curr = v_2_y_curr + (f_2_y_curr / (2*m_2))*dt;
   v_2_z_curr = v_2_z_curr + (f_2_z_curr / (2*m_2))*dt;
   
   v_1_curr = sqrt(v_1_x_curr^2 + v_1_y_curr^2 + v_1_z_curr^2);
   v_2_curr = sqrt(v_2_x_curr^2 + v_2_y_curr^2 + v_2_z_curr^2);
   %Positions
   %New Positions
   % r(new) = r(old) + v(t)*dt + (f(t)/2m) * dt * dt
   r_1_x_curr = r_1_x_curr + v_1_x_curr * dt + ((f_1_x_curr)/(2 * m_1))*dt*dt;
   r_1_y_curr = r_1_y_curr + v_1_y_curr * dt + ((f_1_y_curr)/(2 * m_1))*dt*dt;
   r_1_z_curr = r_1_z_curr + v_1_z_curr * dt + ((f_1_z_curr)/(2 * m_1))*dt*dt;
   
   r_2_x_curr = r_2_x_curr + v_2_x_curr * dt + ((f_2_x_curr)/(2 * m_2))*dt*dt;
   r_2_y_curr = r_2_y_curr + v_2_y_curr * dt + ((f_2_y_curr)/(2 * m_2))*dt*dt;
   r_2_z_curr = r_2_z_curr + v_2_z_curr * dt + ((f_2_z_curr)/(2 * m_2))*dt*dt;
   
   %Calculate kinetic energy and convert from J to eV
   ke_1_curr = 0.5 * m_1 * v_1_curr * v_1_curr * 6.2415e18;
   ke_2_curr = 0.5 * m_2 * v_2_curr * v_2_curr * 6.2415e18;
   
   %save all results in array
   result_velocity_1_x(i) = v_1_x_curr;
   result_velocity_1_y(i) = v_1_y_curr;
   result_velocity_1_z(i) = v_1_z_curr;
   result_velocity_2_x(i) = v_2_x_curr;
   result_velocity_2_y(i) = v_2_y_curr;
   result_velocity_2_z(i) = v_2_z_curr;
   result_force_1_x(i) = f_1_x_curr;
   result_force_1_y(i) = f_1_y_curr;
   result_force_1_z(i) = f_1_z_curr;
   result_force_2_x(i) = f_2_x_curr;
   result_force_2_y(i) = f_2_y_curr;
   result_force_2_z(i) = f_2_z_curr;
   result_position_1_x(i) = r_1_x_curr;
   result_position_1_y(i) = r_1_y_curr;
   result_position_1_z(i) = r_1_z_curr;
   result_position_2_x(i) = r_2_x_curr;
   result_position_2_y(i) = r_2_y_curr;
   result_position_2_z(i) = r_2_z_curr;
   result_energy_1(i) = ke_1_curr;
   result_energy_2(i) = ke_2_curr;
   
   
end
   plot(result_energy_1)
   hold
   plot(result_energy_2)
end
