%Time Step
dt = 1e-15; %In seconds
Total_Steps = 1000;

%Constants
%Avogadro's Constant
avc = 6.0221e23;
%Atomic Mass Unit
amu = 1.6605e-27;
%Atomic Mass
am_c = 12.0107;
am_o = 15.9994;
%Mass in kg

%m_c = (am_c / avc) * 0.001;
%m_o = (am_o / avc) * 0.001;

m_c = mass('C') * amu;
m_o = mass('O') * amu;

%Coulomb's Constant
k = 8.9876e09;
%Elementary Charge
q = 1.602e-19;
%Charges on atoms
q_c = 1*q;
q_o = 1*q;

%Initial Conditions
%Position - Bond Length of CO = 112.8pm
%Assume parallel to x plane and centre of molecule is origin.
%Distance in Meters
r_c_x_curr = 0.498808e-10;
r_c_y_curr = -0.398808e-10;
r_o_x_curr = -0.298808e-10;
r_o_y_curr = 0.398808e-10;
%Velocity
%Assume stationary
v_c_x_curr = 0;
v_c_y_curr = 0;
v_o_x_curr = 0;
v_o_y_curr = 0;
%Force
%Assume 0 at t = 0
f_c_x_curr = 0;
f_c_y_curr = 0;
f_o_x_curr = 0;
f_o_y_curr = 0;

f_c_x_prev = 0;
f_c_y_prev = 0;
f_o_x_prev = 0;
f_o_y_prev = 0;


%Create Result Arrays
result_velocity_c_x = zeros(Total_Steps,1);
result_velocity_c_y = zeros(Total_Steps,1);
result_velocity_o_x = zeros(Total_Steps,1);
result_velocity_o_y = zeros(Total_Steps,1);
result_force_c_x = zeros(Total_Steps,1);
result_force_c_y = zeros(Total_Steps,1);
result_force_o_x = zeros(Total_Steps,1);
result_force_o_y = zeros(Total_Steps,1);
result_position_c_x = zeros(Total_Steps,1);
result_position_c_y = zeros(Total_Steps,1);
result_position_o_x = zeros(Total_Steps,1);
result_position_o_y = zeros(Total_Steps,1);

result_energy_c = zeros(Total_Steps,1);
result_energy_o = zeros(Total_Steps,1);

for i = 1:Total_Steps
   %Velocity Verlet Algorithm
   %1. Given current position and velocity, compute forces.
   %2. Update position
   %3. Repeat
   %Distance Between Atoms
   %Pythagoras
   r_c_o = sqrt((r_c_x_curr - r_o_x_curr)^2 + (r_c_y_curr - r_o_y_curr)^2);
   r_o_c = r_c_o;
   
   %Force Between Atoms
   %Coulombic Forces -> F = (k*q1*q2)/(r^2)
   f_c_o = (k * q_c * q_o) / ((r_c_o)^2);
   f_o_c = f_c_o;
   
   %Resolve Forces into x and y components
   %f_x = f * cos(theta)
   %f_y = f * sin(theta)
   %theta is angle between force and x direction
   f_c_o_x = f_c_o * (r_c_x_curr - r_o_x_curr) / r_c_o;
   f_c_o_y = f_c_o * (r_c_y_curr - r_o_y_curr) / r_c_o;
   
   f_o_c_x = f_o_c * (r_o_x_curr - r_c_x_curr) / r_o_c;
   f_o_c_y = f_o_c * (r_o_y_curr - r_c_y_curr) / r_o_c;
   
   
   
   f_c_x_curr = f_c_o_x;
   f_c_y_curr = f_c_o_y;
   f_o_x_curr = f_o_c_x;
   f_o_y_curr = f_o_c_y;
   
   %Velocities
   % v(new) = v(old) + (f/2m)*dt
   disp(string('Previous: ') + f_c_x_prev)
   disp(string('Current: ') + f_c_x_curr)
   v_c_x_curr = v_c_x_curr + ((f_c_x_curr + f_c_x_prev) / (2*m_c))*dt;
   v_c_y_curr = v_c_y_curr + ((f_c_y_curr + f_c_y_prev) / (2*m_c))*dt;
   %display(v_c_x_curr);
   
   v_o_x_curr = v_o_x_curr + ((f_o_x_curr + f_o_x_prev) / (2*m_o))*dt;
   v_o_y_curr = v_o_y_curr + ((f_o_y_curr + f_o_y_prev) / (2*m_o))*dt;
   
   f_c_x_prev = f_c_x_curr;
   f_c_y_prev = f_c_y_curr;
   f_o_x_prev = f_o_x_curr;
   f_o_y_prev = f_o_y_curr;
   
   v_c_curr = sqrt(v_c_x_curr^2 + v_c_y_curr^2);
   v_o_curr = sqrt(v_o_x_curr^2 + v_o_y_curr^2);
   %Positions
   %New Positions
   % r(new) = r(old) + v(t)*dt + (f(t)/2m) * dt * dt
   r_c_x_curr = r_c_x_curr + v_c_x_curr * dt + ((f_c_x_curr)/(2 * m_c))*dt*dt;
   r_c_y_curr = r_c_y_curr + v_c_y_curr * dt + ((f_c_y_curr)/(2 * m_c))*dt*dt;
   
   r_o_x_curr = r_o_x_curr + v_o_x_curr * dt + ((f_o_x_curr)/(2 * m_o))*dt*dt;
   r_o_y_curr = r_o_y_curr + v_o_y_curr * dt + ((f_o_y_curr)/(2 * m_o))*dt*dt;
   
   %Calculate kinetic energy and convert from J to eV
   ke_c_curr = 0.5 * m_c * v_c_curr * v_c_curr * 6.2415e18;
   ke_o_curr = 0.5 * m_o * v_o_curr * v_o_curr * 6.2415e18;
   
   %save all results in array
   result_velocity_c_x(i) = v_c_x_curr;
   result_velocity_c_y(i) = v_c_y_curr;
   result_velocity_o_x(i) = v_o_x_curr;
   result_velocity_o_y(i) = v_o_y_curr;
   result_force_c_x(i) = f_c_x_curr;
   result_force_c_y(i) = f_c_y_curr;
   result_force_o_x(i) = f_o_x_curr;
   result_force_o_y(i) = f_o_y_curr;
   result_position_c_x(i) = r_c_x_curr;
   result_position_c_y(i) = r_c_y_curr;
   result_position_o_x(i) = r_o_x_curr;
   result_position_o_y(i) = r_o_y_curr;
   result_energy_c(i) = ke_c_curr;
   result_energy_o(i) = ke_o_curr;
   
   
end
