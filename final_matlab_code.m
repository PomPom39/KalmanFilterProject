clc
clear all
close all

%Shape Of earth model
wgs84 = wgs84Ellipsoid('meters');

%Reading input data to test on from a csv file
data = csvread('Formatted_Dataset.csv', 1, 0);

%Extracting data nad storing it in seperate variables
Lat = data(:, 1);
Long = data(:, 2);
Heading = data(:, 3);
SOG = data(:, 5);

format long



%Time Step
dt = 5;

%Time samples
N = 100;


t = 1:N;


%Conversion Factor to convert from knots to m/s
CF = 0.5144


%Initial Parameters for the State Matrix
vel_x = 0;
vel_y = 0;
vel_z = 0;
lat = 76.54565;
lon = -68.9047;
ht = 80;


%Conversion to ECEF, ECEF Vectors, using libraries
[e_lat, e_lon, e_ht] = geodetic2ecef(wgs84, lat, lon, ht, 'degrees');
[e_vel_x, e_vel_y, e_vel_z] = enu2ecefv(vel_x, vel_y, vel_z, lat, lon);

%Initial State Matrix
x1(:, 1) = [e_lat; e_lon; e_ht; e_vel_x; e_vel_y; e_vel_z];

%A Matrix Definition (State Transition Model)

A = [ 1, 0, 0, dt, 0, 0;
      0, 1, 0, 0, dt, 0;
      0, 0, 1, 0, 0, dt;
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0];

%B Matrix(Control Input Model)
B = [0; 0; 0; e_vel_x; e_vel_y; e_vel_z];

%Q Matrix
Q = eye(6)*0.01;
R = eye(6)*1;
P = eye(6)*100;
H = eye(6);
I = eye(6);


for k = 1:N
    %Predicted value
    x2(: ,1) = A*x1(:, 1) + B*1;
    
    %Get Sensor Readings
    lat = Lat(k, :);
    lon = Long(k, :);
    ht = 80;
    speed = CF*SOG(k, :);
    angle = Heading(k, :);
    
    %Split Speed Component into ENU Components
    vel_x = speed*sind(angle);
    vel_y = speed*cosd(angle);
    vel_z = 0;
    
    %Convert to velocity ECEF reference frame, lat/lon to ECEF
    [e_lat, e_lon, e_ht] = geodetic2ecef(wgs84, lat, lon, ht, 'degrees');
    [e_vel_x, e_vel_y, e_vel_z] = enu2ecefv(vel_x, vel_y, vel_z, lat, lon);
    
    %Update The True Measurement Vector
    xt(:, 1) = [e_lat; e_lon; e_ht; e_vel_x; e_vel_y; e_vel_z;];
    
    %Update B Matrix (assumed to be constant)
    B = [0; 0; 0; e_vel_x; e_vel_y; e_vel_z;];
    
    %Calculate Covariance Matrix Values
    P = A*P*A' + Q;
    
    %Calculate Kalman Gain Values
    K = P*H'/(H*P*H' + R);
    
    
    %Corrected Value
    xp = x2 + K*(xt - H*x2);
    
    %Updating Covariance MAtrix
    P = (I - K*H)*P;
    
    x1 = xp;
   
    %Converting to co-ordinates to ECEF
    [k_lat, k_lon, k_ht] = ecef2geodetic(wgs84, x1(1, :), x1(2, :), x1(3, :), 'degrees');
    
    x_original(:, k) = [lat; lon; ht];
    x_predicted(:, k) = [k_lat; k_lon; k_ht];
    Kalman_Gain(:, k) = [K(1, 1);K(2, 2); K(3, 3); K(4, 4); K(5, 5);K(6, 6)];
    
   
end

%Plot graphs for data analysis


%Plot of Original Data vs Kalman Filtered Data
plot(x_original(1, :), x_original(2, :), 'b',x_predicted(1, :), x_predicted(2, :), 'r');
legend('Original',' Kalman Filtered')
xlabel('latitude')
ylabel('longitude')
title('Lat vs Lon')
grid on

%Plot Of Latitude vs Time
figure()
plot(t, x_original(1, :), 'b', t, x_predicted(1, :), 'r+');
legend('Original', 'Kalman Filtered');
xlabel('Time Step');
ylabel('Latitude');
title('Lat vs Time')
grid on

%Plot of Longitude vs Time
figure()
plot(t, x_original(2, :), 'b', t, x_predicted(2, :), 'r+');
legend('Original', 'Kalman Filtered');
xlabel('Time Step');
ylabel('Longitude');
title('Lon vs Time');
grid on

%Plot Of Kalman Gain vs Time
figure()
plot(t, Kalman_Gain(1, :))
xlabel('Time Step');
ylabel('Kalman Gain');
title('Kalman Gain Px');
grid on