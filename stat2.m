clear all
clc
format long

%% Sampling Data
dt = 0.002;
repeat = 100;
time = 0:dt:(repeat-1)*dt;
f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;
f6 = figure;
for i=1:repeat
    v(i,:) = [-20*(dt*i)+50, 15*abs(cos(12*dt*i))-20, 2*(dt*i)+20, pi/2*abs(cos(dt*i)), pi*abs(cos(2*dt*i))];
    sample(i,:) = sampling(v(i,1),v(i,2),v(i,3),v(i,4),v(i,5),8,1);
    error(i,:) = sample(i,:)-v(i,:); 
end

%% Plot Result
figure(f1);
plot(time,sample(:,1));
hold on;
plot(time,v(:,1));
legend('estimated data', 'real data')
title('x-axis')
xlabel('time[s]')
ylabel('x [mm]')
hold off;

figure(f2);
plot(time,sample(:,2));
hold on;
plot(time,v(:,2));
legend('estimated data', 'real data')
title('y-axis')
xlabel('time[s]')
ylabel('y [mm]')

figure(f3);
plot(time,sample(:,3));
hold on;
plot(time,v(:,3));
legend('estimated data', 'real data')
title('z-axis')
xlabel('time[s]')
ylabel('z [mm]')

figure(f4);
plot(time,sample(:,4));
hold on;
plot(time,v(:,4));
legend('estimated data', 'real data')
title('theta')
xlabel('time[s]')
ylabel('theta [rad]')

figure(f5);
plot(time,sample(:,5));
hold on;
plot(time,v(:,5));
legend('estimated data', 'real data')
title('phi')
xlabel('time[s]')
ylabel('phi [rad]')

figure(f6)
for i=1:5
    plot(time,error(:,i));
    hold on;
end
legend('x','y','z','theta','phi')

