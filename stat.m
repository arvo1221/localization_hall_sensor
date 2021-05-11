clear all
clc
format long

%% Init
repeat = 20;

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;

%Method 1 : Interior Point algorithm + 3axis
%Method 2 : Levenberg-Marquardt algorithm + 3axis
%Method 3 : Interior Point algorithm + 2axis
%Method 4 : Levenberg-Marquardt algorithm + 2axis
%Method 5 : Interior Point algorithm + 1axis
%Method 6 : Levenberg-Marquardt algorithm + 1axis


for j=1:6
    %% Sampling Data
    for i=1:repeat
        v(i,:) = [-20 -30 25 pi/4 1.5*pi];%mm mm mm rad rad
        sample(i,:) = sampling(v(i,1),v(i,2),v(i,3),v(i,4),v(i,5),8,j);%x y z theta phi sensor method
        error(i,:) = sample(i,:)-v(i,:); 
    end

    %x = error(:,1); y = error(:,2); z=error(:,3); theta = error(:,4); phi = error(:,5);
    means = [mean(error(:,1)) mean(error(:,2)) mean(error(:,3)) mean(error(:,4)) mean(error(:,5))];
    stds = [std(error(:,1)) std(error(:,2)) std(error(:,3)) std(error(:,4)) std(error(:,5))];
    m(j,:) = means; s(j,:) = stds;
    for i=1:5
        pdf(i,:) = (exp(-(error(:,i)-means(i)).^2./(2*stds(i)^2))./(stds(i)*sqrt(2*pi)))';
    end
    
    %% Plot data
    figure(f1);
    plot(error(:,1),pdf(1,:),'*')
    hold on;
    %plot(time,error(:,1));
    figure(f2);
    %plot(time,error(:,2));
    plot(error(:,2),pdf(2,:),'*')
    hold on;
    figure(f3);
    %plot(time,error(:,3));
    plot(error(:,3),pdf(3,:),'*')
    hold on;
    figure(f4);
    %plot(time,error(:,4));
    plot(error(:,4),pdf(4,:),'*')
    hold on;
    figure(f5);
    %plot(time,error(:,5));
    plot(error(:,5),pdf(5,:),'*')
    hold on;
    
end
%% pdf plot labeling

figure(f1);
legend('x1', 'x2', 'x3', 'x4', 'x5', 'x6')
title('x-axis')
xlabel('x-axis error')
ylabel('probability density function')
figure(f2);
legend('y1', 'y2', 'y3', 'y4', 'y5', 'y6')
title('y-axis')
xlabel('y-axis error')
ylabel('probability density function')
figure(f3);
legend('z1', 'z2', 'z3', 'z4', 'z5', 'z6')
title('z-axis')
xlabel('z-axis error')
ylabel('probability density function')
figure(f4);
legend('theta1', 'theta2', 'theta3', 'theta4', 'theta5', 'theta6')
title('theta')
xlabel('theta-axis error')
ylabel('probability density function')
figure(f5);
legend('phi1', 'phi2', 'phi3', 'phi4', 'phi5', 'phi6')
title('phi')
xlabel('phi-axis error')
ylabel('probability density function')

