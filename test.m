clear all
clc
format long

%% Sampling Data
repeat = 30;

for i=1:repeat
    %error = sampling(-20, -30, 25, pi/4, 1.5*pi,8,1);
    % ith,unknown,method
    % unknown's error absolute value
    for j=1:6
        err(i,:,j) = sampling(-20, -30, 25, pi/4, 1.5*pi,8,j);
    end
end

%% Data classification
%row : method column : n_th data
x = reshape(err(:,1,:),[repeat,6])';
y = reshape(err(:,2,:),[repeat,6])';
z = reshape(err(:,3,:),[repeat,6])';
theta = reshape(err(:,4,:),[repeat,6])';
phi = reshape(err(:,5,:),[repeat,6])';

%column : method
x_mean = mean(x');
y_mean = mean(y');
z_mean = mean(z');
theta_mean = mean(theta');
phi_mean = mean(phi');

x_std = std(x');
y_std = std(y');
z_std = std(z');
theta_std = std(theta');
phi_std = std(phi');

outlierx = zeros(1,6);
outliery = zeros(1,6);
outlierz = zeros(1,6);
outliertheta = zeros(1,6);
outlierphi = zeros(1,6);

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;
f5 = figure;
for i=1:6
    %% Probability Density Function
    pdfx(i,:) = exp(-(x(i,:)-x_mean(i)).^2./(2*x_std(i)^2))./(x_std(i)*sqrt(2*pi));
    pdfy(i,:) = exp(-(y(i,:)-y_mean(i)).^2./(2*y_std(i)^2))./(y_std(i)*sqrt(2*pi));
    pdfz(i,:) = exp(-(z(i,:)-z_mean(i)).^2./(2*z_std(i)^2))./(z_std(i)*sqrt(2*pi));
    pdftheta(i,:) = exp(-(theta(i,:)-theta_mean(i)).^2./(2*theta_std(i)^2))./(theta_std(i)*sqrt(2*pi));
    pdfphi(i,:) = exp(-(phi(i,:)-phi_mean(i)).^2./(2*phi_std(i)^2))./(phi_std(i)*sqrt(2*pi));
    %subplot(1,5,1);
    
    %% Plot P.D.F.
    figure(f1);
    if mod(i,3)==1
        plot(x(i,:),pdfx(i,:),'x')
    elseif mod(i,3)==2
        plot(x(i,:),pdfx(i,:),'o')
    else
        plot(x(i,:),pdfx(i,:),'*')
    end
    hold on;
    
    
    figure(f2);
    if mod(i,3)==1
        plot(y(i,:),pdfy(i,:),'x')
    elseif mod(i,3)==2
        plot(y(i,:),pdfy(i,:),'o')
    else
        plot(y(i,:),pdfy(i,:),'*')
    end
    hold on;
    
    figure(f3);
    if mod(i,3)==1
        plot(z(i,:),pdfz(i,:),'x')
    elseif mod(i,3)==2
        plot(z(i,:),pdfz(i,:),'o')
    else
        plot(z(i,:),pdfz(i,:),'*')
    end
    hold on;
    
    figure(f4);
    if mod(i,3)==1
        plot(theta(i,:),pdftheta(i,:),'x')
    elseif mod(i,3)==2
        plot(theta(i,:),pdftheta(i,:),'o')
    else
        plot(theta(i,:),pdftheta(i,:),'*')
    end
    hold on;
    
    figure(f5);
    if mod(i,3)==1
        plot(phi(i,:),pdfphi(i,:),'x')
    elseif mod(i,3)==2
        plot(phi(i,:),pdfphi(i,:),'o')
    else
        plot(phi(i,:),pdfphi(i,:),'*')
    end
    hold on;
    
    %% if z-score >= 2.33 (1%) then outlier 
    for k=1:repeat
        if abs((x(i,k)-x_mean(i))/x_std(i)) >=2.33
            outlierx(i) = outlierx(i) + 1;
        end
        if abs((y(i,k)-y_mean(i))/y_std(i)) >=2.33
            outliery(i) = outliery(i) + 1;
        end
        if abs((z(i,k)-z_mean(i))/z_std(i)) >=2.33
            outlierz(i) = outlierz(i) + 1;
        end
        if abs((theta(i,k)-theta_mean(i))/theta_std(i)) >=2.33
            outliertheta(i) = outliertheta(i) + 1;
        end
        if abs((phi(i,k)-phi_mean(i))/phi_std(i)) >=2.33
            outlierphi(i) = outlierphi(i) + 1;
        end
    end
    %% Correlation Coefficient
    % row:unknown column:method
    corr_coefficient(:,i) = 100*[(x_std(i)/x_mean(i));(y_std(i)/y_mean(i));(z_std(i)/z_mean(i));(theta_std(i)/theta_mean(i));(phi_std(i)/phi_mean(i))];
end
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
%{
x_mean2 = reshape(mean(err(:,1,:)),[1,6]);
y_mean2 = reshape(mean(err(:,2,:)),[1,6]);
z_mean2 = reshape(mean(err(:,3,:)),[1,6]);
theta_mean2 = reshape(mean(err(:,4,:)),[1,6]);
phi_mean2 = reshape(mean(err(:,5,:)),[1,6]);

x_std2 = reshape(std(err(:,1,:)),[1,6]);
y_std2 = reshape(std(err(:,2,:)),[1,6]);
z_std2 = reshape(std(err(:,3,:)),[1,6]);
theta_std2 = reshape(std(err(:,4,:)),[1,6]);
phi_std2 = reshape(std(err(:,5,:)),[1,6]);
%}


