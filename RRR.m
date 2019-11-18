%% RRR Robot

clear all;
l = [50, 30, 20];

f = 100;
dqmax = 1;
vmax = 1;
ddqmax = 1;
amax = 1;
param = 5/f;

% The following code is for FK and IK verification
% q = [pi/3, pi/4, -pi/2];
% disp(q);
% T = FK(q,l);
% disp(T);
% x = T(1,4);
% y = T(2,4);
% z = T(3,4);
% qi = IK(x, y, z, l);
% disp(qi);

%% Desired points
% P0 --PTP1--> P1 --LIN1--> P2 --PTP2--> P3 --LIN2--> P4

x0 = -15;
y0 = -20;
z0 = -10;
x1 = 15;
y1 = 20;
z1 = 50;
% LIN1: x = y/2 + 5;  z = 50
x2 = 20;
y2 = 30;
z2 = 50;
x3 = 10;
y3 = 20;
z3 = 10;
% LIN2: z = x;  y = 20
x4 = 5;
y4 = 20;
z4 = 5;

%% Trajectory planning

% PTP1
[t1, q1, dq1, ddq1] = PTP(x0, y0, z0, x1, y1, z1, l, dqmax, ddqmax, f);

T = FK(q1(:,1), l);
p1 = T(1:3,4);
for i = 2:length(q1)
    T = FK(q1(:,i), l);
    p1 = [p1, T(1:3,4)];
end
v1 = gradient(p1, 1/f);
a1 = gradient(v1, 1/f);

% LIN1
[t2, q2, dq2, ddq2] = PTP(x1, y1, z1, x2, y2, z2, l, vmax, amax, f);

T = FK(q2(:,1), l);
p2 = T(1:3,4);
for i = 2:length(q2)
    T = FK(q2(:,i), l);
    p2 = [p2, T(1:3,4)];
end
v2 = gradient(p2, 1/f);
a2 = gradient(v2, 1/f);

% Graphs
for i = 1:length(t2)
    t2(:,i) = t2(:,i) + t1(:,length(t1));
end
t = [t1, t2];
q = [q1, q2];
dq = [dq1, dq2];
ddq = [ddq1, ddq2];
p = [p1, p2];
v = [v1, v2];
a = [a1, a2];

figure(1)
plot(t, q(1,:), t, q(2,:), t, q(3,:));
xlabel('time, seconds');
ylabel('q, rads');
title('q1q2q3 trajectory plot');
figure(2)
plot(t, dq(1,:), t, dq(2,:), t, dq(3,:));
xlabel('time, seconds');
ylabel('q dot, rads/second');
title('q1q2q3 velocity plot');
figure(3)
plot(t, ddq(1,:), t, ddq(2,:), t, ddq(3,:));
xlabel('time, seconds');
ylabel('q dot dot, rads/second');
title('q1q2q3 acceleration plot');
figure(4)
plot(t, p(1,:), t, p(2,:), t, p(3,:));
title('XYZ trajectory plot');
figure(5)
plot(t, v(1,:), t, v(2,:), t, v(3,:));
title('XYZ velocity plot');
figure(6)
plot(t, a(1,:), t, a(2,:), t, a(3,:));
title('XYZ acceleration plot');


%% Junction




%% Functions

function T = FK(q,l)
    T = Rz(q(1))*Tz(l(1))*Rx(q(2))*Ty(l(2))*Rx(q(3))*Ty(l(3));
end

function q = IK(x,y,z,l)
    % upper part
    % elbow up configuration
    xi = sqrt(x^2 + y^2);
    z = z-l(1);
    q3 = atan2(real(-sqrt(1-((xi^2+z^2-l(2)^2-l(3)^2)/(2*l(2)*l(3)))^2)), ...
               (xi^2+z^2-l(2)^2-l(3)^2)/(2*l(2)*l(3)));
    phi1 = -atan2(real(l(3)*sin(q3)), real(l(2)+l(3)*cos(q3)));
    phi2 = atan2(z,xi);
    q2 = phi1 + phi2;
    % lower part
    q1 = atan2(-x,y);
    % full solution
    q = [real(q1), real(q2), real(q3)];
end

function [t, q, dq, ddq] = PTP(x0, y0, z0, x1, y1, z1, ...
                                                l, dqmax, ddqmax, f)
    q0 = IK(x0, y0, z0, l);
    q1 = IK(x1, y1, z1, l);
    delta_q = abs(q1 - q0);
    dqmax_ = sqrt(delta_q*ddqmax);
    trapezoid = zeros(1,3);
    t = zeros(1,3);
    T = zeros(1,3);
    dT = zeros(1,3);
    for i = 1:3
        if dqmax_(i) > dqmax
            trapezoid(i) = 1;
            t(i) = dqmax/ddqmax;
            T(i) = delta_q(i)/dqmax;
            dT(i) = T(i)-t(i);
        else
            t(i) = sqrt(abs(delta_q(i))/ddqmax);
        end
    end
    % multiple joints
    time = max(t);
    dTime = max(dT);
    % discretization
    discr = 1/f;
    t_n = ceil(time/discr)*discr;
    dT_n = ceil(dTime/discr)*discr;
    T_n = t_n + dT_n;
    ddqmax_n = delta_q/(T_n*t_n);
    t = linspace(0, 2*t_n+dT_n, (2*t_n+dT_n)*f+1);
    % linear space
    q = q0';
    qi = q0;
    dq = zeros(3,1);
    ddq = zeros(3,1);
    dq0 = zeros(1,3);
    for i = 1:(T_n + t_n)*f
        % figure out the part of the trapezoid
        if i*discr < t_n
            ddq1 = sign(q1 - q0).*ddqmax_n;
        elseif i*discr > T_n
            ddq1 = -sign(q1 - q0).*ddqmax_n;
        else
            ddq1 = zeros(1,3);
        end
        ddq = [ddq, ddq1'];
        dq_ = dq0 + ddq1*discr;
        dq = [dq, dq_'];
        q_ = qi + dq_*discr;
        q = [q, q_'];
        dq0 = dq_;
        qi = q_;
    end 
end

function J = J(q,l)
    Ti = FK(q,l);
    Ti = inv(Ti(1:3,1:3));
    Ti=[Ti,zeros(3,1);0 0 0 1];
    Td = Rzd(q(1))*Tz(l(1))*Rx(q(2))*Ty(l(2))*Rx(q(3))*Ty(l(3))*Ti;
    J1 = Jcol(Td);
    Td = Rz(q(1))*Tz(l(1))*Rxd(q(2))*Ty(l(2))*Rx(q(3))*Ty(l(3))*Ti;
    J2 = Jcol(Td);
    Td = Rz(q(1))*Tz(l(1))*Rx(q(2))*Ty(l(2))*Rxd(q(3))*Ty(l(3))*Ti;
    J3 = Jcol(Td);
    J = [J1, J2, J3];
end

function [t, q, dq, ddq] = LIN1(x0, y0, z0, x1, y1, z1, ...
                                l, vmax, amax, f)
    q1 = IK(x0, y0, z0, l);
    q2 = IK(x1, y1, z1, l);
    dp = [x1-x0, y1-y0, z1-z0]';
    vpmax = [vmax/sqrt(5), vxmax*2, 0]';
    apmax = [amax/sqrt(5), axmax*2, 0]';
    vpmax_ = [sqrt(abs(dx)*axmax), sqrt(abs(dy)*aymax), 0]';
    % in our case we can consider only one axis, because they are dependant
    if vpmax_(1) > vpmax(1)
        t = vpmax(1)/amax;
        T = abs(dp(1))/vpmax(1);
        dT = T-t; 
    else
        t = sqrt(abs(dx)/axmax);
    end
    % discretization
    discr = 1/f;
    t_n = ceil(t/discr)*discr;
    dT_n = ceil(dT/discr)*discr;
    T_n = t_n + dT_n;
%     vmax_n = dp/T_n;
    amax_n = dp/(T_n*t_n);

    % linspace
    p0 = [x0, y0, z0]';
    p = p0;
    v0 = [0, 0, 0]';
    v = v0;
    a = [0, 0, 0]';
    q = q1';
    dq = zeros(3,1);
    for i = 1:(T_n + t_n)*f
        if i*discr < t_n
            a1 = [sign(x1-x0), sign(y1-y0), sign(z1-z0)]'.*amax_n;
        elseif i*discr > T_n
            a1 = -[sign(x1-x0), sign(y1-y0), sign(z1-z0)]'.*amax_n;
        else
            a1 = [0, 0, 0]';
        end
        v1 = v0 + a1/f;
        % compute the inverse Jacobian
        Jac = J(q1,l);
        dq1 = Jac\v1;   % A\b == inv(A)*b
        dq = [dq, dq1];

        p1 = p0 + v1*discr;
        q1 = IK(p1(1), p1(2), p1(3), l);
        q = [q, q1'];
        p0 = p1;
        v0 = v1;
    end
    ddq = gradient(dq, 1/f);
    t = linspace(0, T_n+t_n, length(q));
end