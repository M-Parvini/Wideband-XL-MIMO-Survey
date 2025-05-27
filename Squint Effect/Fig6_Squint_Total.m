clc
clear

%%%%%%%%%%%% far-field narrow %%%%%%%%%%%%%%%%%%%%
Nt = 128;
Nr = 1;
fc = 8e9;
c = 3e8;
lambda_c = c/fc;
d = lambda_c / 2;

num_r = 1000;
r_start = 1;
r_end = 100;


x_list_FN = linspace(r_start, r_end, num_r);
y_list_FN = linspace(r_start, r_end, num_r);
% x_list = x.' * y;
% y_list = r_list.' * sin(theta_list);

g_FN = zeros(num_r, num_r);

r_center = [1, 1];
%%%%%%%%%%%%%%%
r = norm(r_center, 2);
theta = atan(r_center(1)/r_center(2));
f = fc;
H = far_field_manifold(Nt, f, fc, theta);
% w = narrow_focus(H,Nt);
w = H';
for idx_x = 1:num_r
    for idx_y = 1:num_r
        Theta = atan(y_list_FN(idx_y)/x_list_FN(idx_x));
        H = far_field_manifold(Nt, f, fc, Theta);
        g_FN(idx_x, idx_y)=  abs(H*w)/sqrt(Nt);
    end
end
%%%%%%%%%%%%%% far-field wide %%%%%%%%%%%%%%%%%
Nt = 128;
Nr = 1;
fc = 8e9;
B = 1e9;
f = [fc-B/2, fc, fc+B/2];
c = 3e8;
lambda_c = c/fc;
d = lambda_c / 2;

num_r = 1000;
r_start = 1;
r_end = 100;


x_list_FW = linspace(r_start, r_end, num_r);
y_list_FW = linspace(r_start, r_end, num_r);
% x_list = x.' * y;
% y_list = r_list.' * sin(theta_list);

g_FW = zeros(length(f), num_r, num_r);

r_center = [1, 1];
%%%%%%%%%%%%%%%
r = norm(r_center, 2);
theta = atan(r_center(1)/r_center(2));
H = far_field_manifold(Nt, fc, fc, theta);
% w = narrow_focus(H,Nt);
w = H';
for f_idx = 1:length(f)
    for idx_x = 1:num_r
        for idx_y = 1:num_r
            Theta = atan(y_list_FW(idx_y)/x_list_FW(idx_x));
            H = far_field_manifold(Nt, f(f_idx), fc, Theta);
            g_FW(f_idx, idx_x, idx_y)=  abs(H*w)/sqrt(Nt);
        end
    end
end

%%%%%%%%%%%%%%%% neardield narrow
Nt = 128;
Nr = 1;
fc = 8e9;
B = 1e9;
c = 3e8;
lambda_c = c/fc;
d = lambda_c / 2;

num_r = 1000;
r_start = 1;
r_end = 8;


x_list_NN = linspace(r_start, r_end, num_r);
y_list_NN = linspace(r_start, r_end, num_r);
% x_list = x.' * y;
% y_list = r_list.' * sin(theta_list);

g_NN = zeros(num_r, num_r);

r_center = [4, 4];
%%%%%%%%%%%%%%%
r = norm(r_center, 2);
theta = atan(r_center(1)/r_center(2));
H = near_field_manifold(Nt, d, fc, r, theta);
% w = narrow_focus(H,Nt);
w = H';
for idx_x = 1:num_r
    for idx_y = 1:num_r
        R = norm([x_list_NN(idx_x), y_list_NN(idx_y)], 2);
        Theta = atan(y_list_NN(idx_y)/x_list_NN(idx_x));
        H = near_field_manifold(Nt, d, fc, R, Theta);
        g_NN(idx_x, idx_y)=  abs(H*w)/sqrt(Nt);
    end
end

%%%%%%%%%%%%%%% nearfield wide %%%%%%%%%%%%%%%
Nt = 128;
Nr = 1;
fc = 8e9;
B = 1e9;
f = [fc-B/2, fc, fc+B/2];
c = 3e8;
lambda_c = c/fc;
d = lambda_c / 2;


num_r = 1000;
r_start = 1;
r_end = 8;


x_list_NW = linspace(r_start, r_end, num_r);
y_list_NW = linspace(r_start, r_end, num_r);
% x_list = x.' * y;
% y_list = r_list.' * sin(theta_list);

g_NW = zeros(length(f), num_r, num_r);

r_center = [4, 4];
%%%%%%%%%%%%%%%
r = norm(r_center, 2);
theta = atan(r_center(1)/r_center(2));
H = near_field_manifold(Nt, d, fc, r, theta);
w = H';
for f_idx = 1:length(f)
    for idx_x = 1:num_r
        for idx_y = 1:num_r
            R = norm([x_list_NW(idx_x), y_list_NW(idx_y)], 2);
            Theta = atan(y_list_NW(idx_y)/x_list_NW(idx_x));
            H = near_field_manifold(Nt, d, f(f_idx), R, Theta);
            g_NW(f_idx,idx_x, idx_y)=  abs(H*w)/sqrt(Nt);
        end
    end
end
%% plots
figure
subplot(2,2,4)
hold on
mesh(x_list_FN,y_list_FN,g_FN(:,:));
box on;
xlim([1 100])
ylim([1 100])
xlabel('x-axis [m]', Interpreter='latex', FontSize=6)
ylabel('y-axis [m]', Interpreter='latex',FontSize=6)

subplot(2,2,2)
for f_idx = 1:length(f)
    hold on;
    mesh(x_list_FW,y_list_FW,squeeze(g_FW(f_idx,:,:)));
    grid('off');
    box on;
end
xlim([1 100])
ylim([1 100])
xlabel('x-axis [m]', Interpreter='latex', FontSize=6)
ylabel('y-axis [m]', Interpreter='latex',FontSize=6)


subplot(2,2,3)
hold on
mesh(x_list_NN,y_list_NN,g_NN(:,:));
box on;
xlim([1 8])
ylim([1 8])
xlabel('x-axis [m]', Interpreter='latex', FontSize=6)
ylabel('y-axis [m]', Interpreter='latex',FontSize=6)

subplot(2,2,1)
for f_idx = 1:length(f)
    hold on;
    mesh(x_list_NW,y_list_NW,squeeze(g_NW(f_idx,:,:)));
    grid('off');
    box on;
end
xlim([1 8])
ylim([1 8])
xlabel('x-axis [m]', Interpreter='latex', FontSize=6)
ylabel('y-axis [m]', Interpreter='latex',FontSize=6)
