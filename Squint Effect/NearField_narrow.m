clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%

% figure;
% t = tiledlayout(1,1);
% t.TileSpacing = 'compact';
% t.Padding = 'compact';

%% fig 1
% nexttile
figure
hold on;
mesh(x_list_NN,y_list_NN,g_NN(:,:));
grid('off');
box on;
xlabel('x-axis [m]')
ylabel('y-axis [m]')
% xlabel(t,'{\emph{x}}-axis [m]', 'interpreter', 'latex', 'fontsize', 14)
% ylabel(t, '{\emph{y}}-axis [m]', 'interpreter', 'latex', 'fontsize', 14);

% colorbar
% colormap('Jet')