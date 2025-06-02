clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%
% 
% figure;
% t = tiledlayout(1,1);
% t.TileSpacing = 'compact';
% t.Padding = 'compact';

%% fig 1
% nexttile

for f_idx = 1:length(f)
    hold on;
    mesh(x_list_FW,y_list_FW,squeeze(g_FW(f_idx,:,:)));
    grid('off');
    box on;
end
box on;
xlabel('x-axis [m]')
ylabel('y-axis [m]')
% xlabel(t,'{\emph{x}}-axis [m]', 'interpreter', 'latex', 'fontsize', 14)
% ylabel(t, '{\emph{y}}-axis [m]', 'interpreter', 'latex', 'fontsize', 14);
% 
% colorbar
% colormap('Jet')