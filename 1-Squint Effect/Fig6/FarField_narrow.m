clc
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nt = 129;
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

r_center = [1, 1.5];
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
%%%%%%%%%%%%%%%

% figure;
% t = tiledlayout(1,1);
% t.TileSpacing = 'compact';
% t.Padding = 'compact';

%% fig 1
% nexttile
figure
hold on;
mesh(x_list_FN,y_list_FN,g_FN(:,:));
grid('off');
box on;
xlabel('x-axis [m]', Interpreter='latex', FontSize=8)
ylabel('y-axis [m]', Interpreter='latex',FontSize=8)
% xlabel(t,'{\emph{x}}-axis [m]', 'interpreter', 'latex', 'fontsize', 14)
% ylabel(t, '{\emph{y}}-axis [m]', 'interpreter', 'latex', 'fontsize', 14);

% colorbar
% colormap('Jet')