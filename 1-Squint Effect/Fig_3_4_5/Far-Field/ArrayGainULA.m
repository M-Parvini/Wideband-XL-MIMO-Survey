function Array_Gain = ArrayGainULA(M, BFc, theta_F)
theta = -pi/2:0.00001:pi/2;  % AoA/AoD
nN = [-0.5, 0, 0.5];

ULA_Response = @(M,n,b,ang_z,ang_zF) ...
    1/M * abs(sin(((b*n+1)*sin(ang_z)-sin(ang_zF))*M*pi/2) ...
    /sin(((b*n+1)*sin(ang_z)-sin(ang_zF))*pi/2));

%% gain calculation for f = f_c
Array_Gain = zeros(length(nN), length(theta));
for i = 1:length(nN)
    for j=1:length(theta)
        Array_Gain(i, j) = ULA_Response(M, nN(i), BFc, theta(j), theta_F);
    end
end
end