function  at = far_field_manifold( Nt, f, fc, theta0 )
    c = 3e8;
    nn = -(Nt-1)/2:1:(Nt-1)/2;
    at = exp(1j*pi*nn*(f/fc)*sin(theta0))/sqrt(Nt);
end

