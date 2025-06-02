function a = ArrayResponse(angle,N,Chan,lambda)
    
    azimuth = angle;
    d = Chan.lambda/2;
    a = ULAArrayResponse(azimuth,N,d,lambda);
end

function a=ULAArrayResponse(azimuth,N,d,lamada)
    phi = azimuth(1);
    a=[(sqrt(1/N))*(exp(-1i*[0:N-1]*2*pi*d*sin(phi)/lamada)).'];
end
