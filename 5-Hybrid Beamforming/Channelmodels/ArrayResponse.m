function a = ArrayResponse(angle,Entity,Chan,lambda)
    
    azimuth = angle;
    N = Entity.nAntenna;
    d = Chan.lambda/2;
    if N == 1
        a = 1;
    else
        switch Entity.AntennaType
            case 'ULA'
                a = ULAArrayResponse(azimuth,N,d,lambda);
            case 'UPA'
                a = UPAArrayResponse(azimuth,N,d,lambda);
        end
    end
end

function a=ULAArrayResponse(azimuth,N,d,lamada)
    phi = azimuth(1);
    a=[(sqrt(1/N))*(exp(-1i*[0:N-1]*2*pi*d*sin(phi)/lamada)).'];
end

function a=UPAArrayResponse(azimuth,N,d,lamada)
    
    N_H = sqrt(N);
    N_V = sqrt(N);
    % for UPA we assume equal number of elements in azimuth and elevation
    phi = azimuth(1);
    theta = azimuth(2);
    aH=[(sqrt(1/N_H))*(exp(-1i*[0:N_H-1]*2*pi*d*sin(phi)*cos(theta)/lamada)).'];
    aV=[(sqrt(1/N_V))*(exp(-1i*[0:N_V-1]*2*pi*d*sin(phi)/lamada)).'];
    a = kron(aV, aH);
end