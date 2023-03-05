function [ra, dec] = ECI2raDec(R)
% Topocentric to geocentric equatorial
    R_hat = R/norm(R);
    dec = asind(R_hat(3));
    if R(2) > 0
        ra = acosd(R_hat(1) / cosd(dec));
    else
        ra = 360 - acosd((R_hat(1) / cosd(dec)));
    end
end

