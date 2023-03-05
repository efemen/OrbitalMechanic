function xX = t_xX(lat, sid)
% Topocentric to geocentric equatorial
    xX = zeros(3,3);
    xX(1,:) = [-sind(sid) -sind(lat)*cosd(sid) cosd(lat)*cosd(sid)];
    xX(2,:) = [cosd(sid) -sind(lat)*sind(sid) cosd(lat)*sind(sid)];
    xX(3,:) = [0 cosd(lat) sind(lat)];
end

