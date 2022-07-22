Constants
mu = 3.986*10^5; % km^3/s^2
Rearth = 6378; % km
rotrateEarth = 2*pi/86164; % rad/s
ObsLon = 262.25;
ObsLat = 30.3;


Problem 1: Validation
% Geostationary over west coast of Ecuador
inputsEcuador = [42164.1245,10^-6,10^-6,0,0,0,3,0,ObsLon,ObsLat];
[gTEcuador,~,bFEcuador,inertialEcuador,AzElEcuador,TEcuador] = groundTrackTool(inputsEcuador,mu,Rearth,rotrateEarth,'Geostationary Orbit');
% Molniya Orbit
inputsMolniya = [26561.73,0.722,63.4,270,80,0,3,0,ObsLon,ObsLat];
[gTMolniya,animMolniya,bFMolniya,inertialMolniya,AzElMolniya,TMolniya] = groundTrackTool(inputsMolniya,mu,Rearth,rotrateEarth,'Molniya Orbit');
% GPS orbit
inputsGPS = [26562.73,10^-6,55,-90,0,0,3,270,262.25,30.3];
[gTGPS, animGPS,bFGPS,inertialGPS,AzElGPS,TGPS] = groundTrackTool(inputsGPS,mu,Rearth,rotrateEarth,'GPS Orbit');

Problem 2: Examples
% case a
rinputA = [7000,0.01,30,0,0,0,3,0,ObsLon,ObsLat];
[gTA,animA,bFA,inertialA,AzElA,TA] = groundTrackTool(inputA,mu,Rearth,rotrateEarth,'Example A');

% case b
inputB = [7000,0.01,100,0,0,0,3,0,ObsLon,ObsLat];
[gTB,~,bFB,inertialB,AzElB,TB] = groundTrackTool(inputB,mu,Rearth,rotrateEarth,'Example B');

% case c
inputC = [17000 0.6 30 0 0 0 3 0 ObsLon ObsLat];
[gTC,~,bFC,inertialC,AzElC,TC] = groundTrackTool(inputC,mu,Rearth,rotrateEarth,'Example C');

% case d
inputD = [17000 0.6 100 0 0 0 3 0 ObsLon ObsLat];
[gTD,~,bFD,inertialD,AzElD,TD] = groundTrackTool(inputD,mu,Rearth,rotrateEarth,'Example D');

% case e
inputE = [13400 0.5 50 0 0 0 3 0 ObsLon ObsLat];
[gTE,~,bFE,inertialE,AzElE,TE] = groundTrackTool(inputE,mu,Rearth,rotrateEarth,'Example E');

% case f
inputF = [13400 0.5 50 90 0 0 3 0 ObsLon ObsLat];
[gTF,~,bFF,inertialF,AzElF,TF] = groundTrackTool(inputF,mu,Rearth,rotrateEarth,'Example F');

% case g
inputG = [13400 0.5 50 180 0 0 3 0 ObsLon ObsLat];
[gTG,~,bFG,inertialG,AzElG,TG] = groundTrackTool(inputG,mu,Rearth,rotrateEarth,'Example G');

% case h
inputH = [13400 0.5 50 270 0 0 3 0 ObsLon ObsLat];
[gTH,~,bFH,inertialH,AzElH,TH] = groundTrackTool(inputH,mu,Rearth,rotrateEarth,'Example H');


Problem 3: Repeat Ground Track
trepeat = 1*7*24*60*60; % s
i = 85; % degrees
e = 10^-6;
omega = 0;
LAN = 0;
nu0 = 0;
t_possible = 0:1:trepeat;
amin = 2*(300+Rearth);
amax = 2*(305+Rearth);
a_T = (mu*(t_possible/(2*pi)).^2).^(1/3);
for k = 1:length(a_T)
    if (a_T(k)>amin)&&(a_T(k)<amax)
        aRGT = a_T(k)
        break
    end
end

T = 2*pi*sqrt(aRGT^3/mu)

ratio_M_Nsc = rotrateEarth*sqrt(aRGT^3/mu);

M = rotrateEarth*T/(2*pi);
N = M/ratio_M_Nsc

inputsRepeat = [aRGT e i omega LAN nu0 3 0 ObsLon ObsLat];
[gTRepeat,animRepeat,bFRepeat,inertialRepeat,AzElRepeat,TRepeat] = groundTrackTool(inputsRepeat,mu,Rearth,rotrateEarth,'Repeat Ground Track');

[rRGT,~,~] = randv(aRGT,e,i,omega,LAN,nu0,1,mu);

[Xe,Ye,Ze] = sphere;
X2 = Xe*Rearth ;
Y2 = Ye*Rearth;
Z2 = Ze*Rearth;

bFRGT = bodyFixed(rRGT,X2,Y2,Z2,'Repeating Ground Track');
view(2)

h = aRGT/2 - Rearth
alpha = 2*asind(Rearth/(h+Rearth))


Problem 4: Existing Satellite Ground Track
% the ISS
rp = 412 + Rearth;
ra = 421 + Rearth;
a = 0.5*(rp+ra);
M0 = 63.8739;
[E0,~] = solveKepler_E(M0,0.0006808);
nu0 = 2*atand(sqrt((1+e)/(1-e))*tand(E0/2));

inputISS = [a 0.0006808 51.6421 97.6950 143.4753 nu0 3 0 ObsLon ObsLat];
[gTISS,animISS,bFISS,inertialISS,AzElISS,TISS] = groundTrackTool(inputISS,mu,Rearth,rotrateEarth,'ISS Orbit');



function dx = fcn(t,x)
mu = 3.986*10^5;
r = x(1:3);
rmag = norm(r);
v = x(4:6);
dx = zeros(6,1);
dx(1:3) = v;
dx(4:6) = -1*mu*r/rmag^3;
end

function R = R3(thet)
% thet is in degrees
c = cosd(thet);
s = sind(thet);
R = [c,s,0;-s,c,0;0,0,1];
end

function [r, v, t] = randv(a,e,i,omega,LAN,nu0,N,mu)
% calculates the r and v as a function of time for the orbit of a spacecraft around earth
% find r0, v0
p = a*(1-e^2);
r0mag = p/(1+e*cosd(nu0));
% perifocal frame
r0peri = [r0mag*cosd(nu0);r0mag*sind(nu0)];
% Rotation matrix
R11 = cosd(LAN)*cosd(omega) - sind(LAN)*sind(omega)*cosd(i);
R12 = -1*cosd(LAN)*sind(omega) - sind(LAN)*cosd(omega)*cosd(i);
R21 = sind(LAN)*cosd(omega) + cosd(LAN)*sind(omega)*cosd(i);
R22 = -1*sind(LAN)*sind(omega) + cosd(LAN)*cosd(omega)*cosd(i);
R31 = sind(omega)*sind(i);
R32 = cosd(omega)*sind(i);
R = [R11,R12;R21,R22;R31,R32];
% PQW to IJK frame
r0 = R*r0peri;
h = sqrt(mu*a*(1-e^2));
v0 = mu/h*[-R11*sind(nu0)+R12*(e+cosd(nu0));-R21*sind(nu0)+R22*(e+cosd(nu0));-R31*sind(nu0)+R32*(e+cosd(nu0))];
% calculate tf from given initial conditions
tf = N*2*pi*a^(3/2)/sqrt(mu);
tspan = [0:1:tf];
% initial state
x0 = [r0;v0];
x0 = x0';
% ode45 to solve for r(t) and v(t) from t0 to tf
[t,X] = ode113(@(t,x)fcn(t,x),tspan,x0);
r = X(:,1:3);
v = X(:,4:6);
end

function [long, lat] = longlat(r,t,thetaG0,rotrateEarth)
long = zeros(length(t),1);
lat = zeros(length(t),1);

for l = 1:length(t)
    thetaGST = rad2deg(thetaG0*pi/180 + rotrateEarth*l);
    rB = R3(thetaGST)*r(l,:)';
    rBnorm = norm(rB);
    phi = asind(rB(3)/rBnorm);
    lambda = atan2d(rB(2),rB(1));
    lat(l,1) = phi;
    long(l,1) = lambda; 
    
end

end

function [fig1, fig2] = gT(long,lat,caseTitle,Az,El)
% find visible lat/long
visLong = [];
visLat = [];
for k = 1:length(Az)
    if (Az(k,1)>0)&&(El(k,1)<90)
        visLong(k,1) = long(k,1);
        visLat(k,1) = lat(k,1);
    elseif (Az(k,1)>0)&&(El(k,1)>-90)
        visLong(k,1) = long(k,1);
        visLat(k,1) = lat(k,1);
    else 
        continue
    end
end
% plot ground track
fig1 = figure;
geoshow('landareas.shp','FaceColor',[0.5 1 0.5]);
hold on
grid on
plot(long,lat,'.r')
plot(visLong,visLat,'.b')

title(caseTitle)
xlabel('Longitude (\lambda, deg)')
ylabel('Latitude (\phi, deg)')

% ground track animation
fig2 = figure;
geoshow('landareas.shp','FaceColor',[0.5 1 0.5]);
hold on 
grid on
plot(long,lat,'.r')
title(caseTitle)
xlabel('Longitude (\lambda, deg)')
ylabel('Latitude (\phi, deg)')
p = plot(long(1),lat(1),'o','MarkerFaceColor','b');
hold off
axis manual

for k = 2:length(long)
    p.XData = long(k);
    p.YData = lat(k);
    drawnow limitrate
   
end
end

function [orbitBF] = bodyFixed(r,X2,Y2,Z2,caseTitle)
orbitBF = figure;
surf(X2,Y2,Z2)
hold on
plot3(r(:,1),r(:,2),r(:,3),'r')
axis equal
title(caseTitle)
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')

end

function [orbitInertial] = inertial(r,lat,long,X2,Y2,Z2,caseTitle)
xI = zeros(1,length(r));
yI = zeros(1,length(r));
zI = zeros(1,length(r));
for k = 1:length(r)
    rnorm = norm(r(k,:));
    phi = lat(k);
    lambda = long(k);
    xI(1,k) = rnorm*cosd(lambda)*cosd(phi);
    yI(1,k) = rnorm*sind(lambda)*cosd(phi);
    zI(1,k) = rnorm*sind(phi);
    
end

orbitInertial = figure;
surf(X2,Y2,Z2)
hold on 
plot3(xI,yI,zI,'r')
axis equal
title(caseTitle)
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')

end

function [Az,El] = AzEl(long,ObsLon,ObsLat)
Az = zeros(length(long),1);
El = zeros(length(long),1);
for k = 1:length(long)
	lambdas = long(k,1);
	lambda = ObsLon - lambdas;
	Az(k,1) = atan2d(tand(lambda),sind(ObsLat));
	El(k,1) = atan2d((cosd(ObsLat)*cosd(lambda)-0.151),sqrt(1-cosd(ObsLat)^2*cosd(lambda)^2));
end
end

function [AzElplot] = plotAzEl(Az,El,t,caseTitle)
% plot Az and El with time
AzElplot = figure;
plot(t',Az',t',El')
title(caseTitle)
xlabel('Time (s)')
legend('Azimuth (degrees)','Elevation (degrees)')

end

function [figGT,animGT,figBF,figInertial,figAzEl,T] = groundTrackTool(inputs,mu,Rearth,rotrateEarth,caseTitle)
% input vector
% [a, e, i, omega, LAN, nu0, N, thetaG0, ObsLon, ObsLat]
% a = semi-major axis (km)
% e = eccentricity
% i = inclination (deg)
% omega = argument of periapse (deg)
% LAN = longitude of the ascending node (deg)
% nu0 = initial true anomaly (deg)
% N = number of complete spacecraft orbits to propagate
% thetaG0 = greenwich sidereal angle at t = t0 (deg)
% ObsLon = observer longitude (deg)
% ObsLat = observer latitude (deg)

a = inputs(1);
e = inputs(2);
i = inputs(3);
omega = inputs(4);
LAN = inputs(5);
nu0 = inputs(6);
N = inputs(7);
thetaG0 = inputs(8);
ObsLon = inputs(9);
ObsLat = inputs(10);

[r,v,t] = randv(a,e,i,omega,LAN,nu0,N,mu);

r0 = r(1,:);
v0 = v(1,:);
rfinal = r(end,:);
vfinal = v(end,:);

[long,lat] = longlat(r,t,thetaG0,rotrateEarth);
long0 = long(1);
lat0 = lat(1);
longfinal = long(end);
latfinal = lat(end);

% calculate azimuth and elevation as a function of time for the observer
[Az, El] = AzEl(long,ObsLon,ObsLat);
Az0 = Az(1);
El0 = El(1);
Azfinal = Az(end);
Elfinal = El(end);

% plot and animate ground track
[figGT, animGT] = gT(long,lat,caseTitle,Az,El);

% earth ref
[Xe,Ye,Ze] = sphere;
X2 = Xe*Rearth ;
Y2 = Ye*Rearth;
Z2 = Ze*Rearth;

% body fixed frame orbit
[figBF] = bodyFixed(r,X2,Y2,Z2,caseTitle);

% inertial frame orbit
[figInertial] = inertial(r,lat,long,X2,Y2,Z2,caseTitle);

% plot Az and El w time
[figAzEl] = plotAzEl(Az,El,t,caseTitle);

% table of beginning and end values
time = {'Beginning';'End'};
Azimuth = {Az0;Azfinal};
Elevation = {El0;Elfinal};
InertialPosition = {r0;rfinal};
InertialVelocity = {v0;vfinal};
Longitude = {long0;longfinal};
Latitude = {lat0;latfinal};
T = table(time,InertialPosition,InertialVelocity,Azimuth,Elevation,Longitude,Latitude);
T.Properties.VariableNames = {'Time During Simulation','Inertial Position in the Body-Fixed Frame (10^3 km)','Inertial Velocity in the Body-Fixed Frame (km)','Azimuth (degrees)','Elevation (degrees)','Longitude (degrees)','Latitude (degrees)'}
end
