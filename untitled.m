%% Init
clear; clc;

car.v = 1;
car.m = 1;
car.x = 0;

th1 = -30+270; %deg
th2 = 30+270; %deg

l1.l = 1;
l3.l = 1;
l2.l = l1.l; l4.l = l3.l;



%% calc

th4 = (th2-th1)/2;
th3 = (th1+th2)/2;

d1 = l1.l*cosd(th4);
d2 = sqrt( abs(l3.l^2-l1.l^2*sind(th4)^2) );
L = d1+d2;

%% point masses

p1.x = car.x+l1.l*cosd(th1);
p1.y = l1.l*sind(th1);

p2.x = car.x+l2.l*cosd(th2);
p2.y = l2.l*sind(th2);

p3.x = car.x + L*cosd(th3);
p3.y = L*sind(th3);

%% Links

l1.x = car.x + l1.l/2*sind(th1);
l1.y = l1.l/2*cosd(th1);

l2.x = car.x + l2.l/2*sind(th2);
l2.y = l2.l/2*cosd(th2);

BN = sqrt(l3.l^2/4 - d2^2/4);
AN = d1+d2/2;
AB = sqrt( BN^2 +AN^2 );

beta2 = acosd(( l1.l^2 +l3.l^2/4-d2^2/4+(d1+d2/2)^2 -l3.l^2/4 )/(2*(l1.l)*(sqrt(l3.l^2))))

th5 = th1+beta2
l3.x = x+AB*sin(th5);
l3.y = 

figure(1)
clf;
plot(0,0,'o',p1.x,p1.y,'o',p2.x,p2.y,'o',p3.x,p3.y,'o');
grid on;
