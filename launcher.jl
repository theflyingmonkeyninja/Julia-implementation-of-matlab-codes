
radb = 1; N=90;
thetab = 0:2*pi/N:2*pi-2*pi/N;
xb = radb*cos.(thetab);
yb = radb*sin.(thetab);

Nno = collect(1:length(xb));
Nodes =[Nno xb yb];

Nno1 = circshift(Nno,-1);

# connectivity matrix

ConM = [Nno Nno Nno1];

radf = 10; thetaf = 0;
xf = radf*cos(thetaf);
yf = radf*sin(thetaf);

FPnt = [xf yf];

w = 2*pi*collect(1:5:300);

using LinearAlgebra;
using SpecialFunctions;

kx = [0]; kx = transpose(kx);
kx = repeat(kx,inner=(length(w),1))

bemsolver=""; n = 3;  Ub = ones(length(xb),1); boundarypressure=""; route=[]; sampling ="static";
c=340; BEdomain = "";

using Profile

@profile (a,b)=AcBEMsolver(route,c,n,w,kx,Ub,FPnt,boundarypressure,BEdomain,sampling,ConM,Nodes);
