clear all
tic
n = 500;
GD=0.05;
N=49;
sigma=3.5;
%% read simulation data
fid = fopen('poly_4','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
fclose(fid);
box = cell2mat(box1);
xlo = box(1,1);xhi = box(1,2);
ylo = box(2,1);yhi = box(2,2);
zlo = box(3,1);zhi = box(3,2);

fid = fopen('poly_4','r');
for ii= 1:n
    C(ii,:) = textscan(fid,'%f %f %f %f %f %f','Headerlines',9);
end
%C contain contain 1-id 2-mol 3-type 4-x 5-y 6-z
fclose(fid);
%%
%z1 = 36;
ag = cell2mat(C(1,:));
a6 = find(ag(:,3)==8);
z_graft=min(ag(a6,6));
z_ends=ag(a6,6);
z_brush=max(z_ends(z_ends>z_graft));

nbins_inside = 10;  %how many bins in z direction to get the vel profile
nbins_outside = 5;
nbins=nbins_inside+nbins_outside;
nb1 = vertcat(linspace(z_graft,z_brush-10,nbins_inside+1)',linspace(z_brush-10,z_brush+10,nbins_outside+1)'); %divide z into bins
nb1(nbins_inside+1)=[];
%nb2=linspace(25,z1,12)';
%nb=vertcat(nb1,nb2);
nbcent = (nb1(1:nbins,1)+nb1(2:nbins+1,1))*0.5; %get centre of bins for plot
t=0;
num=0;
for jj = 1:n
   t=t+1
   a1 = cell2mat(C(jj,:));
   i=0;
   for kk = 1:nbins
       i=i+1;
       a21 = find(nb1(kk,1)<a1(:,6) & a1(:,6)<nb1(kk+1,1) & a1(:,3)==1);
       a22 = find(nb1(kk,1)<a1(:,6) & a1(:,6)<nb1(kk+1,1) & a1(:,3)==3);
       a23 = find(nb1(kk,1)<a1(:,6) & a1(:,6)<nb1(kk+1,1) & a1(:,3)==77);
       num1=length(a21);
       num2=length(a22);
       num3=length(a23);
       num(i,t)=num1+num2+num3;
   end
end

for i=1:nbins
    delta(i)=nb1(i+1)-nb1(i);
end
a13 = mean(num,2);
%% normalize and plot
delta=delta';
a13=a13./delta;
int=trapz(nbcent,a13)*10^(-10);
a13=a13/int*(N*GD/(sigma*10^(-10))^2);
monomer_dens=a13*(sigma*10^(-10))^3;
plot((nbcent-z_graft)/sigma,monomer_dens,'o')
toc

brush_height=trapz(nbcent,(nbcent-z_graft).*monomer_dens)/trapz(nbcent,monomer_dens);