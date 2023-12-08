clear all
tic
cutoff=12;
sigma=3.5;
nsteps=round((cutoff/0.08)+1);
n = 200;
%% read simulation data
fid = fopen('coords_nvt_6.dump','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
fclose(fid);
box = cell2mat(box1);
xlo = box(1,1);xhi = box(1,2);
ylo = box(2,1);yhi = box(2,2);
zlo = box(3,1);zhi = box(3,2);
%%
delta_x=xhi-xlo;
delta_y=yhi-ylo;
fid = fopen('coords_nvt_6.dump','r');
for ii= 1:n
    C(ii,:) = textscan(fid,'%f %f %f %f %f %f','Headerlines',9);
end
%C contain 1-id 2-x 3-y 4-z 5-vx 6-vy 7-vz
fclose(fid);
ion_neoxy_dis=[];
nparticles=0;
%% find condensed ion
for t=1:n
   %t
   a1 = cell2mat(C(t,:));
   a21 = find(a1(:,3)==6);
   z_graft(t,1) = min(a1(a21,6));
   z_ends=a1(a21,6);
   z_brush(t,1) = min(z_ends(z_ends>z_graft(t,1)));
   a22 = find(a1(:,3)==9 & a1(:,6)<z_brush(t,1)& a1(:,6)>z_graft(t,1));
   a23 = find(a1(:,3)==5);
   for ii=1:length(a22)
       nparticles=nparticles+1;
       k=a22(ii);
       x_ion=a1(k,4);
       y_ion=a1(k,5);
       z_ion=a1(k,6);
       distance_min=linspace(0,0,length(a23));
       for jj=1:length(a23)
           w=a23(jj);
           x_water=a1(w,4);
           y_water=a1(w,5);
           z_water=a1(w,6);
           distance=sqrt((x_ion-x_water)^2+(y_ion-y_water)^2+(z_ion-z_water)^2);
           distance_xhi=sqrt((x_ion-(x_water+delta_x))^2+(y_ion-y_water)^2+(z_ion-z_water)^2);
           distance_xlo=sqrt((x_ion-(x_water-delta_x))^2+(y_ion-y_water)^2+(z_ion-z_water)^2);
           distance_yhi=sqrt((x_ion-x_water)^2+(y_ion-(y_water+delta_y))^2+(z_ion-z_water)^2);
           distance_ylo=sqrt((x_ion-x_water)^2+(y_ion-(y_water-delta_y))^2+(z_ion-z_water)^2);
           distance_xhi_yhi=sqrt((x_ion-(x_water+delta_x))^2+(y_ion-(y_water+delta_y))^2+(z_ion-z_water)^2);
           distance_xhi_ylo=sqrt((x_ion-(x_water+delta_x))^2+(y_ion-(y_water-delta_y))^2+(z_ion-z_water)^2);
           distance_xlo_yhi=sqrt((x_ion-(x_water-delta_x))^2+(y_ion-(y_water+delta_y))^2+(z_ion-z_water)^2);
           distance_xlo_ylo=sqrt((x_ion-(x_water-delta_x))^2+(y_ion-(y_water-delta_y))^2+(z_ion-z_water)^2);
           distance_con=[distance distance_xhi distance_xlo distance_yhi distance_ylo distance_xhi_yhi ...
               distance_xhi_ylo distance_xlo_yhi distance_xlo_ylo];
           distance_min(jj)=min(distance_con);
       end
       ion_neoxy_dis=vertcat(ion_neoxy_dis,min(distance_min));
   end
end
%% normalize
r=linspace(0,cutoff,nsteps);
q=linspace(0,0,nsteps-1);
for i=1:length(ion_neoxy_dis)
bracket=max(find(r(:)<ion_neoxy_dis(i)));
q(bracket)=q(bracket)+1;
end
cdf=linspace(0,0,length(q));
for b=1:length(q)
    cdf(b)=sum(q(1:b))/sum(q);
end
figure(1)
plot(r/sigma,vertcat(0,cdf'))

delta=cutoff/(nsteps-1);
f=zeros(nsteps-1,nsteps-1);
for k=1:nsteps-1
    for m=1:k
    f(k,m)=4*pi*(delta*m)^2*delta;
    end
end
cdf=cdf';
g1=f\cdf;
g=vertcat(0,g1);

check_cdf=zeros(length(r),1);
check_cdf(1)=0;
r=r';
for i=2:length(r)
    check_cdf(i)=trapz(r(1:i),4*pi*r(1:i).^2.*g(1:i));
end

figure(2)
plot(r/sigma,g*sigma^3)

figure (3)
plot(r,check_cdf,r,vertcat(0,cdf))

toc