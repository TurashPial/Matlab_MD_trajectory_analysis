clear all
tic
n = 100;
cutoff=8;
fs1=12;
fs2=20;
nsteps=round((cutoff/0.08)+1);
%% read data
fid = fopen('coords_nvt3.dump','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
fclose(fid);
box = cell2mat(box1);
xlo = box(1,1);xhi = box(1,2);
ylo = box(2,1);yhi = box(2,2);
zlo = box(3,1);zhi = box(3,2);

fid = fopen('coords_nvt3.dump','r');
for ii= 1:n
    C(ii,:) = textscan(fid,'%f %f %f %f %f %f','Headerlines',9);
end
%C contain id mol type x y z
fclose(fid);
%% calculate atom per radius shell
radius_shells=linspace(0,cutoff,nsteps);
particle_shell=linspace(0,0,nsteps);
nparticles=0;
for t=1:n
   t
   a1 = cell2mat(C(t,:));
   a21 = find(a1(:,3)==6);
   z_graft(t,1) = min(a1(a21,6));
   z_ends=a1(a21,6);
   z_brush(t,1) = min(z_ends(z_ends>z_graft(t,1)));
   a22 = find(a1(:,3)==5 & a1(:,6)<z_brush(t,1)-cutoff & a1(:,6)>z_graft(t,1)+cutoff & a1(:,4)<xhi-fs2 & a1(:,4)>xlo+fs2...
& a1(:,5)<yhi-fs2 & a1(:,5)>ylo+fs2);
   a23 = find(a1(:,3)==8 & a1(:,6)<z_brush(t,1) & a1(:,6)>z_graft(t,1) & a1(:,4)<xhi-fs1 & a1(:,4)>xlo+fs1...
& a1(:,5)<yhi-fs1 & a1(:,5)>ylo+fs1);
   for ii=1:length(a22)
       nparticles=nparticles+1;
       k=a22(ii);
       x_ion=a1(k,4);
       y_ion=a1(k,5);
       z_ion=a1(k,6);
       for jj=1:length(a23)
           w=a23(jj);
           x_water=a1(w,4);
           y_water=a1(w,5);
           z_water=a1(w,6);
           distance=sqrt((x_ion-x_water)^2+(y_ion-y_water)^2+(z_ion-z_water)^2);
           distance_xhi=sqrt((x_ion-(2*xhi-x_water))^2+(y_ion-y_water)^2+(z_ion-z_water)^2);
           distance_xlo=sqrt((x_ion-(2*xlo-x_water))^2+(y_ion-y_water)^2+(z_ion-z_water)^2);
           distance_yhi=sqrt((x_ion-x_water)^2+(y_ion-(2*yhi-y_water))^2+(z_ion-z_water)^2);
           distance_ylo=sqrt((x_ion-x_water)^2+(y_ion-(2*ylo-y_water))^2+(z_ion-z_water)^2);
           distance_xhi_yhi=sqrt((x_ion-(2*xhi-x_water))^2+(y_ion-(2*yhi-y_water))^2+(z_ion-z_water)^2);
           distance_xhi_ylo=sqrt((x_ion-(2*xhi-x_water))^2+(y_ion-(2*ylo-y_water))^2+(z_ion-z_water)^2);
           distance_xlo_yhi=sqrt((x_ion-(2*xlo-x_water))^2+(y_ion-(2*yhi-y_water))^2+(z_ion-z_water)^2);
           distance_xlo_ylo=sqrt((x_ion-(2*xlo-x_water))^2+(y_ion-(2*ylo-y_water))^2+(z_ion-z_water)^2);
           shell_number=max(find(radius_shells<distance));
           particle_shell(shell_number)=particle_shell(shell_number)+1;
           shell_number=max(find(radius_shells<distance_xhi));
           particle_shell(shell_number)=particle_shell(shell_number)+1;
           shell_number=max(find(radius_shells<distance_xlo));
           particle_shell(shell_number)=particle_shell(shell_number)+1;
           shell_number=max(find(radius_shells<distance_yhi));
           particle_shell(shell_number)=particle_shell(shell_number)+1;
           shell_number=max(find(radius_shells<distance_ylo));
           particle_shell(shell_number)=particle_shell(shell_number)+1;
           shell_number=max(find(radius_shells<distance_xhi_yhi));
           particle_shell(shell_number)=particle_shell(shell_number)+1;
           shell_number=max(find(radius_shells<distance_xhi_ylo));
           particle_shell(shell_number)=particle_shell(shell_number)+1;
           shell_number=max(find(radius_shells<distance_xlo_yhi));
           particle_shell(shell_number)=particle_shell(shell_number)+1;
           shell_number=max(find(radius_shells<distance_xlo_ylo));
           particle_shell(shell_number)=particle_shell(shell_number)+1;
       end
   end
end
%% normalize
nwater = length(find(a1(:,3)==7 & a1(:,6)<z_brush(n,1) & a1(:,6)>z_graft(n,1) & a1(:,4)<xhi-fs1 & a1(:,4)>xlo+fs1...
& a1(:,5)<yhi-fs1 & a1(:,5)>ylo+fs1));
rho=nwater/(xhi-xlo-2*fs1)/(yhi-ylo-2*fs1)/(z_brush(n,1)-z_graft(n,1));         

particle_shell=particle_shell/nparticles;
particle_shell_cutoff=particle_shell(1:nsteps-1);
for i=1:length(particle_shell_cutoff)
    neighbour(i)=sum(particle_shell_cutoff(1:i));
end
delta=cutoff/(nsteps-1);
f=zeros(nsteps-1,nsteps-1);
n1=neighbour(1:end);
for k=1:nsteps-1
    for m=1:k
    f(k,m)=4*pi*(delta*m)^2*rho*delta; 
    end
end
n1=n1';
g1=f\n1;
g=vertcat(0,g1);
figure (1)
plot(radius_shells,g)
check_n=zeros(length(radius_shells),1);
check_n(1)=0;
radius_shells=radius_shells';
for i=2:length(radius_shells)
    check_n(i)=trapz(radius_shells(1:i),4*pi*radius_shells(1:i).^2*rho.*g(1:i));
end
figure (2)
plot(radius_shells,check_n,radius_shells,vertcat(0,n1))
toc