clear all
tic
n=100;%Provide even number
fs1=5;
fs2=20;
f_per=10;
dt=1000;%In fs
fid = fopen('coords_nvt3.dump','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
D(1,:) = textscan(fid,'%f %f %f %f %f %f','Headerlines',9);
fclose(fid);
box = cell2mat(box1);
xlo = box(1,1);xhi = box(1,2);
ylo = box(2,1);yhi = box(2,2);
zlo = box(3,1);zhi = box(3,2);
%%
fid = fopen('coords_nvt3.dump','r');
for ii= 1:n
    C(ii,:) = textscan(fid,'%f %f %f %f %f %f','Headerlines',9);
end
%C contain 1-id 2-x 3-y 4-z 5-vx 6-vy 7-vz
fclose(fid);
n_origin=n/2;
a3 = cell2mat(D(1,:));
a21 = find(a3(:,3)==6);
z_graft = min(a3(a21,6));
z_ends=a3(a21,6);
z_brush = min(z_ends(z_ends>z_graft));
a1=cell2mat(C(1,:));
a22 = find(a1(:,3)==7 & a1(:,6)>z_graft+fs1+3 & a1(:,6)<z_brush-fs2 & a1(:,4)<xhi-fs2 & a1(:,4)>xlo+fs2...
& a1(:,5)<yhi-fs2 & a1(:,5)>ylo+fs2);
x_dis=linspace(0,0,n_origin);
y_dis=linspace(0,0,n_origin);
z_dis=linspace(0,0,n_origin);
%%
for t=1:n_origin
    t
    initial=cell2mat(C(t,:));
    for i=1:length(a22)
       p_index=a22(i);
       x_correction=0;
       y_correction=0;
       for dis_loop=1:n_origin
           final_prev=cell2mat(C(t+dis_loop-1,:));
           final=cell2mat(C(t+dis_loop,:));
           
           if final(p_index,4)-final_prev(p_index,4)>(xhi-xlo)-f_per
              x_correction=x_correction+(xhi-xlo);
           elseif final(p_index,4)-final_prev(p_index,4)<-(xhi-xlo)+f_per    
              x_correction=x_correction-(xhi-xlo); 
           end
           if final(p_index,5)-final_prev(p_index,5)>(yhi-ylo)-f_per
              y_correction=y_correction+(yhi-ylo);
           elseif final(p_index,5)-final_prev(p_index,5)<-(yhi-ylo)+f_per    
              y_correction=y_correction-(yhi-ylo); 
           end
           
           x_dis(dis_loop)=x_dis(dis_loop)+(final(p_index,4)-initial(p_index,4)-x_correction)^2;
           y_dis(dis_loop)=y_dis(dis_loop)+(final(p_index,5)-initial(p_index,5)-y_correction)^2;
           z_dis(dis_loop)=z_dis(dis_loop)+(final(p_index,6)-initial(p_index,6))^2;
       end
   end
end 
%% MSD
x_dis_av=x_dis/length(a22)/n_origin;
y_dis_av=y_dis/length(a22)/n_origin;
z_dis_av=z_dis/length(a22)/n_origin;
time=linspace(dt,n_origin*dt,n_origin);
plot(time,x_dis_av,'k-',time,y_dis_av,'b-',time,z_dis_av,'r-')

%% diffusion coefficient
Px = polyfit(time(round(n/4):end),x_dis_av(round(n/4):end),1);   
Dx=Px(1)/2*10^(-20)/10^(-15); %In m^2/s
Py = polyfit(time(round(n/4):end),y_dis_av(round(n/4):end),1);   
Dy=Py(1)/2*10^(-20)/10^(-15); %In m^2/s
Pz = polyfit(time(round(n/4):end),z_dis_av(round(n/4):end),1);   
Dz=Pz(1)/2*10^(-20)/10^(-15); %In m^2/s
Diff_constant=(Dx+Dy+Dz)/3
toc

total_dis_av=x_dis_av+y_dis_av+z_dis_av;
Pt=polyfit(time(25:end),total_dis_av(25:end),1);
Dt=Pt(1)/6*10^(-20)/10^(-15);