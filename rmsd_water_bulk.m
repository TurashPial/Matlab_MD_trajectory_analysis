%----------water rmsd in the bulk --------
clear all
tic
n=1061;
fs1=0; 
fs2=3; %cutoff
f_per=10; %averaging
dt=50; %In fs
%% read simulation data
fid = fopen('coords_hbond_long_nvt_relax_5.dump','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
D(1,:) = textscan(fid,'%f %f %f %f %f %f','Headerlines',9);
fclose(fid);
box = cell2mat(box1);
xlo = box(1,1);xhi = box(1,2);
ylo = box(2,1);yhi = box(2,2);
zlo = box(3,1);zhi = box(3,2);

fid = fopen('coords_hbond_long_nvt_relax_5.dump','r');
for ii= 1:n
    C(ii,:) = textscan(fid,'%f %f %f %f %f %f','Headerlines',9);
end
%C contain 1-id 2-x 3-y 4-z 5-vx 6-vy 7-vz
fclose(fid);
n_origin=60;
dis=1001;
a3 = cell2mat(D(1,:));
a21 = find(a3(:,3)==6);
%% define matrix to check rotations
z_graft = min(a3(a21,6));
z_ends=a3(a21,6);
z_brush = min(z_ends(z_ends>z_graft));
a1=cell2mat(C(1,:));
a22 = find(a1(:,3)==7 & a1(:,6)>(-366) & a1(:,6)<(-340) & a1(:,4)<xhi-fs2 & a1(:,4)>xlo+fs2...
& a1(:,5)<yhi-fs2 & a1(:,5)>ylo+fs2);
a22=a22(1:3:end);
flag=0;

Hx_angle_store=zeros(n_origin*length(a22),length(1:10:dis));
Hy_angle_store=zeros(n_origin*length(a22),length(1:10:dis));
Hz_angle_store=zeros(n_origin*length(a22),length(1:10:dis));

Px_angle_store=zeros(n_origin*length(a22),length(1:10:dis));
Py_angle_store=zeros(n_origin*length(a22),length(1:10:dis));
Pz_angle_store=zeros(n_origin*length(a22),length(1:10:dis));

Qx_angle_store=zeros(n_origin*length(a22),length(1:10:dis));
Qy_angle_store=zeros(n_origin*length(a22),length(1:10:dis));
Qz_angle_store=zeros(n_origin*length(a22),length(1:10:dis));

Hx_angle_init=zeros(n_origin*length(a22),1);
Hy_angle_init=zeros(n_origin*length(a22),1);
Hz_angle_init=zeros(n_origin*length(a22),1);

Px_angle_init=zeros(n_origin*length(a22),1);
Py_angle_init=zeros(n_origin*length(a22),1);
Pz_angle_init=zeros(n_origin*length(a22),1);

Qx_angle_init=zeros(n_origin*length(a22),1);
Qy_angle_init=zeros(n_origin*length(a22),1);
Qz_angle_init=zeros(n_origin*length(a22),1);

angle_iter=0;
for t=1:n_origin
    t
    initial=cell2mat(C(t,:));
    for i=1:length(a22)
       angle_iter=angle_iter+1;
       o_index=a22(i);
       h1_index=o_index+1;
       h2_index=o_index+2;
       o_x_correction=0;
       o_y_correction=0;
       h1_x_correction=0;
       h1_y_correction=0;
       h2_x_correction=0;
       h2_y_correction=0;
       
       o_initial=[initial(o_index,4) initial(o_index,5) initial(o_index,6)];
       h1_initial=[initial(h1_index,4) initial(h1_index,5) initial(h1_index,6)];
       h2_initial=[initial(h2_index,4) initial(h2_index,5) initial(h2_index,6)];
       
       if max([norm(o_initial-h1_initial),norm(o_initial-h2_initial),norm(h1_initial-h2_initial)])>f_per
           flag=flag+1;
       end
       
       H_initial=o_initial-(h1_initial+h2_initial)/2;
       H_initial=H_initial/norm(H_initial);
       P_initial=h1_initial-h2_initial;
       P_initial=P_initial/norm(P_initial);
       Q_initial=cross((h1_initial-o_initial),(h2_initial-o_initial));
       Q_initial=Q_initial/norm(Q_initial);

       Hx_angle_init(angle_iter)=H_initial(1);
       Hy_angle_init(angle_iter)=H_initial(2);
       Hz_angle_init(angle_iter)=H_initial(3);

       Px_angle_init(angle_iter)=P_initial(1);
       Py_angle_init(angle_iter)=P_initial(2);
       Pz_angle_init(angle_iter)=P_initial(3);

       Qx_angle_init(angle_iter)=Q_initial(1);
       Qy_angle_init(angle_iter)=Q_initial(2);
       Qz_angle_init(angle_iter)=Q_initial(3);
           
       counter=0;
       for dis_loop=1:10:dis
           counter=counter+1;
           if dis_loop>1
           final_prev=cell2mat(C(t+dis_loop-10,:));
           elseif dis_loop==1
           final_prev=cell2mat(C(t+dis_loop-1,:));
           end
           final=cell2mat(C(t+dis_loop,:));
           
           if final(o_index,4)-final_prev(o_index,4)>(xhi-xlo)-f_per
              o_x_correction=o_x_correction+(xhi-xlo);
           elseif final(o_index,4)-final_prev(o_index,4)<-(xhi-xlo)+f_per    
              o_x_correction=o_x_correction-(xhi-xlo); 
           end
           if final(o_index,5)-final_prev(o_index,5)>(yhi-ylo)-f_per
              o_y_correction=o_y_correction+(yhi-ylo);
           elseif final(o_index,5)-final_prev(o_index,5)<-(yhi-ylo)+f_per    
              o_y_correction=o_y_correction-(yhi-ylo); 
           end
           
           if final(h1_index,4)-final_prev(h1_index,4)>(xhi-xlo)-f_per
              h1_x_correction=h1_x_correction+(xhi-xlo);
           elseif final(h1_index,4)-final_prev(h1_index,4)<-(xhi-xlo)+f_per    
              h1_x_correction=h1_x_correction-(xhi-xlo); 
           end
           if final(h1_index,5)-final_prev(h1_index,5)>(yhi-ylo)-f_per
              h1_y_correction=h1_y_correction+(yhi-ylo);
           elseif final(h1_index,5)-final_prev(h1_index,5)<-(yhi-ylo)+f_per    
              h1_y_correction=h1_y_correction-(yhi-ylo); 
           end
           
           if final(h2_index,4)-final_prev(h2_index,4)>(xhi-xlo)-f_per
              h2_x_correction=h2_x_correction+(xhi-xlo);
           elseif final(h2_index,4)-final_prev(h2_index,4)<-(xhi-xlo)+f_per    
              h2_x_correction=h2_x_correction-(xhi-xlo); 
           end
           if final(h2_index,5)-final_prev(h2_index,5)>(yhi-ylo)-f_per
              h2_y_correction=h2_y_correction+(yhi-ylo);
           elseif final(h2_index,5)-final_prev(h2_index,5)<-(yhi-ylo)+f_per    
              h2_y_correction=h2_y_correction-(yhi-ylo); 
           end
           
           o_final=[final(o_index,4)-o_x_correction final(o_index,5)-o_y_correction final(o_index,6)];
           h1_final=[final(h1_index,4)-h1_x_correction final(h1_index,5)-h1_y_correction final(h1_index,6)];
           h2_final=[final(h2_index,4)-h2_x_correction final(h2_index,5)-h2_y_correction final(h2_index,6)];
           
           H_final=o_final-(h1_final+h2_final)/2;
           H_final=H_final/norm(H_final);
           
           P_final=h1_final-h2_final;
           P_final=P_final/norm(P_final);
           
           Q_final=cross((h1_final-o_final),(h2_final-o_final));
           Q_final=Q_final/norm(Q_final);
           
           Hx_angle_store(angle_iter,counter)=H_final(1);
           Hy_angle_store(angle_iter,counter)=H_final(2);
           Hz_angle_store(angle_iter,counter)=H_final(3);
           
           Px_angle_store(angle_iter,counter)=P_final(1);
           Py_angle_store(angle_iter,counter)=P_final(2);
           Pz_angle_store(angle_iter,counter)=P_final(3);
           
           Qx_angle_store(angle_iter,counter)=Q_final(1);
           Qy_angle_store(angle_iter,counter)=Q_final(2);
           Qz_angle_store(angle_iter,counter)=Q_final(3);
           
       end
   end
end
%%
Hx_angle_rot=zeros(n_origin*length(a22),length(1:10:dis));
Hy_angle_rot=zeros(n_origin*length(a22),length(1:10:dis));
Hz_angle_rot=zeros(n_origin*length(a22),length(1:10:dis));

Px_angle_rot=zeros(n_origin*length(a22),length(1:10:dis));
Py_angle_rot=zeros(n_origin*length(a22),length(1:10:dis));
Pz_angle_rot=zeros(n_origin*length(a22),length(1:10:dis));

Qx_angle_rot=zeros(n_origin*length(a22),length(1:10:dis));
Qy_angle_rot=zeros(n_origin*length(a22),length(1:10:dis));
Qz_angle_rot=zeros(n_origin*length(a22),length(1:10:dis));

for k=1:n_origin*length(a22)
    k       
    for j=1:length(1:10:dis)
        for m=1:j
            if m==1
             H_initial=[Hx_angle_init(k) Hy_angle_init(k) Hz_angle_init(k)];   
             P_initial=[Px_angle_init(k) Py_angle_init(k) Pz_angle_init(k)];
             Q_initial=[Qx_angle_init(k) Qy_angle_init(k) Qz_angle_init(k)];
            elseif m>1
             H_initial=[Hx_angle_store(k,m-1) Hy_angle_store(k,m-1) Hz_angle_store(k,m-1)];
             P_initial=[Px_angle_store(k,m-1) Py_angle_store(k,m-1) Pz_angle_store(k,m-1)];
             Q_initial=[Qx_angle_store(k,m-1) Qy_angle_store(k,m-1) Qz_angle_store(k,m-1)];
            end
             H_final=[Hx_angle_store(k,m) Hy_angle_store(k,m) Hz_angle_store(k,m)];
             H_angle_rot=acos(dot(H_initial,H_final))*cross(H_initial,H_final)/norm(cross(H_initial,H_final));
             Hx_angle_rot(k,j)=Hx_angle_rot(k,j)+H_angle_rot(1);
             Hy_angle_rot(k,j)=Hy_angle_rot(k,j)+H_angle_rot(2);
             Hz_angle_rot(k,j)=Hz_angle_rot(k,j)+H_angle_rot(3);
             
             P_final=[Px_angle_store(k,m) Py_angle_store(k,m) Pz_angle_store(k,m)];
             P_angle_rot=acos(dot(P_initial,P_final))*cross(P_initial,P_final)/norm(cross(P_initial,P_final));
             Px_angle_rot(k,j)=Px_angle_rot(k,j)+P_angle_rot(1);
             Py_angle_rot(k,j)=Py_angle_rot(k,j)+P_angle_rot(2);
             Pz_angle_rot(k,j)=Pz_angle_rot(k,j)+P_angle_rot(3);
                          
             Q_final=[Qx_angle_store(k,m) Qy_angle_store(k,m) Qz_angle_store(k,m)];
             Q_angle_rot=acos(dot(Q_initial,Q_final))*cross(Q_initial,Q_final)/norm(cross(Q_initial,Q_final));
             Qx_angle_rot(k,j)=Qx_angle_rot(k,j)+Q_angle_rot(1);
             Qy_angle_rot(k,j)=Qy_angle_rot(k,j)+Q_angle_rot(2);
             Qz_angle_rot(k,j)=Qz_angle_rot(k,j)+Q_angle_rot(3);
        end
    end
end

 H_angle_mag=zeros(n_origin*length(a22),length(1:10:dis));
 P_angle_mag=zeros(n_origin*length(a22),length(1:10:dis));
 Q_angle_mag=zeros(n_origin*length(a22),length(1:10:dis));
for k=1:n_origin*length(a22)
    k
    for j=1:length(1:10:dis)
        H_angle_mag(k,j)=(norm([Hx_angle_rot(k,j) Hy_angle_rot(k,j) Hz_angle_rot(k,j)]))^2;
        P_angle_mag(k,j)=(norm([Px_angle_rot(k,j) Py_angle_rot(k,j) Pz_angle_rot(k,j)]))^2;
        Q_angle_mag(k,j)=(norm([Qx_angle_rot(k,j) Qy_angle_rot(k,j) Qz_angle_rot(k,j)]))^2; 
    end
end

H_rmsd=zeros(length(1:10:dis),1);
P_rmsd=zeros(length(1:10:dis),1);
Q_rmsd=zeros(length(1:10:dis),1);
for j=1:length(1:10:dis)
    j
    H_rmsd(j)=mean(H_angle_mag(:,j));
    P_rmsd(j)=mean(P_angle_mag(:,j));
    Q_rmsd(j)=mean(Q_angle_mag(:,j));
end
        
time_ps=linspace(dt,dis*dt,length(1:10:dis))/1000;
plot(time_ps,H_rmsd,'k-',time_ps,P_rmsd,'b-',time_ps,Q_rmsd,'r-')

save('rmsd_water_bulk_GD01_300K.mat','time_ps','H_rmsd','P_rmsd','Q_rmsd','flag')