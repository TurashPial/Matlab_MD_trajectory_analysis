clear all
tic
fs2=4;
n = 4000;
oo_crit=3.5; %H-bond defination
oh_crit=1000;
angle_crit=30;
fs3=oo_crit;
f_per=10;
%% read simulation data
fid = fopen('coords_hbond_long_nvt_relax_5.dump','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
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
%%
t1 = cell2mat(C(1,:));
a6 = find(t1(:,3)==6);
z_graft = min(t1(a6,6));
z_ends=t1(a6,6);
z_brush = min(z_ends(z_ends>z_graft));
H_bond_list=zeros([],4);

for t=2:41
   t
   a1 = cell2mat(C(t,:));
   a2 = cell2mat(C(t-1,:));
  
   a21 = find(a1(:,3)==7 & a1(:,6)>(-366) & a1(:,6)<(-340) & a1(:,4)<xhi-fs2 & a1(:,4)>xlo+fs2 ...
       & a1(:,5)<yhi-fs2 & a1(:,5)>ylo+fs2);
   
   for k=1:length(a21)
       x_o=a1(a21(k),4);
       y_o=a1(a21(k),5);
       z_o=a1(a21(k),6);
       x_h1=a1(a21(k)+1,4);
       y_h1=a1(a21(k)+1,5);
       z_h1=a1(a21(k)+1,6);
       x_h2=a1(a21(k)+2,4);
       y_h2=a1(a21(k)+2,5);
       z_h2=a1(a21(k)+2,6);
       
       for m=1:length(a21)
           x_o_ne=a1(a21(m),4);
           y_o_ne=a1(a21(m),5);
           z_o_ne=a1(a21(m),6);
           oo_dis=sqrt((x_o-x_o_ne)^2+(y_o-y_o_ne)^2+(z_o-z_o_ne)^2); 
                if oo_dis==0 || oo_dis>oo_crit
                     continue
                end
       x_h1_ne=a1(a21(m)+1,4);
       y_h1_ne=a1(a21(m)+1,5);
       z_h1_ne=a1(a21(m)+1,6);
       x_h2_ne=a1(a21(m)+2,4);
       y_h2_ne=a1(a21(m)+2,5);
       z_h2_ne=a1(a21(m)+2,6);
       
       o_h1ne_dis=sqrt((x_o-x_h1_ne)^2+(y_o-y_h1_ne)^2+(z_o-z_h1_ne)^2);
       u_o_one_h1ne=[(x_o-x_o_ne) (y_o-y_o_ne) (z_o-z_o_ne)];
       v_o_one_h1ne=[(x_h1_ne-x_o_ne) (y_h1_ne-y_o_ne) (z_h1_ne-z_o_ne)];
       angle_o_one_h1ne=acosd(dot(u_o_one_h1ne,v_o_one_h1ne)/(norm(u_o_one_h1ne)*norm(v_o_one_h1ne))); %hbond angle checking
  
       if o_h1ne_dis<oh_crit & angle_o_one_h1ne<angle_crit
           o_prev=[a2(a21(k),4) a2(a21(k),5) a2(a21(k),6)];
           o_ne_prev=[a2(a21(m),4) a2(a21(m),5) a2(a21(m),6)];
           h_ne_prev=[a2(a21(m)+1,4) a2(a21(m)+1,5) a2(a21(m)+1,6)];
           h_bond_prev=check_hbond(o_prev,o_ne_prev,h_ne_prev);
           if h_bond_prev==0
               H_bond_list=vertcat(H_bond_list,[a21(k) a21(m) a21(m)+1 t]);
           end
       end
       
       o_h2ne_dis=sqrt((x_o-x_h2_ne)^2+(y_o-y_h2_ne)^2+(z_o-z_h2_ne)^2);
       u_o_one_h2ne=[(x_o-x_o_ne) (y_o-y_o_ne) (z_o-z_o_ne)];
       v_o_one_h2ne=[(x_h2_ne-x_o_ne) (y_h2_ne-y_o_ne) (z_h2_ne-z_o_ne)];
       angle_o_one_h2ne=acosd(dot(u_o_one_h2ne,v_o_one_h2ne)/(norm(u_o_one_h2ne)*norm(v_o_one_h2ne)));
       
       if o_h2ne_dis<oh_crit & angle_o_one_h2ne<angle_crit
           o_prev=[a2(a21(k),4) a2(a21(k),5) a2(a21(k),6)];
           o_ne_prev=[a2(a21(m),4) a2(a21(m),5) a2(a21(m),6)];
           h_ne_prev=[a2(a21(m)+2,4) a2(a21(m)+2,5) a2(a21(m)+2,6)];
           h_bond_prev=check_hbond(o_prev,o_ne_prev,h_ne_prev);
           if h_bond_prev==0
               H_bond_list=vertcat(H_bond_list,[a21(k) a21(m) a21(m)+2 t]);
           end
       end
       
       h1_one_dis=sqrt((x_h1-x_o_ne)^2+(y_h1-y_o_ne)^2+(z_h1-z_o_ne)^2);
       u_h1_o_one=[(x_h1-x_o) (y_h1-y_o) (z_h1-z_o)];
       v_h1_o_one=[(x_o_ne-x_o) (y_o_ne-y_o) (z_o_ne-z_o)];
       angle_h1_o_one=acosd(dot(u_h1_o_one,v_h1_o_one)/(norm(u_h1_o_one)*norm(v_h1_o_one)));
       
       if h1_one_dis<oh_crit &  angle_h1_o_one<angle_crit
           o_prev=[a2(a21(m),4) a2(a21(m),5) a2(a21(m),6)];
           o_ne_prev=[a2(a21(k),4) a2(a21(k),5) a2(a21(k),6)];
           h_ne_prev=[a2(a21(k)+1,4) a2(a21(k)+1,5) a2(a21(k)+1,6)];
           h_bond_prev=check_hbond(o_prev,o_ne_prev,h_ne_prev);
           if h_bond_prev==0
               H_bond_list=vertcat(H_bond_list,[a21(m) a21(k) a21(k)+1 t]);
           end
       end
       
       h2_one_dis=sqrt((x_h2-x_o_ne)^2+(y_h2-y_o_ne)^2+(z_h2-z_o_ne)^2);
       u_h2_o_one=[(x_h2-x_o) (y_h2-y_o) (z_h2-z_o)];
       v_h2_o_one=[(x_o_ne-x_o) (y_o_ne-y_o) (z_o_ne-z_o)];
       angle_h2_o_one=acosd(dot(u_h2_o_one,v_h2_o_one)/(norm(u_h2_o_one)*norm(v_h2_o_one)));
       
       if h2_one_dis<oh_crit &  angle_h2_o_one<angle_crit
           o_prev=[a2(a21(m),4) a2(a21(m),5) a2(a21(m),6)];
           o_ne_prev=[a2(a21(k),4) a2(a21(k),5) a2(a21(k),6)];
           h_ne_prev=[a2(a21(k)+2,4) a2(a21(k)+2,5) a2(a21(k)+2,6)];
           h_bond_prev=check_hbond(o_prev,o_ne_prev,h_ne_prev);
           if h_bond_prev==0
               H_bond_list=vertcat(H_bond_list,[a21(m) a21(k) a21(k)+2 t]);
           end
       end
       end    
   end

end

%H_bond_list=H_bond_list(1:4:end,:);
o_a_x_shift=zeros(length(H_bond_list),1);
o_a_y_shift=zeros(length(H_bond_list),1);
o_d_x_shift=zeros(length(H_bond_list),1);
o_d_y_shift=zeros(length(H_bond_list),1);
hyd_x_shift=zeros(length(H_bond_list),1);
hyd_y_shift=zeros(length(H_bond_list),1);

%time_calc=vertcat(linspace(50,1000,20)',linspace(1500,198500,395)')/50;
time_calc=vertcat(linspace(50,1000,20)',linspace(1500,100000,198)')/50;
H_bond_end=zeros(length(H_bond_list),length(time_calc)+1);

for i=1:length(H_bond_list)
    i
    t=H_bond_list(i,4);
    iter=0;
    for j=1:length(time_calc)
        iter=iter+1;
        a3 =cell2mat(C(t+time_calc(iter),:));
        o_a=[a3(H_bond_list(i,1),4) a3(H_bond_list(i,1),5) a3(H_bond_list(i,1),6)];
        o_d=[a3(H_bond_list(i,2),4) a3(H_bond_list(i,2),5) a3(H_bond_list(i,2),6)];
        hyd=[a3(H_bond_list(i,3),4) a3(H_bond_list(i,3),5) a3(H_bond_list(i,3),6)];
        
        if iter>1
        a4 = cell2mat(C(t+time_calc(iter-1),:));
        o_a_prev=[a4(H_bond_list(i,1),4) a4(H_bond_list(i,1),5) a4(H_bond_list(i,1),6)];
        o_d_prev=[a4(H_bond_list(i,2),4) a4(H_bond_list(i,2),5) a4(H_bond_list(i,2),6)];
        hyd_prev=[a4(H_bond_list(i,3),4) a4(H_bond_list(i,3),5) a4(H_bond_list(i,3),6)];
        elseif iter==1
        a4 = cell2mat(C(t,:));
        o_a_prev=[a4(H_bond_list(i,1),4) a4(H_bond_list(i,1),5) a4(H_bond_list(i,1),6)];
        o_d_prev=[a4(H_bond_list(i,2),4) a4(H_bond_list(i,2),5) a4(H_bond_list(i,2),6)];
        hyd_prev=[a4(H_bond_list(i,3),4) a4(H_bond_list(i,3),5) a4(H_bond_list(i,3),6)];
        end
        
        if (o_a(1)-o_a_prev(1))>(xhi-xlo)-f_per
            o_a_x_shift(i)=o_a_x_shift(i)-(xhi-xlo);
        elseif (o_a(1)-o_a_prev(1))<f_per-(xhi-xlo)
            o_a_x_shift(i)=o_a_x_shift(i)+(xhi-xlo);
        end
        
        if (o_a(2)-o_a_prev(2))>(yhi-ylo)-f_per
            o_a_y_shift(i)=o_a_y_shift(i)-(yhi-ylo);
        elseif (o_a(2)-o_a_prev(2))<f_per-(yhi-ylo)
            o_a_y_shift(i)=o_a_y_shift(i)+(yhi-ylo);
        end
        
        if (o_d(1)-o_d_prev(1))>(xhi-xlo)-f_per
            o_d_x_shift(i)=o_d_x_shift(i)-(xhi-xlo);
        elseif (o_d(1)-o_d_prev(1))<f_per-(xhi-xlo)
            o_d_x_shift(i)=o_d_x_shift(i)+(xhi-xlo);
        end
        
        if (o_d(2)-o_d_prev(2))>(yhi-ylo)-f_per
            o_d_y_shift(i)=o_d_y_shift(i)-(yhi-ylo);
        elseif (o_d(2)-o_d_prev(2))<f_per-(yhi-ylo)
            o_d_y_shift(i)=o_d_y_shift(i)+(yhi-ylo);
        end
        
        if (hyd(1)-hyd_prev(1))>(xhi-xlo)-f_per
            hyd_x_shift(i)=hyd_x_shift(i)-(xhi-xlo);
        elseif (hyd(1)-hyd_prev(1))<f_per-(xhi-xlo)
            hyd_x_shift(i)=hyd_x_shift(i)+(xhi-xlo);
        end
        
        if (hyd(2)-hyd_prev(2))>(yhi-ylo)-f_per
            hyd_y_shift(i)=hyd_y_shift(i)-(yhi-ylo);
        elseif (hyd(2)-hyd_prev(2))<f_per-(yhi-ylo)
            hyd_y_shift(i)=hyd_y_shift(i)+(yhi-ylo);
        end
        
        o_a(1)=o_a(1)+o_a_x_shift(i);
        o_a(2)=o_a(2)+o_a_y_shift(i);
        
        o_d(1)=o_d(1)+o_d_x_shift(i);
        o_d(2)=o_d(2)+o_d_y_shift(i);
        
        hyd(1)=hyd(1)+hyd_x_shift(i);
        hyd(2)=hyd(2)+hyd_y_shift(i);
        
        h_bond_exist=check_hbond(o_a,o_d,hyd);
        if h_bond_exist==0
            H_bond_end(i,j+1)=0;
        elseif h_bond_exist==1
            H_bond_end(i,j+1)=1;
        end
    end
end
H_bond_end(:,1)=1;

C_HB=zeros(1,length(time_calc)+1);
for k=1:length(time_calc)+1
    C_HB(k)=mean(H_bond_end(:,k));
end
time_ps=vertcat([0],time_calc)*0.05;

e=2.71828;
tau_C_ps=max(time_ps(find(C_HB>1/e)))
plot(time_ps,C_HB)
%plot(time_ps,C_HB,time_ps,exp(-time_ps/tau_C_ps))
save('hbond_long_timescale_bulk_periodic_GD01_300K.mat','C_HB','time_ps','tau_C_ps')

toc