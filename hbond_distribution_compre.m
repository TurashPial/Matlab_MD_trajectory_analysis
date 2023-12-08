
clear all
tic
fs2=5;
n = 100;
oo_crit=3.4; %h-bond defination 
oh_crit=2.425;
angle_crit=30;
fs3=oo_crit;
%% read data
fid = fopen('coords_nvt_4.dump','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
fclose(fid);
box = cell2mat(box1);
xlo = box(1,1);xhi = box(1,2);
ylo = box(2,1);yhi = box(2,2);
zlo = box(3,1);zhi = box(3,2);
%%
fid = fopen('coords_nvt_4.dump','r');
for ii= 1:n
    C(ii,:) = textscan(fid,'%f %f %f %f %f %f','Headerlines',9);
end

%C contain 1-id 2-x 3-y 4-z 5-vx 6-vy 7-vz
fclose(fid);

t1 = cell2mat(C(1,:));
a6 = find(t1(:,3)==6);
z_graft = min(t1(a6,6));
z_brush = max(t1(a6,6));
brush_avg_height_iter=linspace(0,0,n);
   
z_block=vertcat(linspace(z_graft+fs3,z_brush,6)',linspace(z_brush,zhi-40-fs3,11)');
z_block(6)=[];
H_bond=linspace(0,0,length(z_block)-1);
n_water=linspace(0,0,length(z_block)-1);

for t=1:n
   t
   a1 = cell2mat(C(t,:));
   a6 = find(a1(:,3)==6);
   brush_avg_height_iter(t)=(mean(a1(a6,6))-z_graft)*2;

   for q=1:length(z_block)-1    
   a22 = find(a1(:,3)==7 & a1(:,6)>z_block(q)-fs3 & a1(:,6)<z_block(q+1)+fs3 & a1(:,4)<xhi-fs2+fs3 & a1(:,4)>xlo+fs2-fs3 ...
       & a1(:,5)<yhi-fs2+fs3 & a1(:,5)>ylo+fs2-fs3);    
   a21 = find(a1(:,3)==7 & a1(:,6)>z_block(q) & a1(:,6)<z_block(q+1) & a1(:,4)<xhi-fs2 & a1(:,4)>xlo+fs2 ...
       & a1(:,5)<yhi-fs2 & a1(:,5)>ylo+fs2);
   a24 = find(a1(:,3)==5 & a1(:,6)>z_block(q)-fs3 & a1(:,6)<z_block(q+1)+fs3 & a1(:,4)<xhi-fs2+fs3 & a1(:,4)>xlo+fs2-fs3 ...
       & a1(:,5)<yhi-fs2+fs3 & a1(:,5)>ylo+fs2-fs3);  
   a26 = find(a1(:,3)==11 & a1(:,6)>z_block(q)-fs3 & a1(:,6)<z_block(q+1)+fs3 & a1(:,4)<xhi-fs2+fs3 & a1(:,4)>xlo+fs2-fs3 ...
       & a1(:,5)<yhi-fs2+fs3 & a1(:,5)>ylo+fs2-fs3);  
   n_water(q)=n_water(q)+length(a21);
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
       
       for m=1:length(a22)
           x_o_ne=a1(a22(m),4);
           y_o_ne=a1(a22(m),5);
           z_o_ne=a1(a22(m),6);
           oo_dis=sqrt((x_o-x_o_ne)^2+(y_o-y_o_ne)^2+(z_o-z_o_ne)^2);
                if oo_dis==0 || oo_dis>oo_crit
                     continue
                end
       x_h1_ne=a1(a22(m)+1,4);
       y_h1_ne=a1(a22(m)+1,5);
       z_h1_ne=a1(a22(m)+1,6);
       x_h2_ne=a1(a22(m)+2,4);
       y_h2_ne=a1(a22(m)+2,5);
       z_h2_ne=a1(a22(m)+2,6);
       
       o_h1ne_dis=sqrt((x_o-x_h1_ne)^2+(y_o-y_h1_ne)^2+(z_o-z_h1_ne)^2);
       u_o_one_h1ne=[(x_o-x_o_ne) (y_o-y_o_ne) (z_o-z_o_ne)];
       v_o_one_h1ne=[(x_h1_ne-x_o_ne) (y_h1_ne-y_o_ne) (z_h1_ne-z_o_ne)];
       angle_o_one_h1ne=acosd(dot(u_o_one_h1ne,v_o_one_h1ne)/(norm(u_o_one_h1ne)*norm(v_o_one_h1ne)));
  
       if o_h1ne_dis<oh_crit & angle_o_one_h1ne<angle_crit
           H_bond(q)=H_bond(q)+1;
           continue
       end
       
       o_h2ne_dis=sqrt((x_o-x_h2_ne)^2+(y_o-y_h2_ne)^2+(z_o-z_h2_ne)^2);
       u_o_one_h2ne=[(x_o-x_o_ne) (y_o-y_o_ne) (z_o-z_o_ne)];
       v_o_one_h2ne=[(x_h2_ne-x_o_ne) (y_h2_ne-y_o_ne) (z_h2_ne-z_o_ne)];
       angle_o_one_h2ne=acosd(dot(u_o_one_h2ne,v_o_one_h2ne)/(norm(u_o_one_h2ne)*norm(v_o_one_h2ne)));
       
       if o_h2ne_dis<oh_crit & angle_o_one_h2ne<angle_crit
           H_bond(q)=H_bond(q)+1;
           continue
       end
       
       h1_one_dis=sqrt((x_h1-x_o_ne)^2+(y_h1-y_o_ne)^2+(z_h1-z_o_ne)^2);
       u_h1_o_one=[(x_h1-x_o) (y_h1-y_o) (z_h1-z_o)];
       v_h1_o_one=[(x_o_ne-x_o) (y_o_ne-y_o) (z_o_ne-z_o)];
       angle_h1_o_one=acosd(dot(u_h1_o_one,v_h1_o_one)/(norm(u_h1_o_one)*norm(v_h1_o_one)));
       
       if h1_one_dis<oh_crit & angle_h1_o_one<angle_crit
           H_bond(q)=H_bond(q)+1;
           continue
       end
       
       h2_one_dis=sqrt((x_h2-x_o_ne)^2+(y_h2-y_o_ne)^2+(z_h2-z_o_ne)^2);
       u_h2_o_one=[(x_h2-x_o) (y_h2-y_o) (z_h2-z_o)];
       v_h2_o_one=[(x_o_ne-x_o) (y_o_ne-y_o) (z_o_ne-z_o)];
       angle_h2_o_one=acosd(dot(u_h2_o_one,v_h2_o_one)/(norm(u_h2_o_one)*norm(v_h2_o_one)));
       
       if h2_one_dis<oh_crit & angle_h2_o_one<angle_crit
           H_bond(q)=H_bond(q)+1;
           continue
       end
       end 
       
       for s=1:length(a24) 
           x_ocar_ne=a1(a24(s),4);
           y_ocar_ne=a1(a24(s),5);
           z_ocar_ne=a1(a24(s),6);
           oocar_dis=sqrt((x_o-x_ocar_ne)^2+(y_o-y_ocar_ne)^2+(z_o-z_ocar_ne)^2);
                if oocar_dis>oo_crit
                     continue
                end
           ocarne_h1_dis=sqrt((x_ocar_ne-x_h1)^2+(y_ocar_ne-y_h1)^2+(z_ocar_ne-z_h1)^2);
           u_ocarne_o_h1=[(x_ocar_ne-x_o) (y_ocar_ne-y_o) (z_ocar_ne-z_o)];
           v_ocarne_o_h1=[(x_h1-x_o) (y_h1-y_o) (z_h1-z_o)];
           angle_ocarne_o_h1=acosd(dot(u_ocarne_o_h1,v_ocarne_o_h1)/(norm(u_ocarne_o_h1)*norm(v_ocarne_o_h1)));
       
       if ocarne_h1_dis<oh_crit & angle_ocarne_o_h1<angle_crit
           H_bond(q)=H_bond(q)+1;
           continue
       end
       
           ocarne_h2_dis=sqrt((x_ocar_ne-x_h2)^2+(y_ocar_ne-y_h2)^2+(z_ocar_ne-z_h2)^2);
           u_ocarne_o_h2=[(x_ocar_ne-x_o) (y_ocar_ne-y_o) (z_ocar_ne-z_o)];
           v_ocarne_o_h2=[(x_h2-x_o) (y_h2-y_o) (z_h2-z_o)];
           angle_ocarne_o_h2=acosd(dot(u_ocarne_o_h2,v_ocarne_o_h2)/(norm(u_ocarne_o_h2)*norm(v_ocarne_o_h2)));
       
       if ocarne_h2_dis<oh_crit & angle_ocarne_o_h2<angle_crit
           H_bond(q)=H_bond(q)+1;
           continue
       end
       end
       
       for v=1:length(a26)
           x_cl_ne=a1(a26(v),4);
           y_cl_ne=a1(a26(v),5);
           z_cl_ne=a1(a26(v),6);
           ocl_dis=sqrt((x_o-x_cl_ne)^2+(y_o-y_cl_ne)^2+(z_o-z_cl_ne)^2);
                if ocl_dis>oo_crit
                     continue
                end
           clne_h1_dis=sqrt((x_cl_ne-x_h1)^2+(y_cl_ne-y_h1)^2+(z_cl_ne-z_h1)^2);
           u_clne_o_h1=[(x_cl_ne-x_o) (y_cl_ne-y_o) (z_cl_ne-z_o)];
           v_clne_o_h1=[(x_h1-x_o) (y_h1-y_o) (z_h1-z_o)];
           angle_clne_o_h1=acosd(dot(u_clne_o_h1,v_clne_o_h1)/(norm(u_clne_o_h1)*norm(v_clne_o_h1)));
       
       if clne_h1_dis<oh_crit & angle_clne_o_h1<angle_crit
           H_bond(q)=H_bond(q)+1;
           continue
       end
       
           clne_h2_dis=sqrt((x_cl_ne-x_h2)^2+(y_cl_ne-y_h2)^2+(z_cl_ne-z_h2)^2);
           u_clne_o_h2=[(x_cl_ne-x_o) (y_cl_ne-y_o) (z_cl_ne-z_o)];
           v_clne_o_h2=[(x_h2-x_o) (y_h2-y_o) (z_h2-z_o)];
           angle_clne_o_h2=acosd(dot(u_clne_o_h2,v_clne_o_h2)/(norm(u_clne_o_h2)*norm(v_clne_o_h2)));
       
       if clne_h2_dis<oh_crit & angle_clne_o_h2<angle_crit
           H_bond(q)=H_bond(q)+1;
           continue
       end
       end
       
       end
   end

end
%%
H_bond_per_water=H_bond./n_water;
for i=1:length(z_block)-1
    z_cen(i)=(z_block(i)+z_block(i+1))/2;
end
brush_avg_height=mean(brush_avg_height_iter);
plot(z_cen-z_graft,H_bond_per_water,brush_avg_height,mean(H_bond_per_water),'o')
toc