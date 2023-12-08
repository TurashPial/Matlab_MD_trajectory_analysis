%--------------check bridging interactions ---
clear all
tic
n = 200;
FHS=5.3;
%% read data
fid = fopen('pmf_0_oe_cl_cal.data','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
fclose(fid);
box = cell2mat(box1);
xlo = box(1,1);xhi = box(1,2);
ylo = box(2,1);yhi = box(2,2);
zlo = box(3,1);zhi = box(3,2);

delta_x=xhi-xlo;
delta_y=yhi-ylo;
fid = fopen('pmf_13_oe_cl_cal.data','r');
for ii= 1:n
    C(ii,:) = textscan(fid,'%f %f %f %f %f','Headerlines',9);
end
%C contain 1-id 2-type 3-x 4-y 5-z

%%
nion=0;
wat_count=0;
car_count=0;
same_chain=0;
bridge_chain=0;

for t=1:n
   t;
   clear a22 a24 mol
   a1 = cell2mat(C(t,:));
   a22 = find(a1(:,2)==3);
   a241 = find(a1(:,2)==1);
   a242 = find(a1(:,2)==2);
   a24=vertcat(a241,a242);
   for ii=1:length(a22)
       nion=nion+1;
       k=a22(ii);
       x_ion=a1(k,3);
       y_ion=a1(k,4);
       z_ion=a1(k,5);
       nh=0;
       for kk=1:length(a24)
           p=a24(kk);
           x_car=a1(p,3);
           y_car=a1(p,4);
           z_car=a1(p,5);
           
           distance_car=sqrt((x_ion-x_car)^2+(y_ion-y_car)^2+(z_ion-z_car)^2);
           distance_car_xhi=sqrt((x_ion-(x_car+delta_x))^2+(y_ion-y_car)^2+(z_ion-z_car)^2);
           distance_car_xlo=sqrt((x_ion-(x_car-delta_x))^2+(y_ion-y_car)^2+(z_ion-z_car)^2);
           distance_car_yhi=sqrt((x_ion-x_car)^2+(y_ion-(y_car+delta_y))^2+(z_ion-z_car)^2);
           distance_car_ylo=sqrt((x_ion-x_car)^2+(y_ion-(y_car-delta_y))^2+(z_ion-z_car)^2);
           distance_car_xhi_yhi=sqrt((x_ion-(x_car+delta_x))^2+(y_ion-(y_car+delta_y))^2+(z_ion-z_car)^2);
           distance_car_xhi_ylo=sqrt((x_ion-(x_car+delta_x))^2+(y_ion-(y_car-delta_y))^2+(z_ion-z_car)^2);
           distance_car_xlo_yhi=sqrt((x_ion-(x_car-delta_x))^2+(y_ion-(y_car+delta_y))^2+(z_ion-z_car)^2);
           distance_car_xlo_ylo=sqrt((x_ion-(x_car-delta_x))^2+(y_ion-(y_car-delta_y))^2+(z_ion-z_car)^2);
           if distance_car<FHS
               car_count=car_count+1;
               nh=nh+1;
               mol(ii,nh)=a1(p,1);
           end
           if distance_car_xhi<FHS
               car_count=car_count+1;
               nh=nh+1;
               mol(ii,nh)=a1(p,1);
           end
           if distance_car_xlo<FHS
               car_count=car_count+1;
               nh=nh+1;
               mol(ii,nh)=a1(p,1);
           end
           if distance_car_yhi<FHS
               car_count=car_count+1;
               nh=nh+1;
               mol(ii,nh)=a1(p,1);
           end
           if distance_car_ylo<FHS
               car_count=car_count+1;
               nh=nh+1;
               mol(ii,nh)=a1(p,1);
           end
           if distance_car_xhi_yhi<FHS
               car_count=car_count+1;
               nh=nh+1;
               mol(ii,nh)=a1(p,1);
           end
           if distance_car_xhi_ylo<FHS
               car_count=car_count+1;
               nh=nh+1;
               mol(ii,nh)=a1(p,1);
           end
           if distance_car_xlo_yhi<FHS
               car_count=car_count+1;
               nh=nh+1;
               mol(ii,nh)=a1(p,1);
           end
           if distance_car_xlo_ylo<FHS
               car_count=car_count+1;
               nh=nh+1;
               mol(ii,nh)=a1(p,1);
               
           end
       end
   end

       sz=size(mol); 
       for ix=1:sz(1)
           for iy=1:sz(2)
               if mol(ix,iy)>600
               mol_mod(ix,iy)=2;
                              elseif mol(ix,iy)==0
                   mol_mod(ix,iy)=0;
               elseif 0<mol(ix,iy)<601
               mol_mod(ix,iy)=1;

               end
           end
       end
    for j=1:sz(1) %% check intra vs inter chain
       nonz=find(mol_mod(j,:)>0);
       nonzz=max(nonz);
       sd=std(mol_mod(j,1:nonzz));
       if sd==0
           same_chain=same_chain+1;
       else if sd>0
           bridge_chain=bridge_chain+1;
           
           end
       end
   end
end 
%%
bridge_avg=bridge_chain/(same_chain+bridge_chain);
bridge_chain/n

toc