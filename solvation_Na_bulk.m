%-------------solvation behavior of cation -----------
clear all
tic
n = 200;
FHS=3.2; %cutoff
%% read simulation data
fid = fopen('coords_all_1V_6th.dump','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
fclose(fid);
box = cell2mat(box1);
xlo = box(1,1);xhi = box(1,2);
ylo = box(2,1);yhi = box(2,2);
zlo = box(3,1);zhi = box(3,2);

delta_x=xhi-xlo;
delta_y=yhi-ylo;
fid = fopen('coords_all_1V_6th.dump','r');
for ii= 1:n
    C(ii,:) = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f','Headerlines',9);
end
%C contain 1-id 2-x 3-y 4-z 5-vx 6-vy 7-vz

fclose(fid);
%%
nion=0;
wat_count=0;
car_count=0;
for t=1:n
   t
   a1 = cell2mat(C(t,:));
   a21 = find(a1(:,2)==6);
   alo = find(a1(:,2)==6 & a1(:,5)<110); %lower brush
   z_brush_lo(t,1) = max(a1(alo,5));
   ahi = find(a1(:,2)==6 & a1(:,5)>110); %upper brush
   z_brush_hi(t,1) = min(a1(ahi,5));
   a221 = find(a1(:,2)==9 & a1(:,5)>z_brush_lo(t,1)+8 & a1(:,5)<z_brush_hi(t,1)-8 & a1(:,3)>50);
   a222 = find(a1(:,2)==10 & a1(:,5)>z_brush_lo(t,1)+8 & a1(:,5)<z_brush_hi(t,1)-8 & a1(:,3)>50);
   a22 = vertcat(a221,a222); %ion
   a23 = find(a1(:,2)==7 & a1(:,5)>z_brush_lo(t,1) & a1(:,5)<z_brush_hi(t,1)-5 & a1(:,3)>45); %water
   a24 = find(a1(:,2)==5 ); %coo
   for ii=1:length(a22)
       nion=nion+1;
       k=a22(ii);
       x_ion=a1(k,3);
       y_ion=a1(k,4);
       z_ion=a1(k,5);
       for jj=1:length(a23)
           w=a23(jj);
           x_water=a1(w,3);
           y_water=a1(w,4);
           z_water=a1(w,5);
           distance_water=sqrt((x_ion-x_water)^2+(y_ion-y_water)^2+(z_ion-z_water)^2);
           distance_xhi=sqrt((x_ion-(x_water+delta_x))^2+(y_ion-y_water)^2+(z_ion-z_water)^2);
           distance_xlo=sqrt((x_ion-(x_water-delta_x))^2+(y_ion-y_water)^2+(z_ion-z_water)^2);
           distance_yhi=sqrt((x_ion-x_water)^2+(y_ion-(y_water+delta_y))^2+(z_ion-z_water)^2);
           distance_ylo=sqrt((x_ion-x_water)^2+(y_ion-(y_water-delta_y))^2+(z_ion-z_water)^2);
           distance_xhi_yhi=sqrt((x_ion-(x_water+delta_x))^2+(y_ion-(y_water+delta_y))^2+(z_ion-z_water)^2);
           distance_xhi_ylo=sqrt((x_ion-(x_water+delta_x))^2+(y_ion-(y_water-delta_y))^2+(z_ion-z_water)^2);
           distance_xlo_yhi=sqrt((x_ion-(x_water-delta_x))^2+(y_ion-(y_water+delta_y))^2+(z_ion-z_water)^2);
           distance_xlo_ylo=sqrt((x_ion-(x_water-delta_x))^2+(y_ion-(y_water-delta_y))^2+(z_ion-z_water)^2);
           if distance_water<FHS
               wat_count=wat_count+1;
           end
           if distance_xhi<FHS
               wat_count=wat_count+1;
           end
           if distance_xlo<FHS
               wat_count=wat_count+1;
           end
           if distance_yhi<FHS
               wat_count=wat_count+1;
           end
           if distance_ylo<FHS
               wat_count=wat_count+1;
           end
           if distance_xhi_yhi<FHS
               wat_count=wat_count+1;
           end
           if distance_xhi_ylo<FHS
               wat_count=wat_count+1;
           end
           if distance_xlo_yhi<FHS
               wat_count=wat_count+1;
           end
           if distance_xlo_ylo<FHS
               wat_count=wat_count+1;
           end
       end
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
           end
           if distance_car_xhi<FHS
               car_count=car_count+1;
           end
           if distance_car_xlo<FHS
               car_count=car_count+1;
           end
           if distance_car_yhi<FHS
               car_count=car_count+1;
           end
           if distance_car_ylo<FHS
               car_count=car_count+1;
           end
           if distance_car_xhi_yhi<FHS
               car_count=car_count+1;
           end
           if distance_car_xhi_ylo<FHS
               car_count=car_count+1;
           end
           if distance_car_xlo_yhi<FHS
               car_count=car_count+1;
           end
           if distance_car_xlo_ylo<FHS
               car_count=car_count+1;
           end
       end
   end
end
%% normalize

wat_count_av = wat_count/nion;
car_count_av = car_count/nion;
total_count_av=wat_count_av + car_count_av; 

n_counterions=0;
n_water=0;
for t=1:n
     t
     a1 = cell2mat(C(t,:));
     a261 = length(find(a1(:,2)==9 & a1(:,5)>z_brush_lo(t,1)+5 & a1(:,5)<z_brush_hi(t,1)-5 & a1(:,3)>70));
     a262= length(find(a1(:,2)==10 & a1(:,5)>z_brush_lo(t,1)+5 & a1(:,5)<z_brush_hi(t,1)-5 & a1(:,3)>70));
     a25=a261+a262;
     a26 = length(find(a1(:,2)==7 & a1(:,5)>z_brush_lo(t,1) & a1(:,5)<z_brush_hi(t,1)-5 & a1(:,3)>70));
     n_counterions=n_counterions+a25;
     n_water=n_water+a26;
end

NA=6.022*10^(23);
molarity_av=(n_counterions/NA)/(n)/(xhi-70)/(yhi-ylo)/(z_brush_hi(t,1)-5-z_brush_lo(t,1)-5)*10^(27);
molality_av=n_counterions/(n_water/55.55);
toc