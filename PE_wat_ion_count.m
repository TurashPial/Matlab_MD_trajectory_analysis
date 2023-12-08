%_______________this code count the number of water and PE in the solvation shell of an ion. ___________
clear all
tic
n = 100; %number of snapshot
FHS=3.2; %cutoff
%% read simulation data
fid = fopen('coords_nvt_n9.dump','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
fclose(fid);
box = cell2mat(box1);
xlo = box(1,1);xhi = box(1,2);
ylo = box(2,1);yhi = box(2,2);
zlo = box(3,1);zhi = box(3,2);

delta_x=xhi-xlo;
delta_y=yhi-ylo;
fid = fopen('coords_nvt_n9.dump','r');
for ii= 1:n
    D(ii,:) = textscan(fid,'%f %f %f %f %f %f','Headerlines',9);
end
fclose(fid);

clear C
for jj=001:100
    C(jj-000,:)=D(jj,:);
end
%C contain 1-id 2-x 3-y 4-z 5-vx 6-vy 7-vz
ncar=0;
wat_count=0;
ion_count=0;
%% residence shell
for t=1:100
   t
   a1 = cell2mat(C(t,:));
   a21 = find(a1(:,3)==6);
   z_graft(t,1) = min(a1(a21,6));
   z_ends=a1(a21,6);
   z_brush(t,1) = min(z_ends(z_ends>z_graft(t,1)));
   a22 = find(a1(:,3)==9 & a1(:,6)<z_brush(t,1)-FHS & a1(:,6)>z_graft(t,1)+FHS & a1(:,4)<xhi...
       & a1(:,4)>xlo & a1(:,5)<yhi & a1(:,5)>ylo);
   a23 = find(a1(:,3)==7 & a1(:,6)<z_brush(t,1) & a1(:,6)>z_graft(t,1));
   a24 = find(a1(:,3)==5 & a1(:,6)<z_brush(t,1) & a1(:,6)>z_graft(t,1));
   for ii=1:length(a22)
       ncar=ncar+1;
       k=a22(ii);
       x_car=a1(k,4);
       y_car=a1(k,5);
       z_car=a1(k,6);
       for jj=1:length(a23)
           w=a23(jj);
           x_water=a1(w,4);
           y_water=a1(w,5);
           z_water=a1(w,6);
           distance_water=sqrt((x_car-x_water)^2+(y_car-y_water)^2+(z_car-z_water)^2);
           distance_xhi=sqrt((x_car-(x_water+delta_x))^2+(y_car-y_water)^2+(z_car-z_water)^2);
           distance_xlo=sqrt((x_car-(x_water-delta_x))^2+(y_car-y_water)^2+(z_car-z_water)^2);
           distance_yhi=sqrt((x_car-x_water)^2+(y_car-(y_water+delta_y))^2+(z_car-z_water)^2);
           distance_ylo=sqrt((x_car-x_water)^2+(y_car-(y_water-delta_y))^2+(z_car-z_water)^2);
           distance_xhi_yhi=sqrt((x_car-(x_water+delta_x))^2+(y_car-(y_water+delta_y))^2+(z_car-z_water)^2);
           distance_xhi_ylo=sqrt((x_car-(x_water+delta_x))^2+(y_car-(y_water-delta_y))^2+(z_car-z_water)^2);
           distance_xlo_yhi=sqrt((x_car-(x_water-delta_x))^2+(y_car-(y_water+delta_y))^2+(z_car-z_water)^2);
           distance_xlo_ylo=sqrt((x_car-(x_water-delta_x))^2+(y_car-(y_water-delta_y))^2+(z_car-z_water)^2);
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
           x_ion=a1(p,4);
           y_ion=a1(p,5);
           z_ion=a1(p,6);
           distance_ion=sqrt((x_car-x_ion)^2+(y_car-y_ion)^2+(z_car-z_ion)^2);
           distance_ion_xhi=sqrt((x_car-(x_ion+delta_x))^2+(y_car-y_ion)^2+(z_car-z_ion)^2);
           distance_ion_xlo=sqrt((x_car-(x_ion-delta_x))^2+(y_car-y_ion)^2+(z_car-z_ion)^2);
           distance_ion_yhi=sqrt((x_car-x_ion)^2+(y_car-(y_ion+delta_y))^2+(z_car-z_ion)^2);
           distance_ion_ylo=sqrt((x_car-x_ion)^2+(y_car-(y_ion-delta_y))^2+(z_car-z_ion)^2);
           distance_ion_xhi_yhi=sqrt((x_car-(x_ion+delta_x))^2+(y_car-(y_ion+delta_y))^2+(z_car-z_ion)^2);
           distance_ion_xhi_ylo=sqrt((x_car-(x_ion+delta_x))^2+(y_car-(y_ion-delta_y))^2+(z_car-z_ion)^2);
           distance_ion_xlo_yhi=sqrt((x_car-(x_ion-delta_x))^2+(y_car-(y_ion+delta_y))^2+(z_car-z_ion)^2);
           distance_ion_xlo_ylo=sqrt((x_car-(x_ion-delta_x))^2+(y_car-(y_ion-delta_y))^2+(z_car-z_ion)^2);
           if distance_ion<FHS
               ion_count=ion_count+1;
           end
           if distance_ion_xhi<FHS
               ion_count=ion_count+1;
           end
           if distance_ion_xlo<FHS
               ion_count=ion_count+1;
           end
           if distance_ion_yhi<FHS
               ion_count=ion_count+1;
           end
           if distance_ion_ylo<FHS
               ion_count=ion_count+1;
           end
           if distance_ion_xhi_yhi<FHS
               ion_count=ion_count+1;
           end
           if distance_ion_xhi_ylo<FHS
               ion_count=ion_count+1;
           end
           if distance_ion_xlo_yhi<FHS
               ion_count=ion_count+1;
           end
           if distance_ion_xlo_ylo<FHS
               ion_count=ion_count+1;
           end
       end
   end
end
%% normalize
wat_count_av = wat_count/ncar;
ion_count_av = ion_count/ncar; 

n_counterions=0;
n_water=0;
for t=1:100
     t
     a25 = length(find(a1(:,3)==9 & a1(:,6)<z_brush(1,1) & a1(:,6)>z_graft(1,1)));
     a26 = length(find(a1(:,3)==7 & a1(:,6)<z_brush(1,1) & a1(:,6)>z_graft(1,1)));
     n_counterions=n_counterions+a25;
     n_water=n_water+a26;
end

NA=6.022*10^(23);
molarity_av=(n_counterions/NA)/n/(xhi-xlo)/(yhi-ylo)/(z_brush(1,1)-z_graft(1,1))*10^(27);
molality_av=n_counterions/(n_water/55.55);
toc
