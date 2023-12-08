clear all
tic
%% cutoff distances
n = 5;
FHS1=3;
FHS2=2.4;
%% initialize data
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
n_water=0;
n_bound_OC=0;
n_bound_HO1=0;
n_bound_HO2=0;

%% find bound water 
for t=1:n
   t
   a1 = cell2mat(C(t,:));
   a6 = find(a1(:,3)==6);
   z_graft(t,1) = min(a1(a6,6));
   z_ends=a1(a6,6);
   z_brush(t,1) = min(z_ends(z_ends>z_graft(t,1)));
   a27 = find(a1(:,3)==7 & a1(:,6)>-82 & a1(:,6)<-61 & a1(:,5)>3 & a1(:,5)<73 & a1(:,4)>3 & a1(:,4)<73); %O_water
   %a24 = find(a1(:,3)==9 & a1(:,6)>-84-FHS1 & a1(:,6)<-62+FHS1 & a1(:,5)>3-FHS1 & a1(:,5)<73+FHS1 & a1(:,4)>3-FHS1 & a1(:,4)<73+FHS1); #counterion
   %a28 = find(a1(:,3)==8 & a1(:,6)>-85 & a1(:,6)<-61 & a1(:,5)>2 & a1(:,5)<74 & a1(:,4)>2 & a1(:,4)<74); #H_water
   %a25 = find(a1(:,3)==5 & a1(:,6)>-85-FHS2 & a1(:,6)<-61+FHS2 & a1(:,5)>2-FHS2 & a1(:,5)<74+FHS2 & a1(:,4)>2-FHS2 & a1(:,4)<74+FHS2); #O_Coo
   a281=a27+1;
   a282=a27+2;
   mo1=n_bound_OC;
   mo2=n_bound_HO1;
   mo3=n_bound_HO2;
   for i=1:length(a27)
       %i
   m=a27(i);
   n1=a281(i);
   n2=a282(i);
   n_water=n_water+1;
   x0=a1(m,4);
   y0=a1(m,5);
   z0=a1(m,6);
   x1=a1(n1,4);
   y1=a1(n1,5);
   z1=a1(n1,6);
   x2=a1(n2,4);
   y2=a1(n2,5);
   z2=a1(n2,6);
   a79 = find(a1(:,3)==9  & (a1(:,6)-z0).^2+(a1(:,5)-y0).^2+(a1(:,4)-x0).^2<FHS1^2);
   a815 = find(a1(:,3)==5  & (a1(:,6)-z1).^2+(a1(:,5)-y1).^2+(a1(:,4)-x1).^2<FHS2^2);
   a825 = find(a1(:,3)==5  & (a1(:,6)-z2).^2+(a1(:,5)-y2).^2+(a1(:,4)-x2).^2<FHS2^2);
   
   if length(a79)>0
       n_bound_OC=n_bound_OC+1;
       mol1(n_bound_OC-mo1)=a1(a27(i),2);
   end
   if length(a815)>0
       n_bound_HO1=n_bound_HO1+1;
       mol2(n_bound_HO1-mo2)=a1(a27(i),2);
   end
      if length(a825)>0
       n_bound_HO2=n_bound_HO2+1;
       mol3(n_bound_HO2-mo3)=a1(a27(i),2);
   end
   end
mol=vertcat(mol1',mol2',mol3');
bound_count=unique(mol);
bound(t)=length(bound_count);
end
%% normalize
wat_bound_count=sum(bound);
percent_bound=(wat_bound_count)/n_water*100
percent_free=100-percent_bound
toc