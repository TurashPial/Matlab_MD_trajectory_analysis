%----------------find atoms who form trimer interactions -------------
clear all
tic
cutoff=6;
nsteps=round((cutoff/0.08)+1);
n = 600
%% read simulation data 
fid = fopen('pmf_9_oe_na_cl12.dump','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
fclose(fid);
box = cell2mat(box1);
xlo = box(1,1);xhi = box(1,2);
ylo = box(2,1);yhi = box(2,2);
zlo = box(3,1);zhi = box(3,2);
%% read simulation data 
delta_x=xhi-xlo;
delta_y=yhi-ylo;
fid = fopen('pmf_9_oe_na_cl12.dump','r');
for ii= 1:n
    C(ii,:) = textscan(fid,'%f %f %f %f %f','Headerlines',9);
end
%C contain 1-id 2-x 3-y 4-z 5-vx 6-vy 7-vz
fclose(fid);

%%
an=0;
aaa=0;
Xall=[];
nb = linspace(zlo,zhi,2)'; % use number of bins so that per chuck is around ~~ nm
nbins = length(nb)-1;  %how many bins in required direction to get the profile

for t=1:n
   t
   a1 = cell2mat(C(t,:));
   a221 = find(a1(:,2)==37 & a1(:,3)<299 & a1(:,4)>130 & a1(:,5)<80);
   a222 = find(a1(:,2)==38 & a1(:,3)<299 & a1(:,4)>130 & a1(:,5)<80);
   a231 = find(a1(:,2)==39 & a1(:,3)<299 & a1(:,4)>130 & a1(:,5)<80);
   a232 = find(a1(:,2)==40 & a1(:,3)<299 & a1(:,4)>130 & a1(:,5)<80);
   %a233 = find(a1(:,2)==4 & a1(:,4)<360 & a1(:,4)>120 & a1(:,5)<250);
  % a23=vertcat(a231,a232,a233);
   a22=vertcat(a221,a222);  
   a23=vertcat(a231,a232);
   distance_min=linspace(0,0,length(a232));
   clear distance_min1

   %find Bjerrum pair 
   for ii=1:length(a232)
       k=a232(ii);
       x_ion=a1(k,3);
       y_ion=a1(k,4);
       z_ion=a1(k,5);
     clear distance
       for jj=1:length(a231)
           w=a231(jj);
           x_water=a1(w,3);
           y_water=a1(w,4);
           z_water=a1(w,5);
           distance(jj)=sqrt((x_ion-x_water)^2+(y_ion-y_water)^2+(z_ion-z_water)^2);
       end
       distance_min(1,ii)=find(distance == min(distance));
       distance_min(2,ii) = distance(distance_min(1,ii));
   end
   a_bound_dist = find(distance_min(2,:)<cutoff);
   a_bjurm = (a232(a_bound_dist));
   a_32 = (a232(a_bound_dist));
   a_31 = (a231(distance_min(1,a_bound_dist)));
   a_bj=horzcat(a_31,a_32);
%%%%%%%%%%%
   for iii=1:length(a_bj)
       a_rest=setdiff(a23,a_bj(iii,1:2)'); %get all other ion except this particular pair
       
       
       k_na=a_bj(iii,1); %first do for NA
       x_ion=a1(k_na,3);
       y_ion=a1(k_na,4);
       z_ion=a1(k_na,5);
     clear distance1
       for jjj=1:length(a_rest)
           w_rest=a_rest(jjj);
           x_water=a1(w_rest,3);
           y_water=a1(w_rest,4);
           z_water=a1(w_rest,5);
           distance1(jjj)=sqrt((x_ion-x_water)^2+(y_ion-y_water)^2+(z_ion-z_water)^2);
       end
       if min(distance1)<cutoff
       distance_min1(1,iii)=find(distance1 == min(distance1));
       distance_min1(2,iii) = distance1(distance_min1(1,iii));
       
       
       else
           k_na=a_bj(iii,2); % do for CL
       x_ion=a1(k_na,3);
       y_ion=a1(k_na,4);
       z_ion=a1(k_na,5);
     clear distance1
       for jjj=1:length(a_rest)
           w_rest=a_rest(jjj);
           x_water=a1(w_rest,3);
           y_water=a1(w_rest,4);
           z_water=a1(w_rest,5);
           distance2(jjj)=sqrt((x_ion-x_water)^2+(y_ion-y_water)^2+(z_ion-z_water)^2);
       end
       if min(distance2)<cutoff
       distance_min1(1,iii)=find(distance2 == min(distance2));
       distance_min1(2,iii) = distance2(distance_min1(1,iii));
       end
       end
   
   end
   %dist_min2(:,1:2)=nonzeros(distance_min1);
   a22_free_dist1 = find(0<distance_min1(2,:));
   a22_tri = (a_bj(a22_free_dist1));


% plot distribution
Xtime=[0;0];
   for kk = 1:nbins
       clear loc
   a2c1 = find(a1(:,2)==1 & a1(:,3)<280 & a1(:,4)>150 & a1(:,5)<80);
   a2c2 = find(a1(:,2)==2 & a1(:,3)<280 & a1(:,4)>150 & a1(:,5)<80);
   a2c3 = find(a1(:,2)==3 & a1(:,3)<280 & a1(:,4)>150 & a1(:,5)<80);
   a2c4 = find(a1(:,2)==4 & a1(:,3)<280 & a1(:,4)>150 & a1(:,5)<80);
       a5=a22_tri;
       a2=vertcat(a2c1, a2c2, a2c3, a2c4);

       y_dir=a1(a2,3);
       z_dir=a1(a2,4);
       X=[y_dir,z_dir];
       rng(5);
       GMModel = fitgmdist(X,2,'Replicates',10);
       aaa=aaa+1;
       cent(:,:)=GMModel.mu; % get center 
       cent=sortrows(cent,2);
         d = pdist(cent,'euclidean');
       cent_mean(1)=mean(cent(:,1));
       cent_mean(2)=mean(cent(:,2));

       for ij=1:length(a5)
        loc(1,ij)=a1(a5(ij),3)-cent_mean(1);
        loc(2,ij)=a1(a5(ij),4)-cent_mean(2);
        poly=polyfit(cent(:,1),cent(:,2),1);
        theta0=atan(poly(1));
        theta=-theta0;
%computation of rotation matrix. rotation is about X axis 
        rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        Xrot=rot*loc;   
        Xmean=mean(Xrot,2);
        Xcent=Xrot;
       end
Xtime=horzcat(Xtime,Xcent);
   end
Xall=horzcat(Xall,Xtime);
        
end
%% contour
data_point=81;
data_point_n=(data_point+1)/2;
pts = linspace(-100, 100, data_point);

N=histcounts2(Xall(2,:),Xall(1,:),pts,pts);
N(data_point_n,data_point_n)=(N(data_point_n,data_point_n)-n);
imagesc(pts, pts, N/n/(pts(2)-pts(1))/(pts(2)-pts(1))/(zhi-zlo)*2);
colorbar