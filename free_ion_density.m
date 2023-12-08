%------------non bound ion analysis -------
clear all
tic
cutoff=6.2;
nsteps=round((cutoff/0.08)+1);
n = 300;
%% read simulation data
fid = fopen('prod_2_c12_ion_S4_3.data','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
fclose(fid);
box = cell2mat(box1);
xlo = box(1,1);xhi = box(1,2);
ylo = box(2,1);yhi = box(2,2);
zlo = box(3,1);zhi = box(3,2);

delta_x=xhi-xlo;
delta_y=yhi-ylo;
fid = fopen('prod_2_c12_ion_S4_3_2.data','r');
for ii= 1:n
    C(ii,:) = textscan(fid,'%f %f %f %f %f','Headerlines',9);
end
%C contain 1-id 2-x 3-y 4-z 5-vx 6-vy 7-vz
fclose(fid);

%%
an=0;
aaa=0;
Xall=[];

nb = linspace(xlo,xhi,2)'; % use number of bins so that per chuck is around ~~ nm
nbins = length(nb)-1;  %how many bins in required direction to get the profile
ny=linspace(ylo,yhi,200)';
nybin = length(ny)-1;
nz=linspace(zlo,zhi,200)';
nzbin = length(nz)-1;
%%
for t=1:n
   t
   a1 = cell2mat(C(t,:));
   a22 = find(a1(:,2)==7 & a1(:,4)<360 & a1(:,4)>120 & a1(:,5)<250);
   a222 = find(a1(:,2)==8 & a1(:,4)<360 & a1(:,4)>120 & a1(:,5)<250);
   a231 = find(a1(:,2)==5 & a1(:,4)<360 & a1(:,4)>120 & a1(:,5)<250);
   a232 = find(a1(:,2)==6 & a1(:,4)<360 & a1(:,4)>120 & a1(:,5)<250);
   a233 = find(a1(:,2)==4 & a1(:,4)<360 & a1(:,4)>120 & a1(:,5)<250);
   a23=vertcat(a231,a232,a233);
   %a22=vertcat(a221,a222);
   distance_min=linspace(0,0,length(a22));

   %find free Na+
   for ii=1:length(a22)
       k=a22(ii);
       x_ion=a1(k,3);
       y_ion=a1(k,4);
       z_ion=a1(k,5);
     
       for jj=1:length(a23)
           w=a23(jj);
           x_water=a1(w,3);
           y_water=a1(w,4);
           z_water=a1(w,5);
           distance(jj)=sqrt((x_ion-x_water)^2+(y_ion-y_water)^2+(z_ion-z_water)^2);
       end
       distance_min(1,ii)=find(distance == min(distance));
       distance_min(2,ii) = distance(distance_min(1,ii));
   end
   a22_free_dist = find(distance_min(2,:)>cutoff);
   a22_free = (a22(a22_free_dist));

% plot distribution
Xtime=[0;0];
   for kk = 1:nbins
       clear loc
       a211 = find(a1(:,2)==1 & a1(:,4)<360 & a1(:,4)>120 & a1(:,5)<250);  %as we are using replicate Y 2,      
       a221 = find(a1(:,2)==2 & a1(:,4)<360 & a1(:,4)>120 & a1(:,5)<250);
       a231 = find(a1(:,2)==3 & a1(:,4)<360 & a1(:,4)>120 & a1(:,5)<250);
       a241 = find(a1(:,2)==4 & a1(:,4)<360 & a1(:,4)>120 & a1(:,5)<250);
       a5=a22_free;
       a2=vertcat(a211, a221, a231, a241);

       y_dir=a1(a2,4);
       z_dir=a1(a2,5);
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
        loc(1,ij)=a1(a5(ij),4)-cent_mean(1);
        loc(2,ij)=a1(a5(ij),5)-cent_mean(2);
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
%% contour plot
data_point=101;
data_point_n=(data_point+1)/2;
pts = linspace(-80, 80, data_point);

N=histcounts2(Xall(2,:),Xall(1,:),pts,pts);
N(data_point_n,data_point_n)=N(data_point_n,data_point_n)-n;
imagesc(pts, pts, N);
colorbar