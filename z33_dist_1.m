clear all
%------snapshot numbers
n=300
n1=200
n2=400
%% data initialization 
fid = fopen('prod_20_lmp.data','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
fclose(fid);
box = cell2mat(box1);
xlo = box(1,1);xhi = box(1,2);
ylo = box(2,1);yhi = box(2,2);
zlo = box(3,1);zhi = box(3,2);
%% read data files for various simulation stages
fid = fopen('prod_20_lmp_mdvwhole.data','r');
for ii= 1:n
    E(ii,:) = textscan(fid,'%f %f %f %f %f ','Headerlines',9);
end
fclose(fid);
fid = fopen('prod_21_lmp_mdvwhole.data','r');
for ii= n+1:n+n1
    E(ii,:) = textscan(fid,'%f %f %f %f %f ','Headerlines',9);
end
fclose(fid);
fid = fopen('prod_22_lmp.data','r');
for ii= n+n1+1:n+n1+n2
    E(ii,:) = textscan(fid,'%f %f %f %f %f ','Headerlines',9);
end
fclose(fid);
fid = fopen('prod_20_lmp_z33.data','r');
for ii= 1:n
    Ez(ii,:) = textscan(fid,'%f %f %f %f %f ','Headerlines',9);
end
fclose(fid);
fid = fopen('prod_21_lmp_z33.data','r');
for ii= n+1:n+n1
    Ez(ii,:) = textscan(fid,'%f %f %f %f %f ','Headerlines',9);
end
fclose(fid);
fid = fopen('prod_22_lmp_z33.data','r');
for ii= n+n1+1:n+n1+n2
    Ez(ii,:) = textscan(fid,'%f %f %f %f %f ','Headerlines',9);
end
fclose(fid);
%%
an=0;
aaa=0;
Xall=[];

nb = linspace(zlo,zhi,8)'; % use number of bins so that per chuck is around ~~ nm
nbins = length(nb)-1;  %how many bins in required direction to get the profile
ny=linspace(ylo,yhi,200)'; %as we replicated in x and y 
nybin = length(ny)-1;
nx=linspace(xlo,xhi,200)';
nxbin = length(nx)-1;
dist=0;
for jj = 1:n+n1+n2
    
    jj
   a1 = cell2mat(E(jj,1:5));
   a1z = cell2mat(Ez(jj,1:5));

   %only z33, not VVEE
   for iz=1:80
    for jz=1:37
        if jz>4
            az((iz-1)*33+jz-4,:)=a1z((iz-1)*37+jz,:);
        end
    end
   end
   %
   xlo=min(a1(:,3));
   xhi=max(a1(:,3));
   ylo=min(a1(:,4));
   yhi=max(a1(:,4));
   zlo=min(a1(:,5));
   zhi=max(a1(:,5));

Xtime=[0;0];
% for iz=1:length(az)
%            a1(az(iz),5)=a1(az(iz),5)-173.79;
%        end
   for kk = 1:nbins
       clear loc
       dist=dist+1;
       a211 = find(a1(:,2)==1 & nb(kk,1)<a1(:,5) & a1(:,5)<nb(kk+1,1) );   %C12
       a221 = find(a1(:,2)==2 & nb(kk,1)<a1(:,5) & a1(:,5)<nb(kk+1,1) );    %C12
       %a231 = find(a1(:,2)==3 & nb(kk,1)<a1(:,3) & a1(:,3)<nb(kk+1,1) & a1(:,4)>120 & a1(:,4)<340);
       %a232 = find(a1(:,2)==3 & nb(kk,1)<a1(:,3) & a1(:,3)<nb(kk+1,1) & a1(:,4)<208);
       %a241 = find(a1(:,2)==4 & nb(kk,1)<a1(:,3) & a1(:,3)<nb(kk+1,1) & a1(:,4)>120 & a1(:,4)<340);
       %a242 = find(a1(:,2)==4 & nb(kk,1)<a1(:,3) & a1(:,3)<nb(kk+1,1) & a1(:,4)<208);
       a5 = find(az(:,2)==5 & nb(kk,1)<az(:,5) & az(:,5)<nb(kk+1,1)); %BB
       a6 = find(a1(:,2)==6 & nb(kk,1)<a1(:,5) & a1(:,5)<nb(kk+1,1)); %Na
       a2=vertcat(a211,a221);
       y_dir=a1(a2,4);
       x_dir=a1(a2,3);
       X=[x_dir,y_dir];
       rng(5);
       GMModel = fitgmdist(X,2,'Replicates',5,'RegularizationValue',0.2);
       aaa=aaa+1;
       cent(:,:)=GMModel.mu; % get center 
       cent=sortrows(cent);
         d(dist) = pdist(cent,'euclidean');
       cent_mean(1)=mean(cent(:,1));
       cent_mean(2)=mean(cent(:,2));

        for ij=1:length(a5)   
        loc(1,ij)=az(a5(ij),3)-cent_mean(1);
        loc(2,ij)=az(a5(ij),4)-cent_mean(2);

poly=polyfit(cent(:,1),cent(:,2),1);
theta0=atan(poly(1));
theta=-theta0;
%computation of rotation matrix. rotation is about X axis
rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
Xrot=rot*loc;   
Xmean=mean(Xrot,2); 
Xcent=Xrot;
% plot
%plot(cent1(1,:),cent1(2,:),Xrot(1,:),Xrot(2,:))
        end
Xtime=horzcat(Xtime,Xcent);
   end
Xall=horzcat(Xall,Xtime);

end 
%% contour plot
data_point=101;
data_point_n=(data_point+1)/2;
snap_size=90;
pts = linspace(-snap_size, snap_size, data_point);
hitmap_box_size=abs(pts(2)-pts(1));

N=histcounts2(Xall(2,:),Xall(1,:),pts,pts);
N(data_point_n,data_point_n)=N(data_point_n,data_point_n)-(n+n1+n2);
imagesc(pts, pts, N/(hitmap_box_size*hitmap_box_size*(zhi-zlo)*(n+n1+n2)*33));
axis('square')
colorbar