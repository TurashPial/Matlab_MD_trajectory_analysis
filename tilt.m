%-----calculate how tilted are polymers -------------
clear all
tic
n = 200;
%% read simulation data
fid = fopen('coords_poly_1V_4th.dump','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
fclose(fid);
box = cell2mat(box1);
xlo = box(1,1);xhi = box(1,2);
ylo = box(2,1);yhi = box(2,2);
zlo = box(3,1);zhi = box(3,2);
%%
fid = fopen('coords_poly_1V_4th.dump','r');
for ii= 1:n
    C(ii,:) = textscan(fid,'%f %f %f %f %f %f','Headerlines',9);
end
%% define polymer chain backbone  
onechain=zeros(49,1);
 for i=1:23 
     onechain(2*i+1)=395+8*(i-1);
     onechain(2*i+2)=398+8*(i-1);
 end
onechain(1)=579;
onechain(2)=583;
onechain(49)=588;
lowerchains=zeros(49,54);
 for i=1:6
     for ii=1:9
     lowerchains(:,ii+9*(i-1))=onechain+197*(ii-1)+2955*(i-1);
     end
 end
 
 onechainU=zeros(49,1);
 for i=1:23
     onechainU(2*i+1)=18125+8*(i-1);
     onechainU(2*i+2)=18128+8*(i-1);
 end
onechainU(1)=18309;
onechainU(2)=18313;
onechainU(49)=18318;
upperchains=zeros(49,54);
 for i=1:6
     for ii=1:9
     upperchains(:,ii+9*(i-1))=onechainU+197*(ii-1)+2955*(i-1);
     end
 end

%% average for all polymers, all time steps
clear averagex averagez
for jj = 1:n
   a1 = cell2mat(C(jj,:));
   for i=1:54
       for ii=1:49
       back=find(a1(:,1)==upperchains(ii,i));
       x_coordu(ii,i)=a1(back,4);
       z_coordu(ii,i)=a1(back,6);
       end
       end
      averagexu(:,jj)=mean(x_coordu,2);
      averagezu(:,jj)=mean(z_coordu,2);
      
   for i=1:54
       for ii=1:49
       backl=find(a1(:,1)==lowerchains(ii,i));
       x_coordl(ii,i)=a1(backl,4);
       z_coordl(ii,i)=a1(backl,6);
       end
       end
      averagexl(:,jj)=mean(x_coordl,2);
      averagezl(:,jj)=mean(z_coordl,2);
      
display (jj)
end
%% plot polymer profile
figure
tavgxu=mean(averagexu,2);
tavgzu=mean(averagezu,2);
upx=tavgxu-tavgxu(1);
upz=tavgzu*-1+tavgzu(1);
plot(upx,upz)
figure
tavgxl=mean(averagexl,2);
tavgzl=mean(averagezl,2);
lox=tavgxl-tavgxl(1);
loz=tavgzl-tavgzl(1);
plot(lox,loz)
%%
figure
x=(upx+lox)./2;
z=(upz+loz)./2;
plot(x,z)
toc
% z_graft_bot=41.3231; % check in ovito
% z_graft_top=155.743; % check in ovito