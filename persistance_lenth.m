clear all
tic
n = 8000;
%% read data file
fid = fopen('coords_poly_nvt_n10.dump','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
fclose(fid);
box = cell2mat(box1);
xlo = box(1,1);xhi = box(1,2);
ylo = box(2,1);yhi = box(2,2);
zlo = box(3,1);zhi = box(3,2);

fid = fopen('coords_poly_nvt_n10.dump','r');
for ii= 1:n
    C(ii,:) = textscan(fid,'%f %f %f %f %f %f','Headerlines',9);
end

%C contain 1-id 2-x 3-y 4-z 5-vx 6-vy 7-vz
fclose(fid);
%% calculate distribuion
for jj = 1:n
 a1 = cell2mat(C(jj,:));
 for ii=1:36
     b_c=find(a1(:,1)==(ii-1)*197+185);
     t_c=find(a1(:,1)==(ii-1)*197+194);
     z(ii,1) = a1(t_c,6)-a1(b_c,6);
     x(ii,1) = abs(a1(t_c,4)-a1(b_c,4));
     if x(ii,1)>45
         x(ii,1)=xhi-xlo-x(ii,1);
     end
     y(ii,1) = abs(a1(t_c,5)-a1(b_c,5));
     if y(ii,1)>45
         y(ii,1)=yhi-ylo-y(ii,1);
     end
     dist_i(ii,1)=z(ii,1)^2+x(ii,1)^2+y(ii,1)^2;

 end
 dist(jj,1)=mean(dist_i);

end
distance=mean(dist)
%%
p=10:.1:30;
L=73.392;
for lp=1:size(p,2)
    res(lp)=2*p(lp)*L*(1-p(lp)/L*(1-exp(-L/p(lp))))-distance;
end
plot(p,res)
%% solve for persistance length
ress=@(p) 2*p*L*(1-p/L*(1-exp(-L/p)))-distance;
a=10;
b=20;
iter=1;
x1=abs(a-b);
err=abs(a-b);
if ress(a)*ress(b)>0 
    disp('wrong initial range')
else
    pp = (a + b)/2;
    
    while err > 1e-4 
        iter=iter+1;
   if ress(a)*ress(pp)<0 
       b = pp;
   else
       a = pp;          
   end
    pp = (a + b)/2; 
    x1(iter)=pp;
    err = abs(x1(iter)-x1(iter-1));
    disp([num2str(iter),' : ',num2str(pp),' : ', num2str(err)]);
    end
    
end

    