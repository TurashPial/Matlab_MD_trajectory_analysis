%-----------residance time of water in Na solvation shell----
clear all
tic
n=1000;%Provide even number
fs1=5;
fs2=4;
fhs=2.4; %cutoff
f_per=1;
dt=.05; %In fs
%% read simulation data
fid = fopen('coords_all_1V_12_res.dump','r');
box1 = textscan(fid,'%f %f','Headerlines',5);
D(1,:) = textscan(fid,'%f %f %f %f %f','Headerlines',9); 
fclose(fid);
box = cell2mat(box1);
xlo = box(1,1);xhi = box(1,2);
ylo = box(2,1);yhi = box(2,2);
zlo = box(3,1);zhi = box(3,2);

fid = fopen('coords_all_1V_12_res.dump','r');
for ii= 1:n
    C(ii,:) = textscan(fid,'%f %f %f %f %f','Headerlines',9);
end
fclose(fid);

%% initialization 
n_origin=1000;
z_graft = 44;
z_brush = 88;
Nr=[];
ni=1;
a22=[];
a221=[];
a222=[];
a23=[];
%%
for io=1:ni
    a1=cell2mat(C(io,:));
    io
    na_o=[];
a221 = find(a1(:,2)==9 & a1(:,5)>z_graft+fs1 & a1(:,5)<z_brush-fs1 & a1(:,3)<xhi-10 & a1(:,3)>xlo+10 & a1(:,4)<yhi-10 & a1(:,4)>ylo+10);
%a222 = find(a1(:,2)==10 & a1(:,5)>z_graft+fs1 & a1(:,5)<z_brush-fs1 & a1(:,3)<xhi-40 & a1(:,3)>xlo+30 & a1(:,4)<yhi-15 & a1(:,4)>ylo+15);
%a22=vertcat(a221,a222);
a22=a221;
a23 = find(a1(:,2)==5 & a1(:,5)>z_graft+fs1-fs2 & a1(:,5)<z_brush-fs1+fs2 & ...
    a1(:,3)<xhi-10+fs2 & a1(:,3)>xlo+10-fs2 & a1(:,4)<yhi-10+fs2 & a1(:,4)>ylo+10-fs2);
 for ii=1:length(a22)
       k=a22(ii);
       x_ion=a1(k,3);
       y_ion=a1(k,4);
       z_ion=a1(k,5);
              n=0;
        for jj=1:length(a23)
       w=a23(jj);
       x_water=a1(w,3);
       y_water=a1(w,4);
       z_water=a1(w,5);
       
       distance=sqrt((x_ion-x_water)^2+(y_ion-y_water)^2+(z_ion-z_water)^2);

if distance<fhs
    w_o=a23(jj);
    n=n+1;
    na_o(ii,n)=w_o;
end
        end
 end
 O_n=nnz(na_o);
 nh=length(min(na_o));
 Heavyside_step=[]; %heavyside step function
for t=1:n_origin-ni
    initial=cell2mat(C(t+io-1,:));
            for i=1:length(a22)
             i_index=a22(i);
             for nw=1:nh
                w_index=na_o(i,nw);
                if w_index~=0
    distance=(initial(i_index,3)-initial(w_index,3))^2+(initial(i_index,4)-initial(w_index,4))^2 +(initial(i_index,5)-initial(w_index,5))^2;
                else 
                    distance=12;
                end
                    
       if distance<(fhs*fhs)
               Heavyside_step(i,nw)=1;
            else
                Heavyside_step(i,nw)=0;
       end
             end
            end
            Nr(t,io)=sum(sum(Heavyside_step))/O_n;          
end
end
%%
time=linspace(dt,(n_origin-ni)*dt,(n_origin-ni));
nnr=mean(Nr,2);
plot(time,nnr/nnr(1));
toc