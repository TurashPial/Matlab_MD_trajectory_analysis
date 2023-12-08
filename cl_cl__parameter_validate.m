clear all;
%variables
r= linspace(.0001,13,10000);
eps_cl=.0127850;
r_mn_cl=2.711*2;
eps_na=.3526418;
r_mn_na=1.212*2;
eps_nacl=0.067145;
r_mn_nacl=1.961362122*2;
eps_clo=0.23406;
r_mn_clo=2.988835;
eps_nao=0.0445668;
r_mn_nao=4.487603;
qcl=-1;
qna=1;
qo=-0.8476;


%% Cl
for i=1:size(r,2)
U_v_cl(i)=eps_cl*((r_mn_cl./r(i))^12-2*(r_mn_cl./r(i))^6);
end

for i=1:size(r,2)
U_e_cl(i)= 332.06371*qcl*qcl./r(i).*(1-erf(.2262*r(i)));
end
U_cl=U_v_cl+U_e_cl;
f_cl=-diff(U_cl)/.0013;
f_cl(10000)=0;
Cl_Cl=[1:10000; r; U_cl; f_cl] ;
CC = fopen('Cl_Cl.txt','w');
fprintf(CC, '%1s \n\n', '# Pair potential lj/cut/coul/long for atom types Cl Cl : i,r,energy,forcel');
fprintf(CC, '%1s \n', 'lj_cl');
fprintf(CC,'%1s %1s %1s %1s %1s \n\n','N','10000', 'R', '0.0001', '13.0');
fprintf(CC,'%1u %6.12f %12.12e %12.12e \r\n',Cl_Cl);
 
%% Na
for i=1:size(r,2)
U_v_na(i)=eps_na*((r_mn_na./r(i)).^12-2*(r_mn_na./r(i)).^6);
end
for i=1:size(r,2)
U_e_na(i)= 332.06371*qna*qna./r(i).*(1-erf(.2262*r(i)));
end
U_na=U_v_na+U_e_na;
f_na=-diff(U_na)/.0013;
f_na(10000)=0;
Na_Na=[1:10000; r; U_na; f_na] ;
NN = fopen('Na_Na.txt','w');
fprintf(NN, '%1s \n\n', '# Pair potential lj/cut/coul/long for atom types Na Na : i,r,energy,forcel');
fprintf(NN, '%1s \n', 'lj_na');
fprintf(NN,'%1s %1s %1s %1s %1s \n\n','N','10000', 'R', '0.0001', '13.0');
fprintf(NN,'%1u %6.12f %12.12e %12.12e \r\n',Na_Na);

%% nacl
for i=1:size(r,2)
U_v_nacl(i)=eps_nacl*((r_mn_nacl./r(i))^12-2*(r_mn_nacl./r(i))^6);
end

for i=1:size(r,2)
U_e_nacl(i)= 332.06371*qcl*qna./r(i).*(1-erf(.2262*r(i)));
end
U_nacl=U_v_nacl+U_e_nacl;
f_nacl=-diff(U_nacl)/.0013;
f_nacl(10000)=0;
Na_Cl=[1:10000; r; U_nacl; f_nacl] ;
NC = fopen('Na_Cl.txt','w');
fprintf(NC, '%1s \n\n', '# Pair potential lj/cut/coul/long for atom types Na Cl : i,r,energy,forcel');
fprintf(NC, '%1s \n', 'lj_nacl');
fprintf(NC,'%1s %1s %1s %1s %1s \n\n','N','10000', 'R', '0.0001', '13.0');
fprintf(NC,'%1u %6.12f %12.12e %12.12e \r\n',Na_Cl);

%% nao
for i=1:size(r,2)
U_v_nao(i)=eps_nao*((r_mn_nao./r(i))^12-2*(r_mn_nao./r(i))^6);
end

for i=1:size(r,2)
U_e_nao(i)= 332.06371*qna*qo./r(i).*(1-erf(.2262*r(i)));
end
U_nao=U_v_nao+U_e_nao;
f_nao=-diff(U_nao)/.0013;
f_nao(10000)=0;
Na_O=[1:10000; r; U_nao; f_nao] ;
NO = fopen('Na_O.txt','w');
fprintf(NO, '%1s \n\n', '# Pair potential lj/cut/coul/long for atom types Na O : i,r,energy,forcel');
fprintf(NO, '%1s \n', 'lj_nao');
fprintf(NO,'%1s %1s %1s %1s %1s \n\n','N','10000', 'R', '0.0001', '13.0');
fprintf(NO,'%1u %6.12f %12.12e %12.12e \r\n',Na_O);
%% nacl
for i=1:size(r,2) 
U_v_clo(i)=eps_clo*((r_mn_clo./r(i))^12-2*(r_mn_clo./r(i))^6);
end

for i=1:size(r,2)
U_e_clo(i)= 332.06371*qcl*qo./r(i).*(1-erf(.2262*r(i)));
end
U_clo=U_v_clo+U_e_clo;
f_clo=-diff(U_clo)/.0013;
f_clo(10000)=0;
Cl_O=[1:10000; r; U_clo; f_clo] ;
CO = fopen('Cl_O.txt','w');
fprintf(CO, '%1s \n\n', '# Pair potential lj/cut/coul/long for atom types Cl O : i,r,energy,forcel');
fprintf(CO, '%1s \n', 'lj_clo');
fprintf(CO,'%1s %1s %1s %1s %1s \n\n','N','10000', 'R', '0.0001', '13.0');
fprintf(CO,'%1u %6.12f %12.12e %12.12e \r\n',Cl_O);
