%-------------flow rate----------
a = (-0.7502-0.48293)/2; % average of velocity at points 10 and 15 from 
%0.1V case of velocity from ACS nano
u_es1 = (2.38e-3)/0.1*(-0.6166) ;%gd 0.05 f1 70_82_ns
u_es2 = (3.68e-3)/0.1*(-0.6166) ;%gd 0.05 f2 55_67_ns
u_es3 = (4.34e-3)/0.1*(-0.6166) ;%gd 0.03 f1 17_29_ns
%%
u_T_f1_70_82 = u_pp_f1_70_82+u_EOS_f1_70_82;
u_T_f2_55_67 = u_pp_f2_55_67+u_EOS_f2_55_67;
u_T_f1_17_29 = u_pp_f1_17_29+u_EOS_f1_17_29;
%%
Q_c1_EOS = trapz(y_EOS_f1_70_82,y_EOS_f1_70_82)*1e-10*93.9e-10;
Q_c2_EOS = trapz(y_EOS_f2_55_67,y_EOS_f2_55_67)*1e-10*93.9e-10;
Q_c3_EOS = trapz(y_EOS_f1_17_29,y_EOS_f1_17_29)*1e-10*101.007e-10;
Q_c1_pp = trapz(y_pp_f1_70_82,u_pp_f1_70_82)*1e-10*93.9e-10;
Q_c2_pp = trapz(y_pp_f2_55_67,u_pp_f2_55_67)*1e-10*93.9e-10;
Q_c3_pp = trapz(y_pp_f1_17_29,u_pp_f1_17_29)*1e-10*101.007e-10;
Q_c1_T = trapz(y_EOS_f1_70_82,u_T_f1_70_82)*1e-10*93.9e-10;
Q_c2_T = trapz(y_EOS_f2_55_67,u_T_f2_55_67)*1e-10*93.9e-10;
Q_c3_T = trapz(y_EOS_f1_17_29,u_T_f1_17_29)*1e-10*101.007e-10;
Q_2 = trapz(y_pp_76_88_ns_4,u_pp_76_88_ns_4)*1e-10*93.9e-10;
%%
clear all
load('velocity_gd_05_f1_76_88_ns.mat')
x1 = u; y1 = y;
t = tiledlayout(1,1);
ax1 = axes(t);
plot(ax1,x1,y1,'-b^')
ax1.XColor = 'b';
ax1.YColor = 'k';
ax1.Box = 'off';
%%
load('num_gd_05_f1_76_88_ns.mat')
y2 = y;
x2 = num;
ax2 = axes(t);
plot(ax2,x2*1e3,y2,'-ro')
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.XColor = 'r';
ax2.YColor = 'k';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';