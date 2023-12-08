%--------create data file of water with bond definations, simulate with CNT 
clc
clear all
r = 1.0;
theta = pi*(109.47/180) ;
x_O = -2.50;
x_H_1 = -2.50;
x_H_2 = -2.50;

y_O = -2.50;
y_H_1 = y_O+r*sin(theta/2);
y_H_2 = y_O-r*sin(theta/2);

z_O = 3.0000;
z_H_1 = z_O+r*cos(theta/2);
z_H_2 = z_O+r*cos(theta/2);

xv = 0.0;
yv = 0.0;
zv =  0.0:2.1:195;
X_O = x_O + xv;
Y_O = y_O + yv;
Z_O = z_O + zv;

X_H_1 = x_H_1 + xv;
Y_H_1 = y_H_1 + yv;
Z_H_1 = z_H_1 + zv;

X_H_2 = x_H_2 + xv;
Y_H_2 = y_H_2 + yv;
Z_H_2 = z_H_2 + zv;

W0 = [];
for i = 1:length(X_O)
    for j = 1:length(Y_O)
        for k = 1:length(Z_O)
            W0 = vertcat(W0,[X_O(i) Y_O(j) Z_O(k);X_H_1(i) Y_H_1(j) Z_H_1(k);X_H_2(i) Y_H_2(j) Z_H_2(k)]) ;
        end
    end
end




r = 1.0;
theta = pi*(109.47/180) ;
x_O = 2.500;
x_H_1 = 2.500;
x_H_2 = 2.500;

y_O = -2.500;
y_H_1 = y_O+r*sin(theta/2);
y_H_2 = y_O-r*sin(theta/2);

z_O = 3.0000;
z_H_1 = z_O+r*cos(theta/2);
z_H_2 = z_O+r*cos(theta/2);

xv = 0.0;
yv = 0.0;
zv = 0.0:2.1:195;
X_O = x_O + xv;
Y_O = y_O + yv;
Z_O = z_O + zv;

X_H_1 = x_H_1 + xv;
Y_H_1 = y_H_1 + yv;
Z_H_1 = z_H_1 + zv;

X_H_2 = x_H_2 + xv;
Y_H_2 = y_H_2 + yv;
Z_H_2 = z_H_2 + zv;

W1 = [];
for i = 1:length(X_O)
    for j = 1:length(Y_O)
        for k = 1:length(Z_O)
            W1 = vertcat(W1,[X_O(i) Y_O(j) Z_O(k);X_H_1(i) Y_H_1(j) Z_H_1(k);X_H_2(i) Y_H_2(j) Z_H_2(k)]) ;
        end
    end
end




r = 1.0;
theta = pi*(109.47/180) ;
x_O = -2.500;
x_H_1 = -2.500;
x_H_2 = -2.500;

y_O = 2.500;
y_H_1 = y_O+r*sin(theta/2);
y_H_2 = y_O-r*sin(theta/2);

z_O = 3.0000;
z_H_1 = z_O+r*cos(theta/2);
z_H_2 = z_O+r*cos(theta/2);

xv = 0.0;
yv = 0.0;
zv = 0.0:2.1:195;
X_O = x_O + xv;
Y_O = y_O + yv;
Z_O = z_O + zv;

X_H_1 = x_H_1 + xv;
Y_H_1 = y_H_1 + yv;
Z_H_1 = z_H_1 + zv;

X_H_2 = x_H_2 + xv;
Y_H_2 = y_H_2 + yv;
Z_H_2 = z_H_2 + zv;

W2 = [];
for i = 1:length(X_O)
    for j = 1:length(Y_O)
        for k = 1:length(Z_O)
            W2 = vertcat(W2,[X_O(i) Y_O(j) Z_O(k);X_H_1(i) Y_H_1(j) Z_H_1(k);X_H_2(i) Y_H_2(j) Z_H_2(k)]) ;
        end
    end
end

r = 1.0;
theta = pi*(109.47/180) ;
x_O = 2.500;
x_H_1 = 2.500;
x_H_2 = 2.500;

y_O = 2.500;
y_H_1 = y_O+r*sin(theta/2);
y_H_2 = y_O-r*sin(theta/2);

z_O = 3.0000;
z_H_1 = z_O+r*cos(theta/2);
z_H_2 = z_O+r*cos(theta/2);

xv = 0.0;
yv = 0.0;
zv = 0.0:2.1:195;
X_O = x_O + xv;
Y_O = y_O + yv;
Z_O = z_O + zv;

X_H_1 = x_H_1 + xv;
Y_H_1 = y_H_1 + yv;
Z_H_1 = z_H_1 + zv;

X_H_2 = x_H_2 + xv;
Y_H_2 = y_H_2 + yv;
Z_H_2 = z_H_2 + zv;

W3 = [];
for i = 1:length(X_O)
    for j = 1:length(Y_O)
        for k = 1:length(Z_O)
            W3 = vertcat(W3,[X_O(i) Y_O(j) Z_O(k);X_H_1(i) Y_H_1(j) Z_H_1(k);X_H_2(i) Y_H_2(j) Z_H_2(k)]) ;
        end
    end
end

r = 1.0;
theta = pi*(109.47/180) ;
x_O = 0.00;
x_H_1 = 0.00;
x_H_2 = 0.00;

y_O = 0.00;
y_H_1 = y_O+r*sin(theta/2);
y_H_2 = y_O-r*sin(theta/2);

z_O = 3.0000;
z_H_1 = z_O+r*cos(theta/2);
z_H_2 = z_O+r*cos(theta/2);

xv = 0.0;
yv = 0.0;
zv =  0.0:2.1:195;
X_O = x_O + xv;
Y_O = y_O + yv;
Z_O = z_O + zv;

X_H_1 = x_H_1 + xv;
Y_H_1 = y_H_1 + yv;
Z_H_1 = z_H_1 + zv;

X_H_2 = x_H_2 + xv;
Y_H_2 = y_H_2 + yv;
Z_H_2 = z_H_2 + zv;

W4 = [];
for i = 1:length(X_O)
    for j = 1:length(Y_O)
        for k = 1:length(Z_O)
            W4 = vertcat(W4,[X_O(i) Y_O(j) Z_O(k);X_H_1(i) Y_H_1(j) Z_H_1(k);X_H_2(i) Y_H_2(j) Z_H_2(k)]) ;
        end
    end
end


W = [W0;W1;W2;W3;W4];
id = 3505:1:length(W(:,1))+3504;
a_type = [];
for i = 1:length(id)/3
    a_type = vertcat(a_type,[3;2;2]);
end
q = [];
for i = 1:length(id)/3
    q = vertcat(q,[-0.8476;0.4238;0.4238]);
end
m_ID = [];
for i = 1:length(id)/3
    m_ID = vertcat(m_ID,[i;i;i]);
   
end
m_ID = 1+m_ID;
W_atoms = [id' m_ID a_type q W(:,1) W(:,2) W(:,3)];

bond_id = 1:1:2*length(id)/3 ;

b_O = [];
b_H = [];
for i = 1:3:length(W_atoms(:,1))-2
    b_O = vertcat(b_O,[W_atoms(i,1);W_atoms(i,1)]);
    b_H = vertcat(b_H,[W_atoms(i+1,1);W_atoms(i+2,1)]);
end
bond = [bond_id' ones(length(bond_id),1) b_O b_H];

angle_id = 1:1:length(id)/3 ;
a_O = [];
for i = 1:3:length(W_atoms(:,1))-2
    a_O = vertcat(a_O,[W_atoms(i+1,1) W_atoms(i,1) W_atoms(i+2,1)]);
end
angle = [angle_id' ones(length(angle_id),1) a_O]; 


fileID = fopen('water6.txt','w');

fprintf(fileID,'Atoms\r\n\n');
for i = 1:length(W_atoms(:,1))
    
    fprintf(fileID,'%d %d %d %f %f %f %f\r\n',W_atoms(i,:));
end
fprintf(fileID,'\n');
fprintf(fileID,'Bonds\r\n\n');
for i = 1:length(bond)
    
    fprintf(fileID,'%d %d %d %d\r\n',bond(i,:));
end
fprintf(fileID,'\n');
fprintf(fileID,'Angles\r\n\n');
for i = 1:length(angle)
    
    fprintf(fileID,'%d %d %d %d %d\r\n',angle(i,:));
end

