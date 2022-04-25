clear
t_a=cputime;
% ----------------Generating and numbering mesh-------------%
model = createpde(1);
a=importGeometry(model,'torque_arm.stl');
mesh= generateMesh (model, 'GeometricOrder', 'linear','Hmax',0.75, 'Hmin',0.5);
plot(mesh.Nodes(1,:), ...
     mesh.Nodes(2,:), ...
     'ok','markersize',4,'MarkerFaceColor','g')
hold on
pdemesh(mesh,'NodeLabels','on')
pdeplot(model);
xlim([-10 55])
ylim([-5 25])


%Propeties of the material 
E=200*10^3; %N/mm^2
nu=0.3;
t=10; %mm 

%Plane Strain material property matrix 
D_pstrain=(E/((1+nu)*(1-2*nu)))*[(1-nu) nu 0; nu (1-nu) 0;0 0 ((1-2*nu)/2)];
% D_pstress=(E/(1-nu^2))*[1 nu 0; nu 1 0;0 0 ((1-nu)/2)];

L_elmt=length(mesh.Elements);
L_nd=length(mesh.Nodes);
F=zeros((2*L_nd),1);
K_stress=zeros((2*L_nd),(2*L_nd));

for i=1:(L_elmt)
n_1=mesh.Elements(1,i); 
n_2=mesh.Elements(2,i); 
n_3=mesh.Elements(3,i);

%Values are multiplied by 10 to convert the coordinates from [cm]to [mm]
x1=mesh.Nodes(1,n_1)*10;
y1=mesh.Nodes(2,n_1)*10;
x2=mesh.Nodes(1,n_2)*10;
y2=mesh.Nodes(2,n_2)*10;
x3=mesh.Nodes(1,n_3)*10;
y3=mesh.Nodes(2,n_3)*10;
A_1=[1 x1 y1;
    1 x2 y2;
    1 x3 y3];
Area_1=(1/2)*det(A_1);

%Element Strain Displacement Matrix 
d_psi=(1/(2*Area_1))*[(y2-y3) 0 (y3-y1) 0 (y1-y2) 0;
    0 (x3-x2) 0 (x1-x3) 0 (x2-x1); 
    (x3-x2) (y2-y3) (x1-x3) (y3-y1) (x2-x1) (y1-y2)];

B_1{i,1}=d_psi;
% G_strain{i,1}=t*Area*transpose(d_psi)*D_pstrain*d_psi;
G_stress=t*(Area_1)*transpose(d_psi)*D_pstrain*d_psi;

K_stress(2*(n_1-1)+1,2*(n_1-1)+1)=G_stress(1,1)+K_stress(2*(n_1-1)+1,2*(n_1-1)+1);
K_stress(2*(n_1-1)+1,2*(n_1-1)+2)=G_stress(1,2)+K_stress(2*(n_1-1)+1,2*(n_1-1)+2);
K_stress(2*(n_1-1)+1,2*(n_2-1)+1)=G_stress(1,3)+K_stress(2*(n_1-1)+1,2*(n_2-1)+1);
K_stress(2*(n_1-1)+1,2*(n_2-1)+2)=G_stress(1,4)+K_stress(2*(n_1-1)+1,2*(n_2-1)+2);
K_stress(2*(n_1-1)+1,2*(n_3-1)+1)=G_stress(1,5)+K_stress(2*(n_1-1)+1,2*(n_3-1)+1);
K_stress(2*(n_1-1)+1,2*(n_3-1)+2)=G_stress(1,6)+K_stress(2*(n_1-1)+1,2*(n_3-1)+2);
K_stress(2*(n_1-1)+2,2*(n_1-1)+1)=G_stress(2,1)+K_stress(2*(n_1-1)+2,2*(n_1-1)+1);
K_stress(2*(n_1-1)+2,2*(n_1-1)+2)=G_stress(2,2)+K_stress(2*(n_1-1)+2,2*(n_1-1)+2);
K_stress(2*(n_1-1)+2,2*(n_2-1)+1)=G_stress(2,3)+K_stress(2*(n_1-1)+2,2*(n_2-1)+1);
K_stress(2*(n_1-1)+2,2*(n_2-1)+2)=G_stress(2,4)+K_stress(2*(n_1-1)+2,2*(n_2-1)+2);
K_stress(2*(n_1-1)+2,2*(n_3-1)+1)=G_stress(2,5)+K_stress(2*(n_1-1)+2,2*(n_3-1)+1);
K_stress(2*(n_1-1)+2,2*(n_3-1)+2)=G_stress(2,6)+K_stress(2*(n_1-1)+2,2*(n_3-1)+2);
K_stress(2*(n_2-1)+1,2*(n_1-1)+1)=G_stress(3,1)+K_stress(2*(n_2-1)+1,2*(n_1-1)+1);
K_stress(2*(n_2-1)+1,2*(n_1-1)+2)=G_stress(3,2)+K_stress(2*(n_2-1)+1,2*(n_1-1)+2);
K_stress(2*(n_2-1)+1,2*(n_2-1)+1)=G_stress(3,3)+K_stress(2*(n_2-1)+1,2*(n_2-1)+1);
K_stress(2*(n_2-1)+1,2*(n_2-1)+2)=G_stress(3,4)+K_stress(2*(n_2-1)+1,2*(n_2-1)+2);
K_stress(2*(n_2-1)+1,2*(n_3-1)+1)=G_stress(3,5)+K_stress(2*(n_2-1)+1,2*(n_3-1)+1);
K_stress(2*(n_2-1)+1,2*(n_3-1)+2)=G_stress(3,6)+K_stress(2*(n_2-1)+1,2*(n_3-1)+2);
K_stress(2*(n_2-1)+2,2*(n_1-1)+1)=G_stress(4,1)+K_stress(2*(n_2-1)+2,2*(n_1-1)+1);
K_stress(2*(n_2-1)+2,2*(n_1-1)+2)=G_stress(4,2)+K_stress(2*(n_2-1)+2,2*(n_1-1)+2);
K_stress(2*(n_2-1)+2,2*(n_2-1)+1)=G_stress(4,3)+K_stress(2*(n_2-1)+2,2*(n_2-1)+1);
K_stress(2*(n_2-1)+2,2*(n_2-1)+2)=G_stress(4,4)+K_stress(2*(n_2-1)+2,2*(n_2-1)+2);
K_stress(2*(n_2-1)+2,2*(n_3-1)+1)=G_stress(4,5)+K_stress(2*(n_2-1)+2,2*(n_3-1)+1);
K_stress(2*(n_2-1)+2,2*(n_3-1)+2)=G_stress(4,6)+K_stress(2*(n_2-1)+2,2*(n_3-1)+2);
K_stress(2*(n_3-1)+1,2*(n_1-1)+1)=G_stress(5,1)+K_stress(2*(n_3-1)+1,2*(n_1-1)+1);
K_stress(2*(n_3-1)+1,2*(n_1-1)+2)=G_stress(5,2)+K_stress(2*(n_3-1)+1,2*(n_1-1)+2);
K_stress(2*(n_3-1)+1,2*(n_2-1)+1)=G_stress(5,3)+K_stress(2*(n_3-1)+1,2*(n_2-1)+1);
K_stress(2*(n_3-1)+1,2*(n_2-1)+2)=G_stress(5,4)+K_stress(2*(n_3-1)+1,2*(n_2-1)+2);
K_stress(2*(n_3-1)+1,2*(n_3-1)+1)=G_stress(5,5)+K_stress(2*(n_3-1)+1,2*(n_3-1)+1);
K_stress(2*(n_3-1)+1,2*(n_3-1)+2)=G_stress(5,6)+K_stress(2*(n_3-1)+1,2*(n_3-1)+2);
K_stress(2*(n_3-1)+2,2*(n_1-1)+1)=G_stress(6,1)+K_stress(2*(n_3-1)+2,2*(n_1-1)+1);
K_stress(2*(n_3-1)+2,2*(n_1-1)+2)=G_stress(6,2)+K_stress(2*(n_3-1)+2,2*(n_1-1)+2);
K_stress(2*(n_3-1)+2,2*(n_2-1)+1)=G_stress(6,3)+K_stress(2*(n_3-1)+2,2*(n_2-1)+1);
K_stress(2*(n_3-1)+2,2*(n_2-1)+2)=G_stress(6,4)+K_stress(2*(n_3-1)+2,2*(n_2-1)+2);
K_stress(2*(n_3-1)+2,2*(n_3-1)+1)=G_stress(6,5)+K_stress(2*(n_3-1)+2,2*(n_3-1)+1);
K_stress(2*(n_3-1)+2,2*(n_3-1)+2)=G_stress(6,6)+K_stress(2*(n_3-1)+2,2*(n_3-1)+2);
end 

%Force is assigned to node
F(4,1)=5000; %N
F(477,1)=2800; %N

%Condense form is created according to boundary conditions 
zero=535:1:600;
K_stress(zero,:)=[];
K_stress(:,zero)=[];
K_stress(:,5)=[];
K_stress(5,:)=[];
K_stress(:,6)=[];
K_stress(6,:)=[];
F(5,:)=[];
F(6,:)=[];
F(zero,:)=[];

U=inv(K_stress)*F;

for i=1:2:1730
Ux(i,1)=U(i,1);
end

for i=2:2:1730
Uy(i,1)=U(i,1);
end 

MaxUx=max(Ux);
MaxUy=max(Uy);
fprintf('Maximum x displacement is %.4f [mm]\n',MaxUx)
fprintf('Maximum y displacement is %.4f [mm]\n',MaxUy)

U1=zeros((length(zero)),1);
U=[U(1:4,1);0; 0;U(5:1730,1)];
U=[U(1:534,1);U1;U(535:1732,1)];

for i=1:(L_elmt)
n_1=mesh.Elements(1,i); 
n_2=mesh.Elements(2,i); 
n_3=mesh.Elements(3,i);

C_1=B_1{i,1};
x1=U(2*(n_1-1)+1,1);
y1=U(2*(n_1-1)+2,1);
x2=U(2*(n_2-1)+1,1);
y2=U(2*(n_2-1)+2,1);
x3=U(2*(n_3-1)+1,1);
y3=U(2*(n_3-1)+2,1);
U2=[x1;y1;x2;y2;x3;y3];

%Stress is calculated 
Stress=D_pstrain*C_1*U2;

%Von misses is calculated for each element
vonmisses(i,1)=sqrt((Stress(1,1)^2-Stress(1,1)*Stress(2,1)-Stress(2,1)^2)+3*Stress(3,1)^2);

% Values are multiply by 10, to turn the coordinates [cm] to [mm]
x1=mesh.Nodes(1,n_1)*10;
y1=mesh.Nodes(2,n_1)*10;
x2=mesh.Nodes(1,n_2)*10;
y2=mesh.Nodes(2,n_2)*10;
x3=mesh.Nodes(1,n_3)*10;
y3=mesh.Nodes(2,n_3)*10;

% Deformation Values (U,V) are multiplied with 44 to transform True
% Scale to Auto Scale(2xAuto). This info taken From ANSYS.
x_1=(x1+44*(U(2*(n_1-1)+1,1)))/10; 
x_2=(x2+44*(U(2*(n_2-1)+1,1)))/10;
x_3=(x3+44*(U(2*(n_3-1)+1,1)))/10;
y_1=(y1+44*(U(2*(n_1-1)+2,1)))/10;
y_2=(y2+44*(U(2*(n_2-1)+2,1)))/10;
y_3=(y3+44*(U(2*(n_3-1)+2,1)))/10;
Y=[y_1 y_2 y_3 y_1];
X=[x_1 x_2 x_3 x_1];

figure(2) 
plot(X,Y,'r-o','linewidth',2,'markersize',5,'markeredgecolor','k','markerfacecolor','g')
xlabel('X[cm]', 'fontsize',18)
ylabel('Y[cm]', 'fontsize',18)
title('Deformed Geometry', 'fontsize',18)
xlim([-10 55])
ylim([-5 25])
hold on 
end 

Maxvon=max(vonmisses);
fprintf('Maximum VonMisses Stress %.4f [MPa]',Maxvon)
cputime_1=cputime-t_a