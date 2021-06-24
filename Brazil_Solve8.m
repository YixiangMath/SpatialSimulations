
%%Created on 11/6/2020
%%%Standard incidence
%%%Four equations
%%%Works well with start date March 17, 2020
%%This file solves diffusive coronavirus model using finite element scheme
%%%Initial population data will be used by world data

%%%%%%%%%%%%%%%%%%%%%%%%

%%Run mathematica file  Brazil_Border to get a simple boundary
%%The data will be written in Lat.txt (Latitude data) and Log.txt (Longitude data)
%%Import the data and store them in Brazil_LAT and Brazil_LOG
clear all
clc

fileID=fopen('SimpleLat.txt');
C=textscan(fileID,'%f');
Brazil_LAT=cell2mat(C);
fclose(fileID);


fileID=fopen('SimpleLog.txt');
C=textscan(fileID,'%f');
Brazil_LOG=cell2mat(C);
fclose(fileID);

Num=length(Brazil_LAT);

%%Longitude and Latitude are the boundary data
Longitude =Brazil_LOG(1:1:Num);
Latitude =Brazil_LAT(1:1:Num);


%%Convert Latitude and Longitude into distance with Unit Km (Kilometers)
%%See wikipedia for the formula
e=sqrt(0.00669437999014);
a=6378137.0;

phi=mean(Latitude)/180*pi;
%phi=Latitude/180*pi;

DeltaLatitude=pi*a*(1-e^2)./(180*( 1-e^2 * (sin(phi)).^2   ).^1.5     )/1000; %% 1 Latitude= 1 DeltaLatitude Km
DeltaLongitude= pi*a*cos(phi)./  (180*( 1-e^2 * (sin(phi)).^2   ).^0.5     )/1000; 


ne=length(Longitude);   %% ne = the number of edeges of polygons
X= (Longitude+50).*DeltaLongitude; %% X=Longitude in Km
Y=(Latitude+12).*DeltaLatitude;  %%Y=Latitude in Km

%[X,Y] = grn2eqa(Longitude,Latitude);

%%%%%%%%%%%%%%%%%%%%%%%%
%%Create the geometry

g=decsg([2, ne, X', Y']'); %% Create polygon
pdem=createpde(4) %% 3 is the number of equations of the model

geometryFromEdges(pdem, g);  %% Combine the model with the golygon 
hmax=10; %% Mesh size
%hmax=0.01;
generateMesh(pdem, 'Hmax', hmax);   %%% Generate the mesh 

% h=figure;
% pdeplot(pdem) %% Plot mesh
% saveas(h, 'Mesh', 'epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Solve the model
applyBoundaryCondition(pdem,'neumann', 'Edge', 1:ne, 'q', zeros(3,3),'g', zeros(3,1))

%%Specify coefficients for the SI model
%%c=diffusion rates
%%Matlab must specify a boundary condition, so 
%%set the diffusion rate of S to be a small number

c=[0; 0; 2; 2; 1; 1; 0; 0];

a=0;%%f=char('-0.02*u(1).*u(2)./(1+0.05*u(2))','0.02*u(1).*u(2)./(1+0.05*u(2))-3.5*u(2)');
d=1;

specifyCoefficients(pdem, 'm',0,'d',1,'c', c,'a',0,'f',@fcoeffunction)
pdem.SolverOptions.AbsoluteTolerance=10^1;
pdem.SolverOptions.RelativeTolerance=10^0;
pdem.SolverOptions.ResidualTolerance=10^1;
%pdem.SolverOptions.ResidualNorm=2;


%%Initial condition 

%unit0=@(locations) unit(locations, DeltaLongitude, DeltaLatitude);
setInitialConditions(pdem, @unit1)


%pdem.SolverOptions.ResidualNorm=2;



tlist=0:1:30;   %%Solution time points

sol=solvepde(pdem,tlist);

%save('HMax10Second.mat')



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Right hand side of SI model
function f = fcoeffunction(region,state)

N = 4; % Number of equations
nr = length(region.x); % Number of columns
f = zeros(N,nr); % Allocate f
% tau1=0.05*10^-3;
% tau2=0.5*10^-2;
tau1=0.16;
tau2=0.15;
lambda=1/7;
gamma=1/7;

% Now the particular functional form of f
% f(1,:) = -state.u(1,:)./(state.u(1,:)+ state.u(2,:)+state.u(3,:)+0.0001 ).*(tau1*state.u(2,:)+tau2*state.u(3,:));
% f(2,:) = state.u(1,:)./(state.u(1,:)+ state.u(2,:)+state.u(3,:)+0.0001 ).*(tau1*state.u(2,:)+tau2*state.u(3,:))-lambda*state.u(2,:);
% f(3,:)=lambda*state.u(2,:)-gamma*state.u(3,:);
f(1,:) = -state.u(1,:)/(state.u(1,:)+ state.u(2,:)+state.u(3,:)+state.u(4,:))*(tau1*state.u(2,:)+tau2*state.u(3,:));
f(2,:) = state.u(1,:)/(state.u(1,:)+ state.u(2,:)+state.u(3,:)+state.u(4,:) )*(tau1*state.u(2,:)+tau2*state.u(3,:))-lambda*state.u(2,:);
f(3,:)=lambda*state.u(2,:)-gamma*state.u(3,:);
f(4,:)=gamma*state.u(3,:);
end


%%%%%%%%%%%%%%%%%%%%%%%
%%Initial data

%function init=unit(locations, DeltaLongitude, DeltaLatitude)
function init=unit1(locations)

%%%Import the world population data
[pop20151, R20151] = geotiffread('gpw-v4-population-density-2015/gpw-v4-population-density_2015.tif');

%%%Flip the matrix data
Density=zeros(17400, 43200);
for p=1:17400 
   Density(p,:)=pop20151(17400-p+1,:);
end

%%%%%%%%View the world population density%%%%%%%
% Ref=georasterref('RasterSize', size(pop20151), 'Latlim', [-90 90], 'Lonlim', [0 360]);
% figure
% grid2image(Density, Ref)


%%%%%In the data set, negative values = no one lives
for p=1:17400
for    q=1:43200
    if Density(p, q)<0
        Density(p, q)=0;
    end
end
end

   
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%


fileID=fopen('SimpleLat.txt');
C=textscan(fileID,'%f');
Brazil_LAT=cell2mat(C);
fclose(fileID);


Num=length(Brazil_LAT);

%%Longitude and Latitude are the boundary data
%Longitude =Brazil_LOG(1:1:Num);
Latitude =Brazil_LAT(1:1:Num);


%%Convert Latitude and Longitude into distance with Unit Km (Kilometers)
%%See wikipedia for the formula
e=sqrt(0.00669437999014);
a=6378137.0;

phi=mean(Latitude)/180*pi;
%phi=Latitude/180*pi;

DeltaLatitude=pi*a*(1-e^2)./(180*( 1-e^2 * (sin(phi)).^2   ).^1.5     )/1000; %% 1 Latitude= 1 DeltaLatitude Km
DeltaLongitude= pi*a*cos(phi)./  (180*( 1-e^2 * (sin(phi)).^2   ).^0.5     )/1000;     
    
    
%%%%%%%%%%%%%%%%%%%%%Initial values
StateArea=[164123,27779,142829,1559159,564733,148920,5780,46096,340112,331937,903366,357146,586522,1247955,56470,199308,98148,251578,43780,52811,281730,237591,224301,  95736,248223,21915,277721];
%%%constant state population
%Infected=[30,21,75,407,249,635,231,172,82,162,43,38,281,106,40,252,333,21,1017,121,247,17,34,265,1500,16,9];

%%hmax=20
%Infected=[10,50,200,420,350,720,500,350,0,162,15,60,120,50,80,220,300,5,1000,200,247,20,30,100,2500,5,20];


%%%Asymp=[15,15,100,250,180,600,200,100,50,162,25,25,200,50,30,150,300,10,1000,150,220,15,25,200,2500,25,5];
Asymp=[15,15,100,250,180,600,200,100,50,162,25,25,200,50,30,150,300,10,1000,150,220,15,25,200,2500,25,5];
Infected=5*[0, 0, 0, 1, 3, 5, 20, 2, 6, 0, 0, 4, 7, 0, 0, 6, 16, 0, 30, 1, 10, 0, 0, 7, 160, 4, 0];

%%%%%%State Population Density
fileID=fopen('StatesDensity.txt');
C=textscan(fileID,'%f');
StatesDensity=cell2mat(C);
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%Read state boundary file
BRShape1=shaperead('bra_adm_ibge_2020_shp/bra_admbnda_adm1_ibge_2020.shp');

%%%p=coordinates of notes
nodesX=locations.x;
nodesY=locations.y;

num=length(nodesX);
init=zeros(4, num);

for I=1:num
    %%%%%take a particlar node
   nodesXI=nodesX(I);
   nodesYI=nodesY(I);
    
    MeshLongitude= nodesXI/DeltaLongitude-50;
    MeshLatitude=nodesYI/DeltaLatitude-12;
    p1=round( (MeshLongitude+180)/0.0083333333333, 0);
    p2=round( (MeshLatitude+60)/0.0083333333333, 0);
    init(1, I)=Density(p2, p1);
end


StatePop=zeros(1,27);
for I=1:num
    %%%%%take a particlar node
   nodesXI=nodesX(I);
   nodesYI=nodesY(I);
    %%%%%check which state the node belongs to 
    Kmax=27; %%K=number of states
    for K=1:Kmax

       %%%%Take the coordinates of the K state
      StateX=BRShape1(K).X;
      StateY=BRShape1(K).Y;

      %%Convert Latitudes/Logitudes data into XY coordinates 
      StateX= (StateX+50).*DeltaLongitude; %% X=Longitude in Km
      StateY=(StateY+12).*DeltaLatitude;  %%Y=Latitude in Km
     CheckIn=inpolygon(nodesXI,nodesYI,StateX,StateY); %1=inside,0=outside
     if CheckIn==1
        StatePop(K)=StatePop(K)+init(1, I);
     end
    end
end



for I=1:num
    %%%%%take a particlar node
   nodesXI=nodesX(I);
   nodesYI=nodesY(I);
    %%%%%check which state the node belongs to 
    Kmax=27; %%K=number of states
    for K=1:Kmax

       %%%%Take the coordinates of the K state
      StateX=BRShape1(K).X;
      StateY=BRShape1(K).Y;

      %%Convert Latitudes/Logitudes data into XY coordinates 
      StateX= (StateX+50).*DeltaLongitude; %% X=Longitude in Km
      StateY=(StateY+12).*DeltaLatitude;  %%Y=Latitude in Km

     CheckIn=inpolygon(nodesXI,nodesYI,StateX,StateY); %1=inside,0=outside
     if CheckIn==1

         init(2,I)=Asymp(K)*init(1, I)/(StatesDensity(K)*StateArea(K));
         init(3,I)=Infected(K)*init(1, I)/(StatesDensity(K)*StateArea(K));
     end
    end
end



end


