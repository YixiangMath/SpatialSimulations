
%%This file graphs the density of each state 

f=0.2;
BRShape1=shaperead('bra_adm_ibge_2020_shp/bra_admbnda_adm1_ibge_2020.shp');
BRShape1(1).X=BRShape1(1).X(1, 23:2784);BRShape1(1).Y=BRShape1(1).Y(1, 23:2784);
BRShape1(11).X=BRShape1(11).X(1, 20:7129);BRShape1(11).Y=BRShape1(11).Y(1, 20:7129);
BRShape1(9).X=BRShape1(9).X(1, 26:7582);BRShape1(9).Y=BRShape1(9).Y(1, 26:7582);
BRShape1(13).X=BRShape1(13).X(1, 13:9807);BRShape1(13).Y=BRShape1(13).Y(1, 13:9807);
BRShape1(5).X=BRShape1(5).X(1, 32:9676);BRShape1(5).Y=BRShape1(5).Y(1, 32:9676);
for D=1:27
    len=length(BRShape1(D).X);
    BRShape1(D).X=BRShape1(D).X(1, 1:(len-1));
    BRShape1(D).Y=BRShape1(D).Y(1, 1:(len-1));
end




for I=1:27
   s=num2str(I);
   %filename="BrazilCoronavirus/"+s+".xlsx";
   filename=s+".xlsx";
  % filename="/Users/ywu/MATLAB_Drive/BrazilCoronavirus/BrazilCovidData/"+s+".xlsx";
   T=readtable(filename, 'Range', 'B3:B32', 'ReadVariableNames',false);
   T=table2array(T);
   TotalStatesData(:,I)=T;
end




s=sol.NodalSolution;

%%%%p = the coordinates of the nodes
p=pdem.Mesh.Nodes;
nodesX=p(1,:);
nodesY=p(2,:);

%%%triangles 
t=pdem.Mesh.Elements;

%%%%%%%%%%%%%%
%%Compute the total population in Arecibo District
%%Find the triangles contained in Arecibo
Nt=length(t);

t1=t(1,:); %First point of triangle
t2=t(2,:); % Second point of triangle
t3=t(3,:);
CenterX=(nodesX(t1)+nodesX(t2)+nodesX(t3))/3; %x coordinates of the centers 
CenterY=(nodesY(t1)+nodesY(t2)+nodesY(t3))/3;

TIME=[1  2   4    6    10    16    22]


StateArea=[164123,27779,142829,1559159,564733,148920,5780,46096,340112,331937,903366,357146,586522,1247955,56470,199308,98148,251578,43780,52811,281730,237591,224301,  95736,248223,21915,277721];

%%%%%%State Population Density
fileID=fopen('StatesDensity.txt');
C=textscan(fileID,'%f');
StatesDensity=cell2mat(C);
fclose(fileID);
StatePopulation=StateArea'.*StatesDensity;

for K=1:length(TIME)
    
T=TIME(K); 
%DistrictPopulation=zeros(27,1);%%District D total population
DistrictInfected=zeros(27,1);%District D total infected population



for D=1:27

%%Obtained District boundary data
DistrictX=BRShape1(D).X;
DistrictY=BRShape1(D).Y;

%%Convert District Latitudes/Logitudes data into XY coordinates 
e=sqrt(0.00669437999014);
a=6378137.0;

phi=mean(Latitude)/180*pi;    %%%%Latitude is computed in Brazil_Solve.m

DeltaLatitude=pi*a*(1-e^2)./(180*( 1-e^2 * (sin(phi)).^2   ).^1.5     )/1000; %% 1 Latitude= 1 DeltaLatitude Km
DeltaLongitude= pi*a*cos(phi)./  (180*( 1-e^2 * (sin(phi)).^2   ).^0.5     )/1000; 


DistrictX1= (DistrictX+50)*DeltaLongitude; 
DistrictY1=(DistrictY+12)*DeltaLatitude;  

%DistrictTriangles=zeros(Nt, 1);
DistrictTriangles=inpolygon(CenterX,CenterY,DistrictX1,DistrictY1); %1=inside,0=outside


%%%%%%%%%%%%%
%%Compute the total infected population and total population in Arecibo

for I=1:Nt
    if DistrictTriangles(I)==1
    Pa=t(1, I); Pb=t(2, I); Pc=t(3, I);
    Edge1X=nodesX(Pa);
    Edge1Y=nodesY(Pa);
    Edge2X=nodesX(Pb);
    Edge2Y=nodesY(Pb);
    Edge3X=nodesX(Pc);
    Edge3Y=nodesY(Pc);
    L1=sqrt( (Edge1X-Edge2X)^2+ (Edge1Y-Edge2Y)^2);
    L2=sqrt( (Edge2X-Edge3X)^2+ (Edge2Y-Edge3Y)^2);
    L3=sqrt( (Edge3X-Edge1X)^2+ (Edge3Y-Edge1Y)^2);
    rho=(L1+L2+L3)/2;
    area=sqrt(rho*(rho-L1)*(rho-L2)*(rho-L3));
    DistrictInfected(D)=DistrictInfected(D)+area*(s(Pa, 3, T)+s(Pb, 3, T)+s(Pc, 3, T))/3;
  %  DistrictPopulation(D)=DistrictPopulation(D)+area*(s(Pa, 1, 1)+s(Pb, 1, 1)+s(Pc, 1, 1))/3;
    end
end
end





RelativeDensity=f*DistrictInfected./StatePopulation;
RelativeD=RelativeDensity*1000000;


%RelativeD=DistrictInfected;
Max=max(RelativeD);
%Scale=Max/6;

color=zeros(27,3);
for D=1:27
    
    if  RelativeD(D)<=Max/8
        color(D,:)= [0.9792    0.9375    0.8984];
    end
    if RelativeD(D)>Max/8 && RelativeD(D)<=2*Max/8
        color(D,:)=[0.9083    0.8672    0.6992];
    end
    if RelativeD(D)>2*Max/8 && RelativeD(D)<=3*Max/8
        color(D,:)=[0.85253    0.8072    0.6592];
    end
    if RelativeD(D)>3*Max/8 && RelativeD(D)<=4*Max/8
        color(D,:)=[0.8283    0.7572    0.6092];
    end 
    if RelativeD(D)>4*Max/8 && RelativeD(D)<=5*Max/8
        color(D,:)=[0.8008    0.3594    0.3594];
    end 
    if RelativeD(D)>5*Max/8 && RelativeD(D)<=6*Max/8
        color(D,:)=[0.7008    0.2594    0.2594];
    end 
    if RelativeD(D)>6*Max/8 && RelativeD(D)<=7*Max/8
        color(D,:)=[0.6008    0.1594    0.1594];
    end 
    if RelativeD(D)>7*Max/8
        color(D,:)=[ 0.55         0         0];
    end
    
%     if  RelativeD(D)<=Max/6
%         color(D,:)= [0.9792    0.9375    0.8984];
%     end
%     if RelativeD(D)>Max/6 && RelativeD(D)<=2*Max/6
%         color(D,:)=[0.9583    0.8672    0.6992];
%     end
%     if RelativeD(D)>2*Max/6 && RelativeD(D)<=3*Max/6
%         color(D,:)=[0.8008    0.3594    0.3594];
%     end
%     if RelativeD(D)>3*Max/6 && RelativeD(D)<=4*Max/6
%         color(D,:)=[0.7008    0.2594    0.2594];
%     end 
%     if RelativeD(D)>4*Max/6 && RelativeD(D)<=5*Max/6
%         color(D,:)=[0.6008    0.1594    0.1594];
%     end 
%     if RelativeD(D)>5*Max/6
%         color(D,:)=[ 0.5430         0         0];
%     end
end

h=figure
for D=1:27
    symspec = makesymbolspec('Polygon', {'Default','FaceColor', color(D,:)});
      mapshow(BRShape1(D),'SymbolSpec',symspec)
end 

% map = [0.9792    0.9375    0.8984
%     0.9583    0.8672    0.6992
%     0.8008    0.3594    0.3594
%     0.7008    0.2594    0.2594
%     0.6008    0.1594    0.1594
%     0.55         0         0];
map = [0.9792    0.9375    0.8984
    0.9083    0.8672    0.6992
    0.85253    0.7072    0.6092
    0.8283    0.5572    0.5092
    0.8008    0.3594    0.3594
    0.7008    0.2594    0.2594
    0.6008    0.1594    0.1594
    0.55         0         0];
colormap(map)
colorbar('Ticks',[0, 1/8,2/8,3/8,4/8,5/8, 6/8,7/8,1],...
         'TickLabels',{'0',num2str(round(Max/8)),num2str(round(2*Max/8)),num2str(round(3*Max/8)),num2str(round(4*Max/8)), num2str(round(5*Max/8)), num2str(round(6*Max/8)), num2str(round(7*Max/8)), num2str(round(8*Max/8))});

% colorbar('Ticks',[0, 1/8,2/8,3/8,4/8,5/8, 6/8,7/8,1],...
%          'TickLabels',{'0',num2str(round(Max/6)),num2str(round(2*Max/6)),num2str(round(3*Max/6)),num2str(round(4*Max/6)), num2str(round(5*Max/6)), num2str(round(6*Max/6))});
title(['Day ', num2str(T)])

name=append('Day', num2str(T));
saveas(h, name, 'epsc');




%%%%%%%%Real
RelativeDensity1=TotalStatesData(T+1,:)'./StatePopulation;
RelativeD1=RelativeDensity1*1000000;
%Max1=max(RelativeD1);
Max1=Max;

color1=zeros(27,3);
for D=1:27
    
    if  RelativeD1(D)<=Max1/8
        color1(D,:)= [0.9792    0.9375    0.8984];
    end
    if RelativeD1(D)>Max1/8 && RelativeD1(D)<=2*Max1/8
        color1(D,:)=[0.9083    0.8672    0.6992];
    end
    if RelativeD1(D)>2*Max1/8 && RelativeD1(D)<=3*Max1/8
        color1(D,:)=[0.85253    0.8072    0.6592];
    end
    if RelativeD1(D)>3*Max1/8 && RelativeD1(D)<=4*Max1/8
        color1(D,:)=[0.8283    0.7572    0.6092];
    end 
    if RelativeD1(D)>4*Max1/8 && RelativeD1(D)<=5*Max1/8
        color1(D,:)=[0.8008    0.3594    0.3594];
    end 
    if RelativeD1(D)>5*Max1/8 && RelativeD1(D)<=6*Max1/8
        color1(D,:)=[0.7008    0.2594    0.2594];
    end 
    if RelativeD1(D)>6*Max1/8 && RelativeD1(D)<=7*Max1/8
        color1(D,:)=[0.6008    0.1594    0.1594];
    end 
    if RelativeD1(D)>7*Max1/8
        color1(D,:)=[ 0.55         0         0];
    end  
end

h1=figure
for D=1:27
    symspec = makesymbolspec('Polygon', {'Default','FaceColor', color1(D,:)});
      mapshow(BRShape1(D),'SymbolSpec',symspec)
end 

colormap(map)
colorbar('Ticks',[0, 1/8,2/8,3/8,4/8,5/8, 6/8,7/8,1],...
         'TickLabels',{'0',num2str(round(Max1/8)),num2str(round(2*Max1/8)),num2str(round(3*Max1/8)),num2str(round(4*Max1/8)), num2str(round(5*Max1/8)), num2str(round(6*Max1/8)), num2str(round(7*Max1/8)), append(num2str(round(8*Max1/8)), '+')});
title(['Day ', num2str(T)])

name=append('AcDay', num2str(T));
saveas(h1, name, 'epsc');



end


