
%%%This file is used to graph the infected population in each state

%%Read the boundary data for Puerto Rico
f=0.2;
lambda=1/7;
clear TotalStatesData

BRShape1=shaperead('bra_adm_ibge_2020_shp/bra_admbnda_adm1_ibge_2020.shp');

%%%%%%%%
%%%%%The following codes delete NaN values in the shp file
%%%%Use this to find the NaN values [row, col] = find(isnan(BRShape1(1).X));
BRShape1(1).X=BRShape1(1).X(1, 23:2784);BRShape1(1).Y=BRShape1(1).Y(1, 23:2784);
BRShape1(11).X=BRShape1(11).X(1, 20:7129);BRShape1(11).Y=BRShape1(11).Y(1, 20:7129);
BRShape1(9).X=BRShape1(9).X(1, 26:7582);BRShape1(9).Y=BRShape1(9).Y(1, 26:7582);
BRShape1(13).X=BRShape1(13).X(1, 13:9807);BRShape1(13).Y=BRShape1(13).Y(1, 13:9807);
BRShape1(5).X=BRShape1(5).X(1, 32:9676);BRShape1(5).Y=BRShape1(5).Y(1, 32:9676);


for I=1:27
   s=num2str(I);
   filename=s+".xlsx";
  % filename="/Users/ywu/MATLAB_Drive/BrazilCoronavirus/BrazilCovidData/"+s+".xlsx";
   T=readtable(filename, 'Range', 'B3:B32', 'ReadVariableNames',false);
   T=table2array(T);
   TotalStatesData(:,I)=T;
end
%TotalStatesData=table2array(TotalStatesData);
%TotalData=sum(TotalStatesData,2);


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


%Region=[1 7 19 20 21 24]; 

Region=1:27; 

for K=1:length(Region)

%%Obtained Arecibo boundary data
X=BRShape1(Region(K)).X;
Y=BRShape1(Region(K)).Y;

%%Convert Arecibo Latitudes/Logitudes data into XY coordinates 

X1= (X+50)*DeltaLongitude; 
Y1=(Y+12)*DeltaLatitude;  

%AreciboTriangles=zeros(Nt, 1);
t1=t(1,:); %First point of triangle
t2=t(2,:); % Second point of triangle
t3=t(3,:);
CenterX=(nodesX(t1)+nodesX(t2)+nodesX(t3))/3; %x coordinates of the centers 
CenterY=(nodesY(t1)+nodesY(t2)+nodesY(t3))/3;
AreciboTriangles=inpolygon(CenterX,CenterY,X1,Y1); %1=inside,0=outside


%%%%%%%%%%%%%
%%Graph the total infected population in Arecibo

TimeLength=length(tlist);
TotalInfectedA=zeros(TimeLength,1);

for J=1:TimeLength
    total=0;
for I=1:Nt
    if AreciboTriangles(I)==1
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
    total=total+f*lambda*area*(s(Pa, 2, J)+s(Pb, 2, J)+s(Pc, 2, J))/3;
    end
end
   TotalInfectedA(J)=total;
end

for I=1: TimeLength
    AccumulatedTotalState(I)=sum(TotalInfectedA(1:I));
    %sum(TotalData(((I-1)*7+1):((I-1)*7+7)));
end


h=figure
plot(tlist, AccumulatedTotalState)
hold on
xlabel('Days')
ylabel('Accumulated reported cases')
plot((1:length(TotalStatesData(:, K))), TotalStatesData(:, K), 'o');
%%title(['State ', num2str(K)])
%name=append('state', num2str(K));
%saveas(h, name, 'epsc');
end

