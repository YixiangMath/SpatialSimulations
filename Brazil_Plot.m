%%%Created on 7/25/2020
%%%This file plot the solution of Brazil coronavirus model
%%%The solution was stored in sol
s=sol.NodalSolution;

lambda=1/7;
f=0.2;

%%%Plot initial susceptible data
h=figure;
%pdeplot(pdem,'XYData',s(:,1,1), 'ZData',s(:,1,1))
pdeplot(pdem,'XYData',real(log10(s(:,1,1))),'ZData',s(:,1,1))
view(0, 68)
saveas(h, 'Susceptible', 'epsc');
%title('S_0')


% h=figure;
% %pdeplot(pdem,'XYData',s(:,1,1), 'ZData',s(:,1,1))
% 
% tau1=0.1*10^-3;
% tau2=0.01*10^-2;
% lambda=1/7;
% gamma=1/7;
% %R0=abs((s(:,1,1)*tau1+sqrt(s(:,1,1)/gamma).*sqrt(s(:,1,1)*gamma*tau1^2+4*lambda^2*tau2))/(2*lambda));
% R0=tau1/lambda+tau2/gamma;
% pdeplot(pdem,'XYData',R0, 'ZData', R0)
% saveas(h, 'R0', 'epsc');

%%%Plot infected population data
T=[1 10];
h=figure;
Len=length(T);
for I=1:Len
subplot(1, 2, I)
%pdeplot(pdem,'XYData',log10(s(:,3,T(I))), 'ZData',s(:,3,1))
%pdeplot(pdem,'XYData',real(log10(s(:,3,T(I)))).*sign(real(log10(s(:,1,1)))), 'ZData',s(:,3,T(I)))
pdeplot(pdem,'XYData',real(log10(s(:,3,T(I)))), 'ZData',s(:,3,T(I)))
view(0, 90)
title(['Distribution of infected cases at time t=', num2str(T(I))])
end
saveas(h, 'Infected', 'epsc');


%%%%p = the coordinates of the nodes
p=pdem.Mesh.Nodes;
nodesX=p(1,:);
nodesY=p(2,:);

%%%triangles
t=pdem.Mesh.Elements;

%%%%%%%%%%%%%%
totalS0=0;
Nt=length(t);

ss=zeros(Nt,2);
%ss=u0fun(nodesX', nodesY', phi)';
ss=s(:,:,1);

%%Total population of Brazil
for I=1:Nt
    Pa=t(1, I);
    Pb=t(2, I);
    Pc=t(3, I);
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
    totalS0=totalS0+area*(ss(Pa, 1)+ss(Pb, 1)+ss(Pc, 1))/3;
end
fprintf('Total population of Brazil is %d\n',round(totalS0))

%%%%% 
%%Total population of initial infected population

%%%%

totalI0=0;
for I=1:Nt
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
    totalI0=totalI0+area*(ss(Pa, 2)+ss(Pb, 2)+ss(Pc, 2))/3;
end
fprintf('Initial infected population is %d\n',round(totalI0))


TimeLength=length(tlist);
TotalInfected=zeros(TimeLength,1);


for J=1:TimeLength
    total=0;
for I=1:Nt
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
   TotalInfected(J)=total;
end




for I=1: TimeLength
    AccumulatedTotal(I)=sum(TotalInfected(1:I));
end



h=figure;
plot(tlist, AccumulatedTotal)
xlabel('Days')
ylabel('Total accumulated reported cases')


hold on
clear T
clear TotalStatesData

for I=1:27
   s=num2str(I);
   filename=s+".xlsx";
   T=readtable(filename, 'Range', 'B3:B32', 'ReadVariableNames',false);
   T=table2array(T);
   TotalStatesData(:,I)=T;
end
TotalData=sum(TotalStatesData,2);   

plot((0: (length(TotalData)-1)), TotalData, 'o');
saveas(h, 'Reported', 'epsc');
hold off







