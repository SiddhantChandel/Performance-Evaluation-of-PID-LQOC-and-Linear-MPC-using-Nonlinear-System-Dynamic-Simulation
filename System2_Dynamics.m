%Assignment: 3, Name: Anshumali Jaiswal, Roll nubmer: 213020033 

function fofX= System2_Dynamics(t, X)

global sys

fofX = zeros(3,1) ;

k1=sys.alfa(1)*exp(-5000/X(2));
k2=sys.alfa(2)*exp(-7500/X(2));

fofX(1)=(-sys.alfa(3)*sys.Uk(1)*X(1)/X(3))+k1*(sys.Dk-X(1))-k2*X(1);

fofX(2)=(sys.alfa(3)*sys.Uk(1)*(sys.Uk(2)-X(2))/X(3))+ sys.alfa(4)*(k1*(sys.Dk-X(1))-k2*X(1));

fofX(3)=sys.alfa(3)*sys.Uk(1)- sys.alfa(5)*(X(3)^0.5);

