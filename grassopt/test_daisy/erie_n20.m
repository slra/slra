% Data of a simulation of the western basin of Lake Erie

load erie.dat
U=erie(:,1:20);
Y=erie(:,21:28);
U_erie    =U(:,1:5);
U_erie_n10=U(:,6:10);
U_erie_n20=U(:,11:15);	
U_erie_n30=U(:,16:20);
Y_erie    =Y(:,1:2);
Y_erie_n10=Y(:,3:4);
Y_erie_n20=Y(:,5:6);
Y_erie_n30=Y(:,7:8);
u = U_erie_n20;
y = Y_erie_n20;
clear erie U_erie U_erie_n10 U_erie_n20 U_erie_n30 Y_erie Y_erie_n10 Y_erie_n20 Y_erie_n30