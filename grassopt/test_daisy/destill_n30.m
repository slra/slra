% ethane-ethylene destillation column

load destill.dat
U=destill(:,1:20);
Y=destill(:,21:32);
clear destill

U_dest     = U(:,1:5);
U_dest_n10 = U(:,6:10);
U_dest_n20 = U(:,11:15);	
U_dest_n30 = U(:,16:20);
Y_dest     = Y(:,1:3);
Y_dest_n10 = Y(:,4:6);
Y_dest_n20 = Y(:,7:9);
Y_dest_n30 = Y(:,10:12);

u = U_dest_n30;
y = Y_dest_n30;

clear U_dest U_dest_n10 U_dest_n20 U_dest_n30 Y_dest Y_dest_n10 Y_dest_n20 Y_dest_n30
