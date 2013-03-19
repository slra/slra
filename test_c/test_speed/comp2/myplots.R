resrange <- seq(50,500,50);


resslicot2<-as.matrix(read.table("resslicot.txt", header=TRUE));

resslicot2[,6] <- resslicot2[,6];
xlim <- c(min(as.vector(resslicot2[,4])) , max(as.vector(resslicot2[,6])));
ylim <- c(43, 560);


postscript("testcomp2.eps", horizontal=FALSE, onefile=FALSE, paper="special", width=7.8, height=2.8)
par(mai=c(0.6,0.6,0.4,0.4), mgp=c(1.8,0.6,0));

#abline( v = 0.1 * 2^(1:5), lty = 3, col = colors()[ 440 ] )


plot(resslicot2[,4], resrange, type="b", col="purple", pch=1, xlab="t", xlim=xlim, ylim=ylim, ylab="m", log="xy", asp = 1);
abline( v = 0.002 * 10^(0.5*(0:9)), lty = 3, col = "black" )
abline( h = c(50,50 * sqrt(10),500), lty = 3, col = "black" )

lines(resslicot2[,5], resrange, type="b", col="blue", pch=2);
lines(resslicot2[,6], resrange, type="b", col="black", pch=3);

legendpos <- "bottomright";
legend(legendpos, c("function, bl.", "gradient, bl.", "p.jacobian, bl. "), col = c("purple", "blue", "black"),  pch=1:3, bg="white");
dev.off();

