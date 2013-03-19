resrange <- seq(10,100,10);


resslicot2<-as.matrix(read.table("reselw.txt", header=TRUE));
resslicot2[,6] <-resslicot2[,6];
xlim <- c(min(as.vector(resslicot2[,4])) , max(as.vector(resslicot2[,6])));
ylim <- c(7.5, 135);


postscript("testcomp3.eps", horizontal=FALSE, onefile=FALSE, paper="special", width=7.8, height=3)
par(mai=c(0.6,0.6,0.4,0.4), mgp=c(1.8,0.6,0), yaxs="i");


plot(resslicot2[,4], resrange, type="b", col="purple", pch=1, xlab="t", xlim=xlim, ylim=ylim, ylab="m", log="xy");
abline( v = 0.005 * 10^(0.5*(0:8)), lty = 3, col = "black" )
abline( h = c(10,10*sqrt(10),100), lty = 3, col = "black" )
lines(resslicot2[,5], resrange, type="b", col="blue", pch=2);
lines(resslicot2[,6], resrange, type="b", col="black", pch=3);

legendpos <- "bottomright";
legend(legendpos, c("function, el.", "gradient, el.", "p.jacobian, el."), col = c("purple", "blue", "black"),  pch=1:3, bg="white");
dev.off();

