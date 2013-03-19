resrange <- seq(500,5500,500);
reselw2<-as.matrix(read.table("reselw2.txt", header=TRUE));
resslicot2<-as.matrix(read.table("resslicot2.txt", header=TRUE));

xlim <- c(min(as.vector(resslicot2[,4])) /10, max(as.vector(reselw2[,6])))

postscript("testcomp1.eps", horizontal=FALSE, onefile=FALSE, paper="special", width=8.5, height=2.4);
par(mai=c(0.6,0.6,0.4,0.4), mgp=c(1.8,0.6,0), yaxs="i");
plot(resslicot2[,4], resrange, type="b", col="purple", pch=1, xlab="t", ylab="n", log="xy", xlim = xlim);
abline( v = 1e-04 * 10^(0:4), lty = 3, col = "black" )
abline( h = c(500,5000), lty = 3, col = "black" )

lines(resslicot2[,5], resrange, type="b", col="blue", pch=2);
lines(resslicot2[,6], resrange, type="b", col="black", pch=3);
lines(reselw2[,4], resrange, type="b", col="dark green", pch=4);
lines(reselw2[,5], resrange, type="b", col="red", pch=5);
lines(reselw2[,6], resrange, type="b", col="brown", pch=6);

legendpos <- "bottomright";
legend("topleft", c("function, bl.", "gradient, bl.", "p.jacobian, bl.", "function, el.", "gradient, el.", "p.jacobian, el."), col = c("purple", "blue", "black", "dark green", "red", "brown"),  pch=1:6, bg="white");
dev.off();

