resrange <- seq(500,5500,500);
reselw2<-log10(as.matrix(read.table("reselw2.txt", header=TRUE))) - log10(resrange);
resslicot2<-log10(as.matrix(read.table("resslicot2.txt", header=TRUE))) - log10(resrange);

ylim <- c(min(as.vector(resslicot2[,4])) -2, max(as.vector(reselw2[,6])));

postscript("testcomp.eps", horizontal=FALSE, onefile=FALSE, paper="special", width=6, height=6);
plot(resrange, resslicot2[,4], type="b", col="purple", pch=1, xlab="n", ylim=ylim, ylab="log(t/n)", ylog=TRUE);
lines(resrange, resslicot2[,5], type="b", col="blue", pch=2);
lines(resrange, resslicot2[,6], type="b", col="black", pch=3);
lines(resrange, reselw2[,4], type="b", col="dark green", pch=4);
lines(resrange, reselw2[,5], type="b", col="red", pch=5);
lines(resrange, reselw2[,6], type="b", col="brown", pch=6);

legendpos <- "bottomright";
legend("bottomright", c("function, bl.", "gradient, bl.", "p.jacobian, bl.", "function, el.", "gradient, el.", "p.jacobian, el."), col = c("purple", "blue", "black", "dark green", "red", "brown"),  pch=1:6);
dev.off();

