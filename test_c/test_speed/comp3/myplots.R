resrange <- seq(20,120,10);


resslicot2<-as.matrix(read.table("reselw.txt", header=TRUE));
resslicot2<-(resslicot2/(resrange));
resslicot2[,6] <-resslicot2[,6] / 50;
ylim <- c(min(as.vector(resslicot2[,4])) , max(as.vector(resslicot2[,6])));


postscript("testcomp3.eps", horizontal=FALSE, onefile=FALSE, paper="special", width=6, height=6)
plot(resrange, resslicot2[,4], type="b", col="purple", pch=1, xlab="m", ylim=ylim, ylab="t/m", );
lines(resrange, resslicot2[,5], type="b", col="blue", pch=2);
lines(resrange, resslicot2[,6], type="b", col="black", pch=3);

legendpos <- "topleft";
legend(legendpos, c("function, el.", "gradient, el.", "p.jacobian, el. / 50"), col = c("purple", "blue", "black"),  pch=1:3);
dev.off();

