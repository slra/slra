resrange <- seq(100,200,10);


resslicot2<-as.matrix(read.table("resslicot.txt", header=TRUE))/(resrange);

resslicot2[,6] <- resslicot2[,6]/20;
ylim <- c(min(as.vector(resslicot2[,4])) , max(as.vector(resslicot2[,6])));


postscript("testcomp2.eps", horizontal=FALSE, onefile=FALSE, paper="special", width=6, height=6)
plot(resrange, resslicot2[,4], type="b", col="purple", pch=1, xlab="m", ylim=ylim, ylab="t/m");
lines(resrange, resslicot2[,5], type="b", col="blue", pch=2);
lines(resrange, resslicot2[,6], type="b", col="black", pch=3);

legendpos <- "topleft";
legend(legendpos, c("function, bl.", "gradient, bl.", "p.jacobian, bl. / 20"), col = c("purple", "blue", "black"),  pch=1:3);
dev.off();

