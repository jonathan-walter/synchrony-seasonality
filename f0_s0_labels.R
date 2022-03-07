#quartz(height=6, width=4)
pdf(file="labels_temp.pdf", width=4, height=6)
layout(matrix(c(1, 2, 2, 3,
                4, 4, 5, 5,
                6, 6, 7, 7,
                8, 8, 9, 9), nrow=4, byrow=TRUE))
pal<-colorRampPalette(colors=c("red","white","blue"))
#par(mfrow=c(3,2), mar=c(1,.5,.5,.5), mgp=c(2.7,0.5,0), tcl=-0.3, oma=c(3,3,3,1))
par(mar=c(1,1,1,.5), mgp=c(2.7,0.5,0), tcl=-0.3, oma=c(3,3,.5,1))

plot.new()
#image(f0, prop, results.f, zlim=c(-1,1), col=pal(50),
#      xlab="", ylab="", yaxt="n", cex=1.25)
plot.new()
plot.new()

image(f0, prop, results.f, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", yaxt="n", cex=1.25)
contour(f0, prop, results.f, add=T)
text(0.35,.08,"g)", xpd=NA)
text(1.5, .08, expression(paste('cor(',italic(epsilon['bi']),',', italic(epsilon['bj']),')')), 
     cex=.75, xpd=NA)


image(f0, prop, results.f, zlim=c(-1,1), col=pal(50),
      xlab="", ylab="", yaxt="n", cex=1.25)
contour(f0, prop, results.f, add=T)
text(0.35,.08,"g)", xpd=NA)
text(1.5, .08, expression(paste('cor(',italic(epsilon['wi']),',', italic(epsilon['wj']),')')), 
     cex=.75, xpd=NA)
dev.off()
