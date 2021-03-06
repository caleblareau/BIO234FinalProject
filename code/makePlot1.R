x <- c(1,2,5,10,100)
y1 <- c(326.35,180.49,70.25,32.66,3.33)
y2 <- c(658.23, 905.08,761.23,705.93,663.23)
y3 <- c(146.6,87,40.5,21,2.1)
jpeg("/Users/Ploenzke/Documents/Harvard/Data Structures/Final/paper/plot1.jpg", width=600,height=500)
par(mar=c(5,4,4,5)+.1)
plot(lowess(x,y1,f=.1),type="l",col=2,main="Algorithm Efficiency", xlab="Thinning",ylab="Time (s)",ylim=c(0,1000))
lines(lowess(x,y2,f=.1),type="l",col=3)
par(new=TRUE)
plot(lowess(x,y3,f=.1),type="l",col=4,xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("Size (MB)",side=4,line=3)
legend(60,80,col=c(2,3,4),lty=1,legend=c("Load Time","Compress Time","File Size"))
dev.off()