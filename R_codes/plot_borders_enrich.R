#  script to plot 5 subplots

par(mfrow=c(1,5));

plot(v1/nb_borders,type="o",pch=20,lwd=5.,
     ylab="Proportion",xlab="Distance from CID booundary",col="darkblue"
     , xaxt='n',main=paste("10% most expressed genes\n","p-value=",round(p1,3)))
abline(v = 21, col = "red",lwd=5., lty = 3)
axis(1, at=1, labels= "-100kb" )
axis(1, at=11, labels= "-50kb" )
axis(1, at=21, labels= "0" )
axis(1, at=31, labels= "50kb")
axis(1, at=41, labels= "100kb")

plot(v2/nb_borders,type="o",pch=20,lwd=5.,
     ylab="Proportion",xlab="Distance from CID booundary",col="darkblue"
     , xaxt='n',main=paste("tsEPOD from Vera et al.\n","p-value=",round(p2,3)))
abline(v = 21, col = "red",lwd=5., lty = 3)
axis(1, at=1, labels= "-100kb" )
axis(1, at=11, labels= "-50kb" )
axis(1, at=21, labels= "0" )
axis(1, at=31, labels= "50kb")
axis(1, at=41, labels= "100kb")

plot(v3/nb_borders,type="o",pch=20,lwd=5.,
     ylab="Proportion",xlab="Distance from CID booundary",col="darkblue"
     , xaxt='n',main=paste("SRP genes\n","p-value=",round(p3,4)))
abline(v = 21, col = "red",lwd=5., lty = 3)
axis(1, at=1, labels= "-100kb" )
axis(1, at=11, labels= "-50kb" )
axis(1, at=21, labels= "0" )
axis(1, at=31, labels= "50kb")
axis(1, at=41, labels= "100kb")

plot(v4/nb_borders,type="o",pch=20,lwd=5.,
     ylab="Proportion",xlab="Distance from CID booundary",col="darkblue"
     , xaxt='n',main=paste("10% less expressed genes\n","p-value=",round(p4,3)))
abline(v = 21, col = "red",lwd=5., lty = 3)
axis(1, at=1, labels= "-100kb" )
axis(1, at=11, labels= "-50kb" )
axis(1, at=21, labels= "0" )
axis(1, at=31, labels= "50kb")
axis(1, at=41, labels= "100kb")

plot(v5/nb_borders,type="o",pch=20,lwd=5.,
     ylab="Proportion",xlab="Distance from CID booundary",col="darkblue"
     , xaxt='n',main=paste("SecB genes\n","p-value=",round(p5,3)))
abline(v = 21, col = "red",lwd=5., lty = 3)
axis(1, at=1, labels= "-100kb" )
axis(1, at=11, labels= "-50kb" )
axis(1, at=21, labels= "0" )
axis(1, at=31, labels= "50kb")
axis(1, at=41, labels= "100kb")
