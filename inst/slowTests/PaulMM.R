library(IPEC)

# The velocity of the reaction (counts/min^2) under different substrate concentrations 
#   in parts per million (ppm) (Page 269 of Bates and Watts 1988)

x1 <- c(0.02, 0.02, 0.06, 0.06, 0.11, 0.11, 0.22, 0.22, 0.56, 0.56, 1.10, 1.10)
y1 <- c(76, 47, 97, 107, 123, 139, 159, 152, 191, 201, 207, 200)

# Define the Michaelis-Menten model
MM <- function(theta, x){
    theta[1]*x / ( theta[2] + x )    
}




res0 <- fitIPEC(MM, x=x1, y=y1, ini.val=c(200, 0.05), 
              xlim=c( 0, 1.5 ), ylim=c(0, 250), fig.opt=FALSE)

par1 <- res0$par
print( par1, digits=16 )
# R-3.4.0 (64bits) on Windows 7 (64bits)
# [1] 212.6849086499876194   0.0641242102763406

# The tolerance in these tests may need to be larger for some calculations, 
# you can get slightly different values on different platforms.

if (.Machine$double.eps*1e6 < max(abs(par1/c(1e3, 1) - c(0.2126849086499876194, 0.0641242102763406))))
   stop("MM fitIPEC par comparison test failed.")



   
res1 <- derivIPEC(MM, theta=par1, z=x1[1], method="Richardson",
            method.args=list(eps=1e-4, d=0.11,
	    zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2))

		
print( res1$Jacobian, digits=16 )
# R-3.4.0 (64bits) on Windows 7 (64bits)
# [1]    0.2377436879859202 -601.0694709085224758

if (.Machine$double.eps*1e6 < max(abs(res1$Jacobian/c(1, 1e3) - c(0.2377436879859202, -0.6010694709085224758)
        ))) stop("MM derivIPEC Jacobian comparison test failed.")

		
print( res1$Hessian, digits=16 )
# R-3.4.0 (64bits) on Windows 7 (64bits)
#                       [,1]                  [,2]
#[1,] -2.090700943985578e-14    -2.826103058823163 
#[2,] -2.826103058823163e+00 14290.047274752103476

if (.Machine$double.eps*1e6 < max(abs(res1$Hessian/matrix(c(1, 1e1, 1e1, 1e5), nrow=2, byrow=TRUE) 
        - matrix(c(0, -0.2826103058823163, -0.2826103058823163, 
        0.14290047274752103476), nrow=2, byrow=TRUE)
        ))) stop("MM derivIPEC Hessian comparison test failed.")




		
# To calculate curvatures
res2 <- curvIPEC( MM, theta=par1, x=x1, y=y1, alpha=0.05, method="Richardson",
                  method.args=list(eps=1e-4, d=0.11,
	              zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2)) 

print( res2, digits=16 )
# R-3.4.0 (64bits) on Windows 7 (64bits)
# $rms.ic      0.04542212409811858
# $rms.pec     0.1047135056832769
# $critical.c  0.493694983503531

if (.Machine$double.eps*1e6 < max(abs(c(res2$rms.ic, res2$rms.pec, res2$critical.c) - 
    c(0.04542212409811858, 0.1047135056832769, 0.493694983503531)
        ))) stop("MM curvIPEC  comparison test failed.")



		
# To calculate bias
res3 <- biasIPEC(MM, theta=par1, x=x1, y=y1, tol= 1e-20)

print( res3$bias, digits=16 )
# R-3.4.0 (64bits) on Windows 7 (64bits)
# [1] 0.1900061522450134099 0.0004426350233523202

if (.Machine$double.eps*1e6 < max(abs(res3$bias*c(1, 1e3) - 
    c(0.1900061522450134099, 0.4426350233523202)
        ))) stop("MM biasIPEC  comparison test failed.")


		
		
RNGkind( kind = "Mersenne-Twister",  normal.kind = "Inversion")

res4 <- bootIPEC( MM, x=x1, y=y1, ini.val=par1, target.fun = "RSS", 
                  control=list(reltol=1e-20, maxit=40000), nboot=2000,
		          alpha=0.05, fig.opt=FALSE, seed=123, prog.opt=FALSE )

		   
print(sum(res4$M), digits=16)
# R-3.4.0 (64bits) on Windows 7 (64bits)
#[1] 2466259.66117283
		
if (.Machine$double.eps*1e7 < abs(sum(res4$M)/1e7 - 0.246625966117283)
        ) stop("MM bootIPEC $M comparison test failed.")


print( res4$perc.ci.mat, digits=16 )
# R-3.4.0 (64bits) on Windows 7 (64bits)
#                       Estimate                   SD                Median
# theta[1] 212.68374352538936023 6.573979247635146628 213.14895542046480159
# theta[2]   0.06412128193294594 0.009360044929729298   0.06429521592348557
#                           Mean              perc LCI              perc UCI
# theta[1] 212.30894946410506918 200.93151579045834865 221.10269747553326170
# theta[2]   0.06393910167365632   0.04796804239696822   0.07851481334250451

if (.Machine$double.eps*1e7 < max(abs(res4$perc.ci.mat/matrix(c(
    1e3, 1e1, 1e3, 1e3, 1e3, 1e3, 1, 1, 1, 1, 1, 1), 
	nrow=2, byrow=TRUE) - matrix(c(
    0.21268374352538936023, 0.6573979247635146628, 0.21314895542046480159, 
    0.21230894946410506918, 0.20093151579045834865, 0.22110269747553326170,
    0.06412128193294594, 0.009360044929729298, 0.06429521592348557, 
    0.06393910167365632, 0.04796804239696822, 0.07851481334250451), 
    nrow=2, byrow=TRUE)))) stop("MM bootIPEC $perc.ci.mat comparison test failed.")


print( res4$bca.ci.mat, digits=16 )
# R-3.4.0 (64bits) on Windows 7 (64bits)
#                       Estimate                   SD                Median
# theta[1] 212.68374352538936023 6.573979247635146628 213.14895542046480159
# theta[2]   0.06412128193294594 0.009360044929729298   0.06429521592348557
#                           Mean               bca LCI               bca UCI
# theta[1] 212.30894946410506918 198.15978346528189036 220.11160916334250714
# theta[2]   0.06393910167365632   0.04473035762648244   0.07633008517787113

if (.Machine$double.eps*1e7 < max(abs(res4$bca.ci.mat/matrix(c(
    1e3, 1e1, 1e3, 1e3, 1e3, 1e3, 1, 1, 1, 1, 1, 1), 
	nrow=2, byrow=TRUE) - matrix(c(
    0.21268374352538936023, 0.6573979247635146628, 0.21314895542046480159, 
    0.21230894946410506918, 0.19815978346528189036, 0.22011160916334250714, 
    0.06412128193294594, 0.009360044929729298, 0.06429521592348557, 
    0.06393910167365632, 0.04473035762648244, 0.07633008517787113), 
    nrow=2, byrow=TRUE)))) stop("MM bootIPEC $ca.ci.mat comparison test failed.")


print( res4$covar.mat, digits=16 )
# R-3.4.0 (64bits) on Windows 7 (64bits)
# [1,] 43.21720314833756049 5.101175298380493e-02
# [2,]  0.05101175298380493 8.761044108655115e-05

if (.Machine$double.eps*1e10 < max(abs(res4$covar.mat/matrix(c(1e2, 1, 1, 1),
    nrow=2, byrow=TRUE) - matrix(c( 
    0.4321720314833756049, 5.101175298380493e-02,
    0.05101175298380493, 8.761044108655115e-05), nrow=2, byrow=TRUE)
        ))) stop("MM bootIPEC $covar.mat comparison test failed.")


print( res4$cor.mat, digits=16 )
# R-3.4.0 (64bits) on Windows 7 (64bits)
# [1,] 1.0000000000000000 0.8290180487847936
# [2,] 0.8290180487847936 1.0000000000000000

if (.Machine$double.eps*1e10 < max(abs(res4$cor.mat - matrix(c( 
    1, 0.8290180487847936, 0.8290180487847936, 1), nrow=2, byrow=TRUE)
        ))) stop("MM bootIPEC $cor.mat comparison test failed.")



		
# To calculate skewness
res5 <- skewIPEC(MM, theta=par1, x=x1, y=y1, tol= 1e-20)

print( res5, digits=16 )
# R-3.4.0 (64bits) on Windows 7 (64bits)
# $skewness
# [1] 0.09605656534651008 0.32070003611121428		

if (.Machine$double.eps*1e6 < max(abs(res5$skewness - 
    c(0.09605656534651008, 0.32070003611121428)
        ))) stop("MM skewIPEC $skewness comparison test failed.")		
		