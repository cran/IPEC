biasIPEC <- function (expr, theta, x, y, tol = 1e-16, method = "Richardson", 
    method.args = list(eps = 1e-04, d = 0.11, zero.tol = sqrt(.Machine$double.eps/7e-07), 
    r = 6, v = 2, show.details = FALSE), side = NULL){

    x <- rbind(x)
    y <- as.vector(y)
    if (min(dim(x))[1] == 1) 
        x <- cbind(x)
    if (nrow(x) != length(y)) 
        x <- t(x)
    n <- nrow(x)
    p <- length(theta)
    yhat <- expr(theta, x)
    sigma <- sqrt(sum((y - yhat)^2)/(n - p))
    sp <- sigma * sqrt(p)
    v1 <- matrix(NA, nrow = n, ncol = p)
    v2 <- array(NA, c(p, p, n))
    for (k in 1L:n) {
        resu <- derivIPEC(expr, theta, x[k, ], method = method, 
            method.args = method.args, side = side)
        v1[k, ]   <- resu$Jacobian
        v2[, , k] <- resu$Hessian
    }
    a <- array(NA, c(n, p, p))
    for (k in 1L:p) {
        temp.mat <- matrix(NA, nrow = 1, ncol = p)
        for (q in 1L:n) temp.mat <- rbind(temp.mat, v2[k, , q])
        temp.mat <- temp.mat[-1, ]
        a[, , k] <- temp.mat
    }
    D <- v1
    for (j in 1L:p) {
        D <- cbind(D, a[, 1L:p, j])
    }
    qrd <- qr(D, tol)
    Q   <- qr.Q(qrd)
    rnk <- qrd$rank
    if (rnk <= p) {
        warning("regression apparently linear")
    }
    Q1 <- Q[, 1L:rnk]
    C  <- array(0, c(rnk, p, p))
    for (j in 1L:p) {
        C[, , j] <- crossprod(Q1, a[, , j])
    }
    C     <- aperm(C, c(2, 3, 1))
    r11i  <- ginv(qr.R(qrd)[1L:p, 1L:p], tol=tol)
    sumaT <- c()
    for (j in 1L:p) {
        C[, , j] <- crossprod(r11i, C[, , j]) %*% r11i * sp 
        sumaT    <- c( sumaT, sum(diag(C[, , j])) ) 
    }
    bias <- -1/(2*p) * r11i %*% t(t(sumaT)) * sp
    bias <- array( bias )
    percent.bias <- bias/theta * 100
    list(bias=bias, percent.bias=percent.bias)

}

