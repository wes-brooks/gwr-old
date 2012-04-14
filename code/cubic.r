rk <- function(x, z) {
    ((z-0.5)**2 - 1/12) * ((x-0.5)**2 - 1/12) / 4 - 
    ((abs(x-z) - 0.5)**4 - (abs(x-z)-0.5)**2 / 2 + 7/240) / 24
}


spl.X <- function(x, xk) {
    q = length(xk) + 2
    n = length(x)
    X = matrix(1, n, q)
    X[,2] = x
    X[,3:q] = outer(x, xk, FUN=rk)

    return(X)
}


spl.S <- function(xk) {
    q = length(xk) + 2
    S = matrix(0, q, q)
    S[3:q, 3:q] = outer(xk, xk, FUN=rk)

    return(S)
}


mat.sqrt <- function(S) {
    d = eigen(S, symmetric=TRUE)
    ev = ifelse(d$values>0, d$values, 0)
    rS = d$vectors %*% diag(ev**0.5) %*% t(d$vectors)
}

prs.fit <- function(y, x, xk, lambda) {
    q = length(xk) + 2
    n = length(x)
    
    Xa = rbind(spl.x(x, xk), mat.sqrt(spl.S(xk))*sqrt(lambda))
    y[(n+1):(n+q)] = 0
    lm(y~Xa-1)
}


am.setup <- function(x.list) {
    S = list()
    q.tot = sum(sapply(x.list, function(x) {length(unique(x))}))
    q.tot = q.tot + length(x.list)

    n = length(x.list[[1]])
    X = matrix(1, n, q.tot+1)

    indx = 2
    for (i in 1:length(x.list)) {
        x = x.list[[i]]
        xk = unique(x)
        q.x = length(xk)
        
        #Penalty matrix:
        S[[i]] = matrix(0, q.tot+1, q.tot+1)
        S[[i]][indx:(indx+q.x),indx:(indx+q.x)] = spl.S(xk)[-1,-1]

        #Data matrix:
        X[,indx:(indx+q.x)] = spl.X(x, xk)[,-1]

        indx = indx + q.x + 1
    }

    return(list(X=X, S=S))
}


fit.am <- function(y, X, S, sp) {
    dS = dim(S[[1]])
    S.tot = matrix(0, dS[1], dS[2])
    for (i in 1:length(S)) {
        S.tot = S.tot + sp[i]*S[[i]]
    }

    rS = mat.sqrt(S.tot)
    q.tot = ncol(X)
    n = nrow(X)
    X1 = rbind(X, rS)
    y1 = c(y, rep(0, q.tot))
    b = lm(y1 ~ X1-1)

    trA = sum(influence(b)$hat[1:n])
    norm = sum((y-fitted(b)[1:n])**2)
    
    return(list(model=b, gcv=norm*n / (n-trA)**2, sp=sp))
}
