earth.dist = function (lats, dist = TRUE) 
{
    if (class(lats) == "SpatialPoints") 
        lats <- coordinates(lats)
    name <- list(rownames(lats), rownames(lats))
    n <- nrow(lats)
    z <- matrix(0, n, n, dimnames = name)
    z <- outer(1:n, 1:n, function(i, j) deg.dist(long1 = lats[i, 
        1], lat1 = lats[i, 2], long2 = lats[j, 1], lat2 = lats[j, 
        2]))
    if (dist == TRUE) 
        z <- as.dist(z)
    return(z)
}