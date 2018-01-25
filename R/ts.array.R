ts.array <-
function (n, n.sim, R, l) {
    endpt <- n - l + 1 ##total number of observer blocks of length l
    cont <- TRUE
    {
        nn <- ceiling(n.sim/l)
        lens <- c(rep(l, nn - 1), 1 + (n.sim - 1)%%l)
        st <- matrix(sample.int(endpt, nn * R, replace = TRUE), 
            R)
    }
    list(starts = st, lengths = lens)
}
