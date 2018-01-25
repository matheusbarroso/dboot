ts.return <-
function (t0, t, R, tseries, seed, stat, sim, endcorr, n.sim, l, ran.gen, ran.args, call, norm,blocks) {
    out <- list(t0 = t0, t = t, R = R, data = tseries, seed = seed, 
        statistic = stat, sim = sim, n.sim = n.sim, call = call,blocks=blocks)
    if (sim == "scramble") 
        out <- c(out, list(norm = norm))
    else if (sim == "model") 
        out <- c(out, list(ran.gen = ran.gen, ran.args = ran.args))
    else {
        out <- c(out, list(l = l, endcorr = endcorr))
        if (!is.null(call$ran.gen)) 
            out <- c(out, list(ran.gen = ran.gen, ran.args = ran.args))
    }
    class(out) <- "boot"
    out
}
