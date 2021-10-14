wtd_var = function(x, w, na.rm = TRUE) {

    if (na.rm) {
        s <- !is.na(x + w)
        x <- x[s]
        w <- w[s]
    }

    Weighted.Desc.Stat::w.var(x, w)
}

wtd_skewness = function(x, w, na.rm = TRUE) {

    if (na.rm) {
        s <- !is.na(x + w)
        x <- x[s]
        w <- w[s]
    }

    Weighted.Desc.Stat::w.skewness(x, w)
}

wtd_kurtosis = function(x, w, na.rm = TRUE) {

    if (na.rm) {
        s <- !is.na(x + w)
        x <- x[s]
        w <- w[s]
    }

    Weighted.Desc.Stat::w.kurtosis(x, w)
}
