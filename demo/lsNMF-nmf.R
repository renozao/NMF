# random data
y <- rmatrix(100, 20)
# add missing values
NA_values <- sample(length(y), 5)
y[ NA_values ] <- NA

x <- y
# Now a trick: as fixed dummy value (because NA values break other stuffs)
x[ NA_values ] <- 123456789

# run ls-nmf using weights that cancel out the missing values
w <- matrix(1, nrow(x), ncol(x))
w[ NA_values ] <- 0

res <- nmf(x, 3, 'ls-nmf', weight = w)

# The result can be used to input missing values
x[ NA_values ] <- fitted(res)[ NA_values ]

# NOTE (to convince yourself that the missing/dummy values are not used)
# use a common seed (only fixing RNG does not work here because the range of values in the target matrix affects the initial seed)
s <- rnmf(3, y)
res <- nmf(x, 3, 'ls-nmf', weight = w, seed = s)

# use another dummy value
x2 <- y
x2[ NA_values ] <- 987654321
res2 <- nmf(x2, 3, 'ls-nmf', weight = w, seed = s)

# results are identical
nmf.equal(res, res2)
