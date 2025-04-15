X = as.matrix(iris[,-5])
muhat = colMeans(X)
X = sweep(X, 2, muhat)
res_svd = svd(X)
col3 =  hcl.colors(3)
png("scatter_u.png", width = 900, height = 900)
plot(res_svd$u, col =col3[iris$Species], pch=16)
legend("topright", legend=levels(iris$Species), col=col3, pch=16)
dev.off()

res_eigen = eigen(t(X)%*%X)

round(t(res_eigen$vectors)%*%res_eigen$vectors, 5)


P1 = sweep(X%*%res_eigen$vectors, 2, sqrt(res_eigen$values), FUN = "/")

print(res_svd$d)
print(sqrt(res_eigen$values))
# > print(res_svd$d)
# [1] 25.099960  6.013147  3.413681  1.884524
# > print(sqrt(res_eigen$values))
# [1] 25.099960  6.013147  3.413681  1.884524

print(res_svd$v)
print(res_eigen$vectors)
# > print(res_svd$v)
#             [,1]        [,2]        [,3]       [,4]
# [1,]  0.36138659 -0.65658877  0.58202985  0.3154872
# [2,] -0.08452251 -0.73016143 -0.59791083 -0.3197231
# [3,]  0.85667061  0.17337266 -0.07623608 -0.4798390
# [4,]  0.35828920  0.07548102 -0.54583143  0.7536574
# > print(res_eigen$vectors)
#             [,1]        [,2]        [,3]       [,4]
# [1,]  0.36138659 -0.65658877  0.58202985  0.3154872
# [2,] -0.08452251 -0.73016143 -0.59791083 -0.3197231
# [3,]  0.85667061  0.17337266 -0.07623608 -0.4798390
# [4,]  0.35828920  0.07548102 -0.54583143  0.7536574

png("comparison_u.png", width = 900, height = 900)
plot(P1, res_svd$u)
abline(0, 1, col="royalblue")
dev.off()
###

X = as.matrix(iris[,-5])
res_pca = prcomp(X)

muhat = colMeans(X)
X = sweep(X, 2, muhat)

res_svd = svd(X)
PC = sweep(res_svd$u, 2, res_svd$d, FUN="*")

png("comparison_pc.png", width = 900, height = 900)
plot(PC, res_pca$x)
abline(0, 1, col="royalblue")
dev.off()

