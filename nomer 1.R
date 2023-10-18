#Latihan modul
#Nomer 1 Tomat
# dathe variance not equal.
fieldA=c(1.3)
fieldB=c(1.6)
nA=22
nB=24
t.test(fieldA, fieldB,n=c(nA,nB)
       alternative = "two.sided",
       mu = 0, paired = FALSE, var.equal = TRUE,
       conf.level = 0.95)

