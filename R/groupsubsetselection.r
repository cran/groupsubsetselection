groupsubsetselection <- function(y,x,nvarmax,nbest,nb,consind,conslb,ngv=rep(2,30))
# GSS with weights constraints R wrapper
# Working with simplified output version of the original GSS. Only 4 outputs: 
# groups, RSS, wts, and computation time (comptime) which is not necessary. 
# Input is also simplied a little, e.g. ng is deleted since ngv has already the 
# number of groups. 
#
# Input arguments:
# y: dependent variable
# x: explanatory variables, data in column in groups
# nvarmax: max size of combinations to be searched
# nbest: number of best combinations to be searhced
# nb: number of fixed variables 
# consind: the indicator vector of weight constraints
# conslb: the lower bound of weight if it is constrained
# ng: number of groups (excluding fixed variables)! this one is obsolete and removed. 
# ngv: number of variables in groups
# disp: flag of whether the result is shown after running ! Removed for simplicity.
#
# Output arguments: 
# gss outputs a list as below
# rss: the nbest rss for combinations of size up to nvarmax
# groups: the nbest combinations of size up to nvarmax found 
# wts: the weights of variables for each combination 
# wtsvars: the variables selected on variable level
# wtsvnvars: number of variables selected actually for each combination
# comptime: the total computation time for this task

# NOTE: 
# The format of explanatory variables in x:
# x = [f1, f2, ... , fn, g11, g12, ... , g1n_1, g21, g22, ... , g2n_2, .... ]
# i.e. fixed variables are in first n columns, then group1, group 2, and so on. 
# Group variables cluster together. There should be nb columns of fixed variables, 
# sum(ngv) number of grouped variables and therefore nb+sum(ngv) columns in total.
# ngv is a vector of showing the nubmer of variables for each group. 
# ngv = [n_1, n_2, ....] of length ng (ng is no more necessary). 
#
# The format of constriant vectors consind and conslb:
# They control which variable should have weight constraint and hence they have 
# the number of elements of nb+sum(ngv), i.e. the number of columns in x. 
# consind consists of only 1 or 0 meaning constrained or not constrained.
# conslb is a sequence of real numbers which are the lower bounds of the weights. 
# Example:
# consind = [0 0 0 0 0 0 1 1 0 1 0 0 1 1 ......]
# conslb = [0 1 0 0 0 0 0.1 0.2 0 -0.1 0 0 0.02 -0.05 ......]
# They show that the first 6 variables are not constrained, 7-th, 8-th , 10-th, 
# 13-th, 14-th, are constrained ..., and the 7-th variables's weight should be 
# at least 0.1 (>=), 8-th should be >= 0.2, 10-th >= -0.1, 13-th >= 0.02, 
# 14-th >= -0.05, ....
# If consind[i] is 0, then conslb[i] is irrelavant, e.g. consind[2] is 0, conslb[2]
# is not used. 
# Every variable can be weight constrained separately and the lower bounds can 
# also be tuned separately. 
#
# Output format 
# rss = [1st best mix of 1, 2nd best mix of 1, nbest-th best mix of 1, 1st best 
#        mix of 2, 2nd best mix of 2, ..., 1st best mix of nvarmax, ..., 
#        nbest-th best mix of nvarmax]
# rss is a vector of size nbest x nvarmax. If there is no combination is found, 
# then the corresponding element is a huge number. 
#
# groups = [1st best mix of 1, 2nd bset mix of 1, ..., nbest-th best mix of 1,
#         1st best mix of 2, 2nd bset mix of 2, ..., nbest-th best mix of 2,
#         ...,
#         1st best mix of nvarmax, 2nd bset mix of nvarmax, ..., nbest-th best 
#         mix of nvarmax]
# groups is a vector of size sum(1:nvarmax) x nbest. If there is no combination is
# found, then corresponding elements are 0. 
#
# wts = [weights of mix of 1
#        weights of mix of 2
#              ... ...
#        weights of mix of nvarmax]
# NB. each line is a matrix of nbest columns. If there is no variable at the 
# position, the weight is 0. The columns is the rank accordingly, i.e. the first 
# column is the 1st best, and the second column is the 2nd best and so on. There
# are various number of variables (a group is a collection of varaibles). Because
# of this, m (maximum number of variables across all groups) places is preserved 
# for weights. If the a group has less than m variables, the rest places will be
# filled with 0. The format of the weights for each combination is as follows:
# c(wt_fix1, ... ,wt_fixn, wt_var1, wt_var2,...)
# So basically, the weights of all fixed variables, then followed by weights for 
# variables in first group and so on which may be followed by zeros. The weights 
# for variables are in such an order that their corresponding variables are in 
# increasing order regarding their places in x (explanatory variable matrix, ref. 
# above) and of course the groups presented in output groups are also sorted in 
# increasing order. For example, if the second best combination of size 2 is 
# group 5 and 12 with 2 and 1 variables respectively, and the index of the variables 
# in x is 12, 13 for group 5 and 34 for group 12, the weights would be 
# c((wt_fix1, ... ,wt_fixn, wt_var12(var1 of group5), wt_var12(var2 of group5)
# , wt_var34(var1 of group12),0,0,0)
# supppose m = 3.  
# The corresponding variables are in in wtsvars in the same format as wts. 
#
# vars = [variables of mix of 1
#            variables of mix of 2
#                   ... ...
#            variables of mix of nvarmax]
# NB. each line is a matrix of nbest columns. If there is no variable at the 
# position, then it is 0. Variables include background (fixed variables, always 
# in front). The indices of the variables are the locations of the variables in 
# explanatory matrix x, in which fixed variables (background) are always in 
# front columns. Both wts and wtsvars are matrices of size 
# (nb * nvarmax + sum(1:nvarmax) * max(ngv)) * nbest
# The reason is it is not possible to determine how many variables are selected
# at runtime which is connected to memory allocation. 
#
# nvars = [ number of variables of mix of 1
#               number of variables of mix of 2
#                         ... ...
#               number of variables of mix of nvarmax]
# 
# NB. each line is a vector showing the number of variables for each combination
# selected. 
# nvars is size nvarmax x nbest.
{
  if (length(consind)!=length(conslb)) stop('Check constraint indicator vector and lower bound vector!')
  if (length(consind)!=dim(x)[2]) stop('Check consind or library!')
  ng <- length(ngv)
  ncase <- length(y)
  ntotalvars <- dim(x)[2] - nb
  if (sum(ngv)!=ntotalvars) stop('Arguments fault!')
  all <- (nvarmax+1)*nvarmax/2
  wtslen <- (all*max(ngv)+nb*nvarmax)*nbest

  
  out <- .Fortran("gss",as.double(y),as.double(t(x)),as.integer(consind), as.double(conslb),as.integer(ncase),as.integer(ntotalvars),as.integer(nb),as.integer(nvarmax),as.integer(nbest),as.integer(ng),as.integer(ngv),groups = integer(all*nbest),rss = double(nvarmax*nbest),wts = double(wtslen), as.integer(wtslen), comptime=double(1))
  
  vars <- matrix(0,all*max(ngv)+nb*nvarmax,nbest)
  nvars <- matrix(0,nvarmax,nbest)
  l1 <- c(0,cumsum((1:nvarmax)*max(ngv)+nb))+1
  for (i in 1:nvarmax) 
  {
    for (j in 1:nbest)
    {
      g <- out$groups[(sum(0:(i-1))*nbest+(j-1)*i+1):(sum(0:(i-1))*nbest+j*i)]
      g <- g[which(g>0)]
      if(length(g)>0) {
      g <- sort(g,decreasing=F)
      tmpvars <- sort(getvars.gss_ro(ngv,g),decreasing=F)
      vars[l1[i]:(l1[i]+nb+length(tmpvars)-1),j] <- c(1:nb,tmpvars+nb)
      nvars[i,j] <- nb+length(tmpvars)}
    }
  }
  r <- list(groups = out$groups,rss = out$rss, coef = matrix(out$wts,ncol=nbest), vars=vars, nvars=nvars, comptime = out$comptime)
  r
}

# Auxillary function used for gss_ro
getvars.gss_ro <- function(nos,idx)
{
t <- c(0,cumsum(nos))+1
s <- NULL
for(i in idx) s<-c(s,t[i]:(t[i+1]-1))
s
}
