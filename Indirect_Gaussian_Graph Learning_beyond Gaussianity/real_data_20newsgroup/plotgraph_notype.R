rm(list=ls())



library(GGally)
library(network)
library(sna)
library(ggplot2)
library(RColorBrewer)
#Need two files: WEst.txt, and varnames.txt in the right subfolder
mypath = paste(readLines("mypath.txt"), collapse=" ")
setwd(mypath)
varnames = read.table(paste(mypath, "varnames.txt", sep=""), sep="\n")
varnames = unlist(varnames)
names(varnames) = NULL
varnames = as.character(varnames)

files=readLines("plotfile.txt")
W_est = read.table(paste(mypath, files, sep=""), sep=",")
#W_est = (W_est+t(W_est))/2

m = nrow(W_est)
diag(W_est)<-rep(0, m)
# W_est[which(W_est != 0, arr.ind=TRUE)] = 1



myseed = 44
#myseed=myseed+1
set.seed(myseed) 
g = network.initialize(m) # initialize the network
net = network.adjacency(W_est, g)
network.vertex.names(net) = varnames 

k = 50


# define variable types
#net %v% "type" <-c (rep("Gaussian",k),rep("Count",k),rep("Binary",k))

# #mycolors = rainbow(21)
# labelcolor = "firebrick1" #"mediumvioletred" # "darkorchid" #"pink3" # mediumorchid2" #"black" #  "blue" # red" # mycolors[1] # #royalblue4"
# nodecolor = "phono" # mycolors[2] # 
# edgecolor = "lightsteelblue1" # white"#"bisque3"#gray90" #  "pink" # mycolors[3] #
# ggnet2(net,label=TRUE, node.size=10, 
#        label.size=1.5, label.color=labelcolor, 
#        color=nodecolor, alpha="phono", alpha.palette=0.15,  palette = "Set1",
#        edge.color=edgecolor, edge.size=.1) + guides(color=FALSE)

# ggnet2(net,label=TRUE, node.size=10, label.size=1.5, label.color="gray30", 
#        color="phono", alpha="phono", alpha.palette=0.25,  palette = "Set2",
#        edge.color="deeppink", edge.size=.15) + guides(color=FALSE)

#ggsave(paste(mypath, gsub('.txt', '_graph.pdf', files[i]), sep=""), width = 10, height = 10)


### plot a subgraph

if (unlist(strsplit(files, "_"))[2]=="gaussian.txt" ||  unlist(strsplit(files, "_"))[2]=="tukey.txt"){
  varnames = varnames[1:k]
}

if (unlist(strsplit(files, "_"))[2]=="poisson.txt"){
  varnames = varnames[(k+1):(2*k)]
}

if (unlist(strsplit(files, "_"))[2]=="bernoulli.txt" ||  unlist(strsplit(files, "_"))[2]=="lorenz.txt"){
  varnames = varnames[(2*k+1):(3*k)]
}

inds = which(colSums(W_est)==0)
sub_W_est = W_est[-inds, -inds]
sub_varnames = varnames[-inds]
sub_m = length(sub_varnames)
myseed = 20
#myseed = myseed+1
set.seed(myseed) 
sub_g = network.initialize(sub_m) # initialize the network
sub_net = network.adjacency(sub_W_est, sub_g)
network.vertex.names(sub_net) = sub_varnames 

#mycolors = rainbow(21)
labelcolor =  "firebrick2" #"mediumvioletred" "darkorchid" #"pink3" # mediumorchid2" #"black" #  "blue" # red" # mycolors[1] # #royalblue4"
nodecolor = "phono" # mycolors[2] # 
edgecolor = "gold"# white"#"bisque3"#gray90" #  "pink" # mycolors[3] #


if (dim(W_est)[1]==3*k){
sub_net %v% "Type" <-c (rep("Gaussian",k),rep("Count",k),rep("Binary",k))[ -inds]
sub_net %v% "Color" <-c (rep("phono",k),rep("steelblue",k),rep("grey",k))[ -inds]

ggnet2(sub_net,label=TRUE, node.size=degree(sub_net), 
       label.size=1.5, label.color=labelcolor, max_size = max(degree(sub_net))/5,
       color=nodecolor, alpha="phono", alpha.palette=0.7,  palette = "Set2",
       edge.color=edgecolor, edge.size=0.1,  legend.size = 6)+ guides(color=FALSE,size=FALSE)
}   

if (dim(W_est)[1]==k){
  ggnet2(sub_net,label=TRUE, node.size=degree(sub_net), 
         label.size=1.5, label.color=labelcolor, max_size = max(degree(sub_net))/5,
         color=nodecolor, alpha="phono", alpha.palette=0.7,  palette = "Set2",
         edge.color=edgecolor, edge.size=0.4 ,  legend.size = 8) + guides(shape=F,size=FALSE)+theme(panel.margin = unit(0, 'mm'))
} 
    


# ggnet2(sub_net,label=TRUE, node.size=10, label.size=1.5, label.color="gray10",
#        color="phono", alpha="phono", alpha.palette=0.25,  palette = "Set2",
#        edge.color="deeppink", edge.size=.15) + guides(color=FALSE)

ggsave(paste(mypath, gsub('.txt', '_subgraph.pdf', files), sep=""), width = 5, height = 5)

