x = [1,2,3]
y = [1,3,2]


X = [1 1
     2 3
     3 2
     4 4
     5 4
     4 5]


Xc = X .- mean(X, dims=1)

correlation(X[:,1],X[:,2],:spearman)
correlation(X[:,1],X[:,2],:pearson)
correlation(X[:,1],X[:,2],:cosine_similarity)

correlation(Xc[:,1],Xc[:,2],:spearman)
correlation(Xc[:,1],Xc[:,2],:pearson)
correlation(Xc[:,1],Xc[:,2],:cosine_similarity)
