rm(list = ls())


#' Fills NA values with the mean.
#'
#' `replace_with_mean()` look for variables "x_varname"and "y_varname" in DF and
#' estimates the partial correlation adjusting by Z.
#'
#' @param x a numeric vector with some NA.
#' @returns A numeric vector without NA.
#' @examples
#' replace_with_mean(x)
#' @export
replace_with_mean <- function(x){

  mean_x <- mean(x,na.rm=T)
  x[is.na(x)] <- mean_x
  return(x)
}


#' Fills NA values with the mean by groups.
#'
#' `replace_with_mean_by_group()` Fills NA values with the mean by groups (specified by z). If
#' the group is missing, the NA value is replaced with the overall mean.
#'
#' @param x a numeric vector with some NA.
#' @param z a numeric vector with group information. No NA alloewd.
#' @returns A numeric vector without NA.
#' @examples
#' replace_with_tratified_mean(x,z)
#' @export
replace_with_tratified_mean <- function(x,z){

  # check if z has NA
  if(sum(is.na(z))!=0){
    stop("z cannot have NA")
  }

  # make a df with x and z
  XZ <- cbind(x,z)
  colnames(XZ) <- paste0("Var_",c(1:ncol(XZ)))
  groupped_means <- aggregate(Var_1 ~ ., XZ, mean)

  result <- merge(XZ, groupped_means, by = colnames(XZ)[-1] , all.x = TRUE)
  result$Var_1 <- ifelse(is.na(result$Var_1.x), result$Var_1.y, result$Var_1.x)
  result$Var_1 <- ifelse(is.na(result$Var_1), mean(XZ[,1],na.rm=TRUE), result$Var_1)
  x_new <- result$Var_1
  return(x_new)

}


#' Fills or drops NA values from a DF.
#'
#' `NA_replace()` Fills or drops NA values from a DF.
#'
#' @param DF a matrix or data frame with some NA.
#' @param replace_method replacement method for NA values: either drop or mean.
#' @returns A data frame without NA.
#' @examples
#' NA_replace(DF, "mean")
#' NA_replace(DF, Z, "mean")
#' @export
NA_replace <- function(DF,Z=NULL,replace_method = "mean"){

  # check if method was correctly specified
  if(!replace_method%in%c("mean","drop","stratified_mean")){
    stop("replace_method must be mean or drop")
  }

  # check if there is Z for stratified mean replacement
  if(missing(Z) & replace_method == "stratified_mean"){
    stop("No Z matrix was provided.")
  }


  if(replace_method == "mean"){
    na_ratio <- mean(is.na(DF))
    warning(paste0(round(na_ratio,4)*100,"% of the data was replaced."))
    new_DF <- apply(DF,2,replace_with_mean)
  }

  if(replace_method == "drop"){
    na_ratio <- mean(rowSums(is.na(DF))!=0)
    warning(paste0(round(na_ratio,4)*100,"% of the rows were dropped"))
    new_DF <- DF[rowSums(is.na(DF))==0,]
  }

  if(replace_method == "stratified_mean"){
    na_ratio <- mean(rowSums(is.na(DF))!=0)
    warning(paste0(round(na_ratio,4)*100,"% of the data was replaced."))
    new_DF <- apply(DF,2,replace_with_tratified_mean,z=Z)
  }

  return(new_DF)
}


#' Fits a PCA and estimates the errors using 95% of the explained variance.
#'
#' `pca_errors()` Fits a PCA and estimates the errors using 90% of the
#' explained variance.
#'
#' @param X a matrix or data frame of numeric variables.
#' @returns A reduced dataframe.
#' @examples
#' pca_errors(X)
#' @export
pca_errors <- function(X){

  ## Check for NA
  if(sum(is.na(X))!=0){
    stop("Cannot have NA.")
  }

  ## PCA Errors
  my_pca <- prcomp(X, scale = TRUE,
                   center = TRUE, retx = T)
  my_pca_var <- my_pca$sdev ^ 2
  my_pca_var_relative <- my_pca_var/sum(my_pca_var)
  my_pca_var_relative_cum <- cumsum(my_pca_var_relative)
  last_pc <- max(which(my_pca_var_relative_cum<=0.95))
  xhat <- t(t(my_pca$x[,1:last_pc] %*% t(my_pca$rotation[,1:last_pc])) * my_pca$scale + my_pca$center)
  pca_errors <- X-xhat
  return(pca_errors)
}


#' Identifies outliers using the specified method.
#'
#' `outlier_detection_3sigma()` Identifies outliers which are 3 std. dev.
#' higher than the mean.
#'
#' @param x a matrix or data frame of numeric variables.
#' @returns A Boolean vector indicating the outliers.
#' @examples
#' outlier_detection_3sigma(X)
#' @export
outlier_detection_3sigma <- function(x){

  sigma <- sqrt(var(x))
  lower <- mean(x)-3*sigma
  upper <- mean(x)+3*sigma
  outliers <- x>upper|x<lower
}


#' Identifies outliers using the specified method.
#'
#' `outlier_detection_iqr()` Identifies outliers which are 1.5 IQR higher than
#' the mean.
#'
#' @param x a matrix or data frame of numeric variables.
#' @returns A Boolean vector indicating the outliers.
#' @examples
#' outlier_detection_iqr(X)
#' @export
outlier_detection_iqr <- function(x){

  qtiles <- quantile(x,c(0.25,0.75))
  iqr <- qtiles[2]-qtiles[1]
  lower <- qtiles[1]-1.5*iqr
  upper <- qtiles[2]+1.5*iqr
  outliers <- x>upper|x<lower
}


#' Identifies outliers using the specified method.
#'
#' `detect_outliers()` Identifies outliers using either 3 sigma or IQR. By
#' default, a PCA is estimated first and the outliers are identified in the
#' error term. This option can be ignored if PCA is set to FALSE.
#'
#' @param df a matrix or data frame of numeric variables.
#' @param outlier_method a matrix or data frame of numeric variables.
#' @param drop_method a matrix or data frame of numeric variables.
#' @param method a matrix or data frame of numeric variables.
#' @returns A reduced dataframe.
#' @examples
#' detect_outliers(df=df,outlier_method="IQR",drop_method="obs", PCA=TRUE)
#' @export
detect_outliers <- function(df,outlier_method="IQR",drop_method="obs", PCA=TRUE){

  # check if outlier and drop method are correct
  if(!outlier_method%in%c("IQR","3Sigma")){
    stop("outlier_method must be IQR or 3Sigma")
  }

  if(!drop_method%in%c("row","obs")){
    stop("drop_method must be row or obs")
  }
  # save df
  new_df <- df

  # if PCA is true, get the errors
  if(PCA==TRUE){
  df <- pca_errors(df)
  }

  # detect outliers
  if(outlier_method=="IQR"){
    Outliers <- apply(df,2,outlier_detection_iqr)
  }

  if(outlier_method=="3Sigma"){
    Outliers <- apply(df,2,outlier_detection_3sigma)
  }

  if(drop_method=="obs"){
    new_df[Outliers] <- NA

  }

  if(drop_method=="row"){
    Outliers <- rowSums(Outliers)>0
    new_df[Outliers==TRUE,] <- NA
  }

  return(list(Outliers=Outliers,new_df=new_df))

}


#' Estimates the partial correlation of X and Y given Z.
#'
#' `partial_cor_fun()` look for variables "x_varname"and "y_varname" in DF and
#' estimates the partial correlation adjusting by Z.
#'
#' @param x_varname a character string indicating the name of X variable.
#' @param y_varname a character string indicating the name of Y variable.
#' @param Z a data frame with the conditioning variables.
#' @param DF the data frame with X and Y variables.
#' @param cor_method a character string indicating which correlation coefficient (or covariance) is to be computed.
#' One of pearson (default), spearman or kendall.
#' @returns A data frame with the estimated partial correlation and the P. Value.
#' @examples
#' partial_cor_fun("x_varname", "y_varname", Z, DF, "pearson")
#' @export
partial_cor_fun <- function(x_varname,y_varname,Z,DF,cor_method){

  # Get X and Y variables with the index or variable names
  X <- as.matrix(DF[,x_varname])
  Y <- as.matrix(DF[,y_varname])

  # Drop rows with any NA because pcor.test does not allow missing values
  na_ind <- rowSums(is.na(cbind(X,Y,Z))) != 0
  X_notNA <- X[!na_ind]
  Y_notNA <- Y[!na_ind]
  Z_notNA <- Z[!na_ind,]


  # Test for equal variable
  if (identical(X_notNA, Y_notNA)) {
    rho <- 1
    pval <- 0
  } else {

    # Estimate the partial correlation if it's not the same variable
    PcorTestOutput <- ppcor::pcor.test(X_notNA, Y_notNA, Z_notNA, method = cor_method)
    rho <- PcorTestOutput$estimate
    pval <- PcorTestOutput$p.value
  }

  # Build output data frame
  OutPut <- data.frame(rho=rho,pval=pval)
  return(OutPut)
}


#' Estimates the partial correlation matrix of a data frame DF.
#'
#' `mat_partial_cor_fun()` estimates the partial correlation matrix of DF
#' adjusting by Z variables.
#'
#' @param DF the data frame with X and Y variables.
#' @param Z a data frame with the conditioning variables.
#' @param cor_method a character string indicating which correlation coefficient (or covariance) is to be computed.
#' One of pearson (default), spearman or kendall.
#' @returns A data frame with the estimated partial correlation and the P. Value.
#' @seealso [partial_cor_fun].
#' @examples
#' mat_partial_cor_fun(Z, DF, "pearson")
#' @export
mat_partial_cor_fun <- function(DF,Z,cor_method){

  # Create data frame with variable names
  PCor_mat <- expand.grid(x_varname=colnames(DF),y_varname=colnames(DF))

  # Map PCor
  PCor_mat$PCor <- purrr::map2(.x=PCor_mat$x_varname,.y=PCor_mat$y_varname,
                        .f=partial_cor_fun,Z=Z,DF=DF,
                        cor_method=cor_method)

  # Get rho and p. values
  PCor_mat$rho <- purrr::map_dbl(.x=PCor_mat$PCor,.f="rho")
  PCor_mat$pval <- purrr::map_dbl(.x=PCor_mat$PCor,.f="pval")

  # Get the Rho Matrix
  Rho <- tidyr::pivot_wider(data=PCor_mat[,c("x_varname","y_varname","rho")],
                     names_from=y_varname,values_from=rho)
  Rho <- as.matrix(Rho[, -1])  # Remove the Xindx column
  rownames(Rho) <- colnames(Rho)

  # Get the P. Values Matrix
  PVal <- tidyr::pivot_wider(data=PCor_mat[,c("x_varname","y_varname","pval")],
                      names_from=y_varname,values_from=pval)
  PVal <- as.matrix(PVal[, -1])  # Remove the Xindx column
  rownames(PVal) <- colnames(PVal)

  # Return output
  OutList <- list(Rho=Rho,PVal=PVal)
  return(OutList)
}


#' Hierarchical cluster analysis using TOM Distance.
#'
#' `hc_fun()` does a hierarchical cluster analysis on the TOM distance matrix
#' D and chooses the best number of clusters K by maximizing the silhouette
#' values for each possible number of clusters.
#'
#' @param D the distance matrix of the Lower Level Network calculated with TOM.
#' @param linkage type of linkage for the cluster algorithm.
#' @returns A list with the new group index, the silhouette plot, and the best K.
#' @examples
#' hc_fun(D)
#' @export
hc_fun <- function(D,linkage="ward.D"){

  # Max possible number of clusters p-1
  k_max <- ncol(D)-1

  # Hierarchical Cluster
  clusters <- stats::hclust(as.dist(D), method = linkage)

  # Calculate Silhouette values for all k values
  mean_sil <- rep(0,k_max) # To store all mean(si)

  # Calculate Silhouette value for each k above 1
  for(i in 2:k_max){
    sil_values <- cluster::silhouette(stats::cutree(clusters,i),as.dist(D))
    mean_sil[i] <- as.data.frame(sil_values)%>%dplyr::group_by(cluster)%>%
      dplyr::summarise(si=mean(sil_width))%>%ungroup()%>%
      dplyr::summarise(si=mean(si))%>%dplyr::pull(si)
  }

  # Find K that maximizes sk
  best_K <- which.max(mean_sil)

  # Drop first value (when k=1) because we don't calculate it
  mean_sil[1] <- NA

  # Make the Silhouette plot
  SilhouettePlot <-
    ggplot2::ggplot()+
    ggplot2::geom_line(aes(x=c(2:length(mean_sil)),y=mean_sil[-1]),color = "red")+
    ggplot2::geom_vline(xintercept = best_K, linetype = "dotted", color="black")+
    ggplot2::labs(x="# of clusters", y="Silhouette Value")+
    ggplot2::theme_minimal()+
    ggplot2::scale_x_continuous(breaks = c(seq(from=0,to=k_max,by=floor(k_max/4)),3),
                       limits= c(0,k_max))

  # Get the Best Clusters
  ClGroupIndex <- stats::cutree(clusters,best_K)

  # Return
  OutList <- list(ClGroupIndex = ClGroupIndex,
                  SilhouettePlot = SilhouettePlot,
                  K = best_K)
  return(OutList)
}


#' K medoids cluster analysis using TOM Distance.
#'
#' `pam_fun()` does a K medoids cluster analysis on the TOM distance matrix
#' D and chooses the best number of clusters K by maximizing the silhouette
#' values for each possible number of clusters.
#'
#' @param D the distance matrix of the Lower Level Network calculated with TOM.
#' @param linkage type of linkage for the cluster algorithm.
#' @returns A list with the new group index, the silhouette plot, and the best K.
#' @examples
#' pam_fun(D)
#' @export
pam_fun <- function(D,linkage="euclidean"){

  # Set seed for replicability
  set.seed(42)

  # Max possible number of clusters p-1
  k_max <- ncol(D)-1
  mean_sil <- rep(0,k_max) # To store all mean(si)

  for(i in 2:k_max){
    # Runs PAM for each k above 1
    Clusters <-
      cluster::pam(x=as.dist(D),metric=linkage,k=i,diss=TRUE)
    summaryPAM <- summary(Clusters)

    # Keeps Silhouette value
    mean_sil[i]  <- mean(summaryPAM$silinfo$clus.avg.widths)
  }

  # Find K that maximizes sk
  best_K <- which.max(mean_sil)

  # Drop first value (when k=1) because we don't calculate it
  mean_sil[1] <- NA

  # Make the Silhouette plot
  SilhouettePlot <-
    ggplot2::ggplot()+
    ggplot2::geom_line(aes(x=c(2:length(mean_sil)),y=mean_sil[-1]),color = "red")+
    ggplot2::geom_vline(xintercept = best_K, linetype = "dotted", color="black")+
    ggplot2::labs(x="# of clusters", y="Silhouette Value")+
    ggplot2::theme_minimal()+
    ggplot2::scale_x_continuous(breaks = c(seq(from=0,to=k_max,by=floor(k_max/4)),3),
                       limits= c(0,k_max))

  # Get the Best Clusters
  ClGroupIndex <- cluster::pam(x=as.dist(D),metric=linkage,k=best_K)$clustering

  # Return
  OutList <- list(ClGroupIndex = ClGroupIndex,
                  SilhouettePlot = SilhouettePlot,
                  K = best_K)
  return(OutList)
}


#' Cluster analysis on the distance matrix of TOM.
#'
#' `clusters_fun()` does a cluster analysis on the TOM distance matrix
#' D and chooses the best number of clusters K by maximizing the silhouette
#' values for each possible number of clusters.
#'
#' @param W the Lower Level Network calculated with TOM.
#' @param ClGroup group of variables that will be clustered.
#' If not specified, it will use all variables.
#' @param ClType one of "HC" (Hierarchical Clustering), "PAM"
#' (K Medoids clustering) or "NOCL" (No Cluster).
#' @param GroupIndx_df a data frame with the groups information.
#' @returns A list with the new group index, the silhouette plot, and the best K.
#' @examples
#' clusters_fun(W,ClGroup,ClType,GroupIndx_df)
#' @export
clusters_fun <- function(W,ClGroup,ClType,GroupIndx_df){

  # Distance Matrix
  D <- 1-W

  # Names of Variables inside the Cluster Group
  ClusteredVarNames <- GroupIndx_df$VarName[GroupIndx_df$GroupName==ClGroup]

  # Names of Variables outside the Cluster Group
  NotClusteredVarNames <- GroupIndx_df$VarName[GroupIndx_df$GroupName!=ClGroup]

  # Subset of W for all variables inside the Cluster Group
  ClDistanceMatrix <- D[ClusteredVarNames,ClusteredVarNames]

  # Get the clusters by Hierarchical Clustering (HC)
  if(ClType == "HC"){Clusters <- hc_fun(D=ClDistanceMatrix)}

  # Get the clusters by Partitioning Around Medoids (PAM)
  if(ClType == "PAM"){Clusters <- pam_fun(D=ClDistanceMatrix)}

  # Cluster Index
  ClGroupIndex <- paste0(ClGroup,Clusters$ClGroupIndex)

  # Get the cluster variable names and its corresponding cluster group
  ClustersIndex <- data.frame(VarName = ClusteredVarNames,
                              cluster = ClGroupIndex)

  # Replace old groups
  GroupIndx_New <- merge(GroupIndx_df,ClustersIndex,by = "VarName",all.x = T)
  GroupIndx_New$GroupName <- ifelse(GroupIndx_New$GroupName == ClGroup,
                                    GroupIndx_New$cluster,
                                    GroupIndx_New$GroupName)
  GroupIndx_New$cluster <- NULL
  GroupIndx_New <- GroupIndx_New[order(GroupIndx_New$GroupName),]

  # Return
  Out <- list(GroupIndex_New = GroupIndx_New,
              SilhouettePlot = Clusters$SilhouettePlot,
              K = Clusters$K)
  return(Out)
}


#' Upper Level Network.
#'
#' `uln_fun()` calculates the Upper Level Network using the mean w_{ij} values
#' of all the units in the groups.
#'
#' @param W the Lower Level Network calculated with TOM.
#' @param groupsID data frame with the groups information.
#' @returns A matrix with the weights W of the Upper Level NetworK.
#' @examples
#' uln_fun(groupsID,W)
#' @export
uln_fun <- function(W,groupsID){

  # Create the matrix for the upper level network
  nG <- length(unique(groupsID$GroupName))
  W_ULN <- matrix(NA,nG,nG)
  colnames(W_ULN) <- unique(unique(groupsID$GroupName))
  rownames(W_ULN) <- unique(unique(groupsID$GroupName))

  # Calculate the upper level w_ij as the average W_ij of the groups
  for(col in colnames(W_ULN)){
    for (row in rownames(W_ULN)){
      colName <- col
      rowName <- row
      colVarNames <- groupsID$VarName[groupsID$GroupName==colName]
      rowVarNames <- groupsID$VarName[groupsID$GroupName==rowName]
      W_ULN[rowName,colName] <- mean(W[rowVarNames,colVarNames])
    }
  }

  # Replace diagonal with 1's
  W_ULN <- W_ULN - diag(diag(W_ULN)) + diag(ncol(W_ULN))
  return(W_ULN)
}


#' Multi Layer Network Analysis
#'
#' `multi_network()` This code applies the Multivariate Network Analysis
#' methodology used in Guan, Cheng, Koo.
#'
#'First, a lower-level network based on correlation matrix is constructed to
#'represent interaction patterns between each pair of variables. Soft
#'thresholding (power function and Topological Overlap Measure) was applied to
#'the matrix to emphasize strong connections and remove spurious connections.
#'
#'Then, to improve visualization of extensive network structure from the high
#'dimensional data, biologically predefined domains (i.e., blood lipid,
#'cytokines, etc.) and hierarchical cluster analysis is used to define modules
#'of potentially highly interconnected variables. The resulted upper-level
#'network represents average connectivity between modules.
#'
#' @param DF a data frame used to estimate the networks.
#' @param Groups a vector with the corresponding group of each variable. If
#' not specified, the groups will be optimized with a clustering algorithm.
#' @param Z a data frame for adjusting the partial correlations. If not
#' specified, we will use correlation instead of partial correlation.
#' @param cor_method a character variable that specifies the correlation
#' method. Default is "pearson".
#' @param power the power function for the soft threshold. Default is 6.
#' @param TOMType type of TOM function. Default is unsigned.
#' @param TOMDenom type of function in the denominator. Default is min.
#' @param ClGroup a specific group to do the clustering analysis.
#' @param ClType one of "HC" (Hierarchical Clustering), "PAM"
#' (K Medoids clustering) or "NOCL" (No Cluster). Default if "HC".
#' @param writeCSV set TRUE if you want to save your output into CSV files. The
#' default is FALSE.
#' @returns A list with the Upper and Lower Level Networks, the correlation
#' matrix, a data frame with the groups index, and the best K and the
#' silhouette plot if any clustering algorithm was used.
#' @examples
#' multi_network(DF,Groups,Z,ClGroup="Group 1",ClType="HC")
#' multi_network(DF,ClType="HC")
#' multi_network(DF,Groups,ClType="NOCL")
#' @export
multi_network <- function(DF,Groups=NULL,Z=NULL,
                              cor_method="pearson",
                              power=6,
                              TOMType="unsigned",
                              TOMDenom="min",
                              ClGroup="G1",ClType="HC",
                              writeCSV=FALSE){

  # Check if all variables in DF are numeric
  all_numeric <- sapply(DF, is.numeric)
  if (!all(all_numeric)) {
    stop("All columns must be numeric.")
  }

  # Check Groups info was not provided
  if(is.null(Groups)){
    GroupIndx = data.frame(GroupName="G1",
                           VarName=colnames(DF))
  } else {

    # Check if dimensions of DF and Groups are the same
    if(length(Groups)!=ncol(DF)){
      stop("The group information does not match the dimensions of the data frame. Each variable must be uniquely placed within a group.")
    }

    # Create groups data frame using group info
    GroupIndx = data.frame(GroupName=Groups,
                           VarName=colnames(DF))
  }

  # Check if Cor_method was correctly specified
  if(!cor_method%in%c("pearson","kendall","spearman")){
    stop("Correlation method must be pearson, kendall or spearman.")
  }

  # Check if Z was provided
  if(missing(Z)){
    message("No Z matrix was provided. Correlation will be used instead of Partial Correlation.")
    Rho = cor(DF,method=cor_method,use="pairwise.complete.obs")
  } else {

    # If Z was provided use the partial correlation
    Rho <- mat_partial_cor_fun(DF=DF,Z=Z,cor_method=cor_method)$Rho
  }

  # Power Function
  A <- (abs(Rho))^power

  # TOM
  W <-
    WGCNA::TOMsimilarity(adjMat=A,
                  TOMType = TOMType,
                  TOMDenom = TOMDenom)
  colnames(W) <- colnames(A)
  rownames(W) <- rownames(A)

  # Check for cluster type

  if(!ClType%in%c("HC","PAM","NOCL")){
    stop("Cluster type must be HC, PAM or NOCL. NOCL is for No Cluster.")
  }

  if(ClType%in%c("HC","PAM")){

    # Check if ClGroup was correctly specified
    if(!ClGroup%in%GroupIndx$GroupName){
      stop("Cluster group is not in groups.")
    }

    # Get new W Matrix with clusters
    Out <-
      clusters_fun(W=W,
                   GroupIndx_df=GroupIndx,
                   ClGroup=ClGroup,ClType=ClType)
    GroupIndx <- Out$GroupIndex_New

  }

  # Check if there are at least 2 groups
  ng <- length(unique(GroupIndx$GroupName))
  if(ng<=1){stop("Must be at least 2 groups")}

  # Upper Level Network
  ULN <- uln_fun(W=W,groupsID=GroupIndx)

  # Return objects if no CL
  if(ClType=="NOCL"){
    ReturnList <- list(Rho=Rho,LLN=W,ULN=ULN,GroupIndx=GroupIndx)
  }

  # Return objects if CL
  if(ClType%in%c("HC","PAM")){
    ReturnList <- list(Rho=Rho,LLN=W,ULN=ULN,GroupIndx=GroupIndx,
                       Kbest=Out$K,SilhouettePlot=Out$SilhouettePlot)
  }

  # write CSV of outputs
  if(writeCSV){write.csv(W,"LLN.csv"); write.csv(ULN,"ULN.csv")}

  # Return objects
  return(ReturnList)
}


#' Distance Transformation
#'
#' `distance_transformation()` calculates the distance using Inverse or
#' Logit transformations.
#'
#' @param W a matrix with the weights of the Network calculated with TOM.
#' @param transform Either Inv for inverse, or Log for -log10 transformation.
#' @returns A new Distance Matrix.
#' @examples
#' distance_transformation(W,"Inv")
#' distance_transformation(W,"Log")
#' @export
distance_transformation <- function(W,transform="Inv"){
  # check if transform is correctly specified
  if(!transform%in%c("Inv","Log")){stop("Transform must be Inv of Log.")}

  # Inverse transformation
  if(transform=="Inv"){W_transformed = 1/W}

  # Log transformation
  if(transform=="Log"){W_transformed = -log10(W)}

  # Diagonal to 0
  W_transformed <- W_transformed-diag(diag(as.matrix(W_transformed)))
  return(W_transformed)
}


#' Floyd Warshall
#'
#' `floyd_warshall()` calculates the minimum distance path between the nodes.
#' Outputs two matrices: min_Distance is the new Matrix containing the minimum
#' distance, and Next contains a matrix with the minimum distance neighbor
#'
#' @param Distance a V*V Distance Matrix.
#' @returns A list with the minimum distance matrix, and the Next matrix with
#' the closest path.
#' @examples
#' floyd_warshall(Distance)
#' @export
floyd_warshall <- function(Distance) {
  # Number of vertices in the graph
  V <- nrow(Distance)

  # Initialize the distance and next matrices
  min_Distance <- Distance
  Next <- matrix(0, nrow = V, ncol = V)

  # Find the shortest path for all pairs of vertices
  for (k in 1:V) {
    for (i in 1:V) {
      for (j in 1:V) {
        if (min_Distance[i, k] + min_Distance[k, j] < min_Distance[i, j]) {
          min_Distance[i, j] <- min_Distance[i, k] + min_Distance[k, j]
          Next[i, j] <- k
        }
      }
    }
  }

  # Return the distance and next matrices
  return(list(min_Distance = min_Distance, Next = Next))
}


#' Get Path
#'
#' `get_path()` gets the best path between two nodes.
#'
#' @param Next a Matrix with the next closest neighbor between two nodes.
#' @param start the starting node.
#' @param end the ending node.
#' @returns a string with the best path from the starting node to the ending node.
#' @examples
#' floyd_warshall(Distance)
#' @export
get_path <- function(Next, start, end) {
  initial_hop <- Next[start, end]
  if(initial_hop!=0){
    hop = initial_hop
    hops <- c(hop)
    # Hops To the end
    while (hop != 0) {
      new_hop <- Next[hop, end]
      if(new_hop!=0){
        hops <- c(hops,new_hop)
      }
      hop <- new_hop
    }
    # Hops from the start
    hop <- hops[1]
    while (hop != 0) {
      new_hop <- Next[start, hop]
      if(new_hop!=0){
        hops <- c(new_hop, hops)
      }
      hop <- new_hop
    }
    path <- c(start,hops,end)
    path <- paste(path, collapse = ", ")
  }
  if(initial_hop==0){
    path <- c(start,end)
    path <- paste(path, collapse = ", ")
  }
  return(paste0("[",path,"]"))
}


#' Path Finder
#'
#' `path_finder()` gets the best path between each pair of nodes.
#'
#' @param Next a Matrix with the next closest neighbor between two nodes.
#' @returns a Matrix with the best paths for each pair of nodes.
#' @examples
#' path_finder(Next)
#' @export
path_finder <- function(Next){
  V = nrow(Next)
  Path_Mat = matrix(NA,V,V)
  for(i in 1:V){
    for(j in 1:V){
      if(j!=i){
        Path_Mat[i,j] = get_path(Next,start=i,end=j)
      }
      if(j==i){
        Path_Mat[i,j] = 0
      }
    }
  }
  return(Path_Mat)
}


#' Split Path
#'
#' `split_path()` split path into pairs of nodes.
#'
#' @param path a path resulting from get_path().
#' @returns a string with the hops, separated by a comma.
#' @examples
#' split_path(path)
#' @export
split_path <- function(path){
  numbers <- str_extract_all(path, "\\d+")[[1]]
  n <- length(numbers)
  pairs <- paste(numbers[c(1:(n-1))],numbers[2:n], sep=",")
  pairs <- paste0("[",pairs,"]")

  return(pairs)
}


#' Network Paths
#'
#' `networks_paths()` calculates a distance matrix using the specified
#' transformation, then calculates the minimum distance matrix and the optimum
#' paths for each pair of nodes.
#'
#' @param W the Lower Level Network calculated with TOM.
#' @param transform Either Inv for inverse 1/W or Log for -log10 transformation.
#' @returns a list with distance matrix, the min distance matrix, a Next
#' matrix with the closest neighbor, and a path matrix with the best paths for e
#' ach pair of nodes.
#' @examples
#' split_path(path)
#' @export
networks_paths <- function(W,transform="Inv"){
  Distance <- distance_transformation(W,"Inv")
  FW_output <- floyd_warshall(Distance)
  Next = FW_output$Next
  min_Distance = FW_output$min_Distance
  Paths <- path_finder(Next)
  return(list(Distance=Distance,
              min_Distance=min_Distance,
              Next=Next,
              Paths=Paths))

}


#' Key Paths
#'
#' `key_paths()` finds the key paths in a vector of paths by breaking them by
#' pairs and counting the amount of times each pair appears.
#'
#' @param path a path resulting from get_path().
#' @returns a data frame with all the hops with the number of times it has been
#' used, and a plot.
#' @examples
#' split_path(path)
#' @export
key_paths <- function(paths){
  # split all paths from Paths Matrix
  path_pairs_vector <- unlist(map(paths,split_path))

  # Drops paths from diagonal
  pairs_with_NA <- grep(pattern="NA", x=path_pairs_vector)
  pairs_with_0_start <- grep(pattern="\\[0", x=path_pairs_vector)
  pairs_with_0_end <- grep(pattern=",0]", x=path_pairs_vector)
  pairs_with_0_or_NA <- unique(c(pairs_with_NA,pairs_with_0_end,pairs_with_0_start))
  if(length(pairs_with_0_or_NA)>0){
    path_pairs_vector <- path_pairs_vector[-pairs_with_0_or_NA]
  }

  # Turn into a data frame
  key_paths <- as.data.frame(table(path_pairs_vector))
  colnames(key_paths) <- c("key_path","N")

  # Make a bar plot
  key_paths_plot <-
    ggplot(key_paths)+
    geom_col(aes(y = reorder(key_path, -N),x=N,fill=N),alpha=.7)+
    labs(y="",x="",fill="")+
    theme_minimal()+
    scale_fill_gradient(low = "black", high = "darkred")

  return(list(key_paths=key_paths,key_paths_plot=key_paths_plot))

}


#' Multi Layer Network Analysis
#'
#' `multi_network_boot()` This code applies the Multivariate Network Analysis
#' methodology used in Guan, Cheng, Koo, using bootstrap to generate a
#' distribution for each weight in the Upper and Lower Level Networks.
#'
#'First, a lower-level network based on correlation matrix is constructed to
#'represent interaction patterns between each pair of variables. Soft
#'thresholding (power function and Topological Overlap Measure) was applied to
#'the matrix to emphasize strong connections and remove spurious connections.
#'
#'Then, to improve visualization of extensive network structure from the high
#'dimensional data, biologically predefined domains (i.e., blood lipid,
#'cytokines, etc.) and hierarchical cluster analysis is used to define modules
#'of potentially highly interconnected variables. The resulted upper-level
#'network represents average connectivity between modules.
#'
#' @param DF a data frame used to estimate the networks.
#' @param Groups a vector with the corresponding group of each variable. If
#' not specified, the groups will be optimized with a clustering algorithm.
#' @param Z a data frame for adjusting the partial correlations. If not
#' specified, we will use correlation instead of partial correlation.
#' @param cor_method a character variable that specifies the correlation
#' method. Default is "pearson".
#' @param B the number of bootstrap samples.
#' @param prop proportion of the bootstrap sample.
#' @param replace set TRUE if sample with replacement is preferred.
#' @param seed the number of bootstrap samples.
#' @param power the power function for the soft threshold. Default is 6.
#' @param TOMType type of TOM function. Default is unsigned.
#' @param TOMDenom type of function in the denominator. Default is min.
#' @param ClGroup a specific group to do the clustering analysis.
#' @param ClType one of "HC" (Hierarchical Clustering), "PAM"
#' (K Medoids clustering) or "NOCL" (No Cluster). Default if "HC".
#' @param writeCSV set TRUE if you want to save your output into CSV files. The
#' default is FALSE.
#' @returns A list with the Upper and Lower Level Networks (both the mean
#' matrix and bootstrap 3D matrix), a 3D bootstrap correlation matrix,
#' a data frame with the groups index, and the best K and the silhouette plot
#' if any clustering algorithm was used.
#' @examples
#' multi_network(DF,Groups,Z,ClGroup="Group 1",ClType="HC")
#' multi_network(DF,ClType="HC")
#' multi_network(DF,Groups,ClType="NOCL")
#' @export
multi_network_boot <- function(DF,Groups=NULL,Z=NULL,
                               B=1000,seed=NULL,prop=0.8,replace=F,
                               cor_method="pearson",
                               power=6,
                               TOMType="unsigned",
                               TOMDenom="min",
                               ClGroup="G1",ClType="HC",
                               writeCSV=FALSE){

  # Check if b is numeric
  if (!is.numeric(B)) {
    stop("b must be numeric.")
  }

  # Check if all variables in DF are numeric
  all_numeric <- sapply(DF, is.numeric)
  if (!all(all_numeric)) {
    stop("All columns must be numeric.")
  }

  # Check Groups info was not provided
  if(is.null(Groups)){
    GroupIndx = data.frame(GroupName="G1",
                           VarName=colnames(DF))
  } else {

    # Check if dimensions of DF and Groups are the same
    if(length(Groups)!=ncol(DF)){
      stop("The group information does not match the dimensions of the data frame. Each variable must be uniquely placed within a group.")
    }

    # Create groups data frame using group info
    GroupIndx = data.frame(GroupName=Groups,
                           VarName=colnames(DF))
  }

  # Check if Cor_method was correctly specified
  if(!cor_method%in%c("pearson","kendall","spearman")){
    stop("Correlation method must be pearson, kendall or spearman.")
  }

  # Gen seeds for bootstrap
  ## Check if seed was provided
  if(missing(seed)){
    boot_seeds = sample(c(1:92173),10*B)

  } else {
    set.seed(seed)
    boot_seeds = sample(c(1:92173),10*B)
  }

  # Gen matrices for bootstrapped samples
  n_df = nrow(DF)
  p = ncol(DF)
  Rho_boostraps = array(NA,dim = c(p,p,B))
  A_boostraps = array(NA,dim = c(p,p,B))
  TOM_boostraps = array(NA,dim = c(p,p,B))
  dimnames(TOM_boostraps)[[1]] <- colnames(DF)
  dimnames(TOM_boostraps)[[2]] <- colnames(DF)

  # Check if Z was provided
  b <- 0
  cont <- 0
  if(missing(Z)){
    message("No Z matrix was provided. Correlation will be used instead of Partial Correlation.")
    while (cont < B){
      b <- b+1
      cont <- cont+1

      set.seed(boot_seeds[b])
      indx = sample(c(1:n_df),floor(prop*n_df),replace = replace)
      # Check for warnings
      tt <- tryCatch(mat_partial_cor_fun(DF[indx,],Z[indx,],"pearson"),error=function(e) e, warning=function(w) w)$message
      if(!is.null(tt)){
        cont <- cont - 1
        next
      }
      # Correlation
      Rho_boostraps[,,cont] = cor(DF,method=cor_method,use="pairwise.complete.obs")
      # Power Function
      A_boostraps[,,cont] = abs(Rho_boostraps[,,cont])^power
      # TOM
      TOM_boostraps[,,cont] = WGCNA::TOMsimilarity(adjMat=A_boostraps[,,cont],
                                                TOMType = TOMType,
                                                TOMDenom = TOMDenom,
                                                verbose = 0)
    }

  } else {
    while (cont < B){
    b <- b+1
    cont <- cont+1

    # If Z was provided use the partial correlation
    set.seed(boot_seeds[b])
    indx = sample(c(1:n_df),floor(prop*n_df),replace = replace)
    # Check for warnings
    tt <- tryCatch(mat_partial_cor_fun(DF[indx,],Z[indx,],"pearson"),
                   error=function(e) e,
                   warning=function(w) w)$message
    if(!is.null(tt)){
      cont <- cont - 1
      next
    }
    # Correlation
    Rho_boostraps[,,cont] = mat_partial_cor_fun(DF[indx,],Z[indx,],cor_method)$Rho
    # Power Function
    A_boostraps[,,cont] = abs(Rho_boostraps[,,cont])^power
    # TOM
    TOM_boostraps[,,cont] = WGCNA::TOMsimilarity(adjMat=A_boostraps[,,cont],
                                              TOMType = TOMType,
                                              TOMDenom = TOMDenom,
                                              verbose = 0)
    }
  }


  # TOM
  W <- apply(TOM_boostraps, c(1,2), mean)

  # Check for cluster type

  if(!ClType%in%c("HC","PAM","NOCL")){
    stop("Cluster type must be HC, PAM or NOCL. NOCL is for No Cluster.")
  }

  if(ClType%in%c("HC","PAM")){

    # Check if ClGroup was correctly specified
    if(!ClGroup%in%GroupIndx$GroupName){
      stop("Cluster group is not in groups.")
    }

    # Get new W Matrix with clusters
    Out <-
      clusters_fun(W=W,
                   GroupIndx_df=GroupIndx,
                   ClGroup=ClGroup,ClType=ClType)
    GroupIndx <- Out$GroupIndex_New

  }

  # Check if there are at least 2 groups
  ng <- length(unique(GroupIndx$GroupName))
  if(ng<=1){stop("Must be at least 2 groups")}

  # Upper Level Network
  ULN_boostraps = array(NA,dim = c(ng,ng,B))
  dimnames(ULN_boostraps)[[1]] <- unique(GroupIndx$GroupName)
  dimnames(ULN_boostraps)[[2]] <- unique(GroupIndx$GroupName)
  for (b in 1:B){
    ULN_boostraps[,,b] = uln_fun(groupsID=GroupIndx,W=TOM_boostraps[,,b])
  }
  ULN <- apply(ULN_boostraps, c(1,2), mean)

  # Return objects if no CL
  if(ClType=="NOCL"){
    ReturnList <- list(Rho_b=Rho_boostraps,
                       LLN=W,ULN=ULN,GroupIndx=GroupIndx,
                       LLN_b=TOM_boostraps,ULN_b=ULN_boostraps)
  }

  # Return objects if CL
  if(ClType%in%c("HC","PAM")){
    ReturnList <- list(Rho_b=Rho_boostraps,
                       LLN=W,ULN=ULN,GroupIndx=GroupIndx,
                       LLN_b=TOM_boostraps,ULN_b=ULN_boostraps,
                       Kbest=Out$K,SilhouettePlot=Out$SilhouettePlot)
  }

  # write CSV of outputs
  if(writeCSV){write.csv(W,"LLN.csv"); write.csv(ULN,"ULN.csv")}

  # Return objects
  return(ReturnList)
}


#' Wilcoxon Paired Test
#'
#' `Wilcoxon_paired_test()` performs a Wilcoxon test on the median of the
#' bootstrapped distributions of the weights, and compares the Networks of
#' two groups.
#'
#' @param X a Bootstrapped Network Matrix.
#' @param Y a Bootstrapped Network Matrix.
#' @param paired either TRUE of FALSE for Wilcoxon Paired Test.
#' @returns A list with the output of the wilcoxon test for each weight, and a
#' plot with the confidence interval and the median.
#' @examples
#' split_path(path)
#' @export
Wilcoxon_paired_test <- function(X,Y,paired=T){

  # Check if X and Y have the same dimensions
  if(min(dim(X)==dim(Y))==0){
    stop("Both Networks must have the same dimmensions")
  }

  # Check if X and Y have 3 Dimmensions
  if(length(dim(X))!=3){
    stop("Networks must have 3 dimmensions")
  }

  # matrix parameters and variable names
  p=dim(X)[2]
  n=dim(X)[3]
  var_names = dimnames(X)[[1]]

  # Wilcoxon test results matrix
  Wilcoxon = data.frame(from = rep("FROM",p^2),
                        to = rep("TO",p^2),
                        V = -9,
                        Pvalue = -9,
                        Median = -9,
                        Upper = -9,
                        Lower = -9)

  # Fill the matrix
  row_indx = 0
  for(i in 1:(p-1)){
    for(j in 2:p){
      if(i<j){
        row_indx = row_indx + 1
        #Test[i,j]=1
        x = X[i,j,]
        y = Y[i,j,]
        w.xy = wilcox.test(x=x,y=y, paired=paired,conf.int = 0.95)
        Wilcoxon$from[row_indx] = var_names[i]
        Wilcoxon$to[row_indx] = var_names[j]
        Wilcoxon$V[row_indx] = w.xy$statistic
        Wilcoxon$Pvalue[row_indx] = w.xy$p.value
        Wilcoxon$Median[row_indx] = w.xy$estimate
        Wilcoxon$Upper[row_indx] = w.xy$conf.int[2]
        Wilcoxon$Lower[row_indx] = w.xy$conf.int[1]
      }

    }
  }

  # Drop unnecesary rows
  Wilcoxon <- Wilcoxon%>%filter(from!="FROM")

  # Plot
  Wilcoxon_median_plot <-
    Wilcoxon%>%
    ggplot(aes(x=from,y=Median))+
    geom_point(position=position_dodge(width=.75))+
    geom_errorbar(aes(ymin=Lower, ymax=Upper),
                  lwd=1, width=0, position=position_dodge(width=.75),
                  color="navy",alpha=0.7)+
    geom_hline(aes(yintercept=0),color="red")+
    facet_wrap(~to,scales = "free")+
    coord_flip()+
    labs(x="From Node", y="",color="",title="Median Difference and CI")

  return(list(Wilcoxon=Wilcoxon,
              Wilcoxon_plot = Wilcoxon_median_plot))

}
















