\name{HiCseg_linkC_R}
\alias{HiCseg_linkC_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Link between C and R
}
\description{
This function makes the link between C language and the R software. It consists 
in a two-dimensional segmentation approach.
}
\usage{
HiCseg_linkC_R(taille_mat, nb_rupt_max, distrib, mat_donnees, modele)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taille_mat}{
Size of the data matrix
}
  \item{nb_rupt_max}{
Maximal number of change-points
}
  \item{distrib}{
Distribution of the data: "B" is for Negative Binomial distribution, "P" is for the Poisson distribution and "G" is for the Gaussian distribution.
}
  \item{mat_donnees}{
Matrix of data
}
  \item{modele}{
Type of model: "D" for block-diagonal and "Dplus" for the extended block-diagonal model.
}
}

\value{
Contains three attributes : 
\item{tchap}{Contains the estimated change-points}
\item{J}{Values of the log-likelihood for different number of change-points up to some constants} 
\item{t_est_mat}{It gives the matrix of the estimated change-points for different number of change-points: in the first line when there is no change-point,
in the second line when there is one change-point, in the third line when there are two change-points....}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{ The method developped in this package is described in the 
paper "Two-dimensional segmentation for analyzing HiC data" by C. Levy-Leduc, M. Delattre, T. Mary-Huard and S. Robin, submitted to ECCB 2014.
}
\author{
Celine Levy-Leduc
}


\examples{

## The function is currently defined as
HiCseg_linkC_R <-
function(taille_mat,nb_rupt_max,distrib,mat_donnees,modele)
{
  K=nb_rupt_max^2
  
  tmp=.C("Fonction_HiC_R",as.integer(taille_mat),as.integer(nb_rupt_max),
         as.character(distrib),as.double(as.vector(mat_donnees)),
         tchap=as.integer(rep(0,nb_rupt_max)),J=as.double(rep(0.0,nb_rupt_max)),
         t_est=as.integer(rep(0,K)),as.character(modele))
  
  t_est_mat=matrix(tmp$t_est,ncol=nb_rupt_max,byrow=T)

  return(list(tchap=tmp$tchap,J=tmp$J,t_est_mat=t_est_mat))
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ segmentation }
\keyword{ two-dimensional }% __ONLY ONE__ keyword per line
