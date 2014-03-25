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

