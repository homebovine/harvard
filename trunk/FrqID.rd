\name{FrqID}
\alias{FrqID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function to implement Xu et al. (2010)'s method to fit semi-parametric regression models to semi-competing risks data.
}
\description{
The function to implement Xu et al. (2010)'s method to fit semi-parametric regression models to semi-competing risks data. The function hands the semi-competing risk data where the  covariates are the same for both the non-terminated and terminated event time. The function also provides the information matrix for the finite dimensional paramters.  
}
\usage{
FrqID(survData,  startValues, stheta, wtheta, hessian, mitr, tol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survData}{
  The data frame containing semi-competing risks data with covariate matrices from \code{n} subjects. See *Examples*. 
}
 
 
  \item{startValues}{
	a vector of length 3 * p where p is the dimensions of the covariates. The first, second, and third  p elements are the starting values for \eqn{\beta_1}, \eqn{\beta_2}, \eqn{\beta_3}, respectively. 
}	
\item{stheta}{
	a small interval for theta
}
\item{wtheta}{
	a wider interval for theta
}
	\item{hessian}{If true, the function return the information matrix, whose inverse is the estimation variance-covariance matrix. 
}	
  \item{miter}{ The  inner iterative estimation  stops  when the iteration time reaches miter
}
\item{tol}{
	The inner estimation procedure stops when the maximal absolute difference between the current and the last finite dimentional estimators smaller than tol. 
}
}
\details{
In order to get a stablize the estimation procedure, we use the EM algorithm for the parameter estimations. The detailed procedures are describled in Fei Jiang (2013). To get a good estimator, we suggest to normlize the continous covariates before applying the function. 
}


\value{
\code{thbb} a vector containing the estimators for \eqn{\theta}, and $\eqn{\beta_1}$, $\eqn{\beta_2}$, $\eqn{\beta_3}$.  \cr
\code{vv1, vv2, vv3} matricies for estimated baseline hazard functions. Each matrix contains two columns, while the second column contains the event time.\cr
\code{hessian} when hessian is TRUE, return a  hessian matrix at estimated values.   
}
\references{
Jinfen Xu and John D. Kalbfleisch and Beechoo Tai (2010).
Statistical Analysis of Illness–Death Processes and Semicompeting Risks Data,  Biometrics
}
\referneces{
Fei Jiang (2013).\cr
An EM Algorithm for Semi-competing Risk data, Technique report.
}
\author{
Kyu Ha Lee and Sebastien Haneuse\cr
Maintainer: Kyu Ha Lee <klee@hsph.harvard.edu>
}
\note{
The pos	# prior parameter for J3
c_lam1 	<- 1	# prior parameter for MVN-ICAR specification
c_lam2 	<- 1	
c_lam3 	<- 1
psi	<- 0.7	# prior parameters for 1/theta
omega	<- 0.7	

nGam_save <- 5  # the number of subjects whose gamma parameters are being saved

hyperParams <- c(a1, b1, a2, b2, a3, b3, alpha1, alpha2, alpha3, c_lam1, c_lam2, c_lam3, 
psi, omega)

#########################
# setting starting values
  	  	   
s1_max			<- max(scrData$time1[scrData$event1 == 1])
s2_max			   <- max(scrData$time2[scrData$event1 == 0 & scrData$event2 == 1])
s3_max			      <- max(scrData$time2[scrData$event1 == 1 & scrData$event2 == 1])

beta1			      	 <- c(0.1, (unique(scrData$time2[scrData$event1 == 0 & scrData$event2 == 1])))
s_propBI2			 <- unique(s_propBI2[s_propBI2 < s2_max])
s_propBI3			 <- floor(sort(unique(scrData$time2[scrData$event1 == 1 & scrData$event2 == 1])))
s_propBI3			 <- unique(s_propBI3[s_propBI3 < s3_max])

num_s_propBI1			 <- length(s_propBI1)
num_s_propBI2			 <- length(s_propBI2)
num_s_propBI3			 <- length(s_propBI3)
J1_max 				    <- 50
J2_max 				       <- 50
J3_max 				       	  <- 50
time_lambda1				  <- 1:36
time_lambda2				  <- 1:36
time_lambda3				  <- 1:36
nTime_lambda1				  <- length(time_lambda1)
nTime_lambda2				  <- length(time_lambda2)
nTime_lambda3				  <- length(time_lambda3)

mhProp_theta_var <- 0.05

mcmcParams <- c(C1, C2, C3, delPert1, delPert2, delPert3, num_s_propBI1, num_s_propBI2, 
num_s_propBI3, J1_max, J2_max, J3_max, s1_max, s2_max, s3_max, nTime_lambda1, 
nTime_lambda2, nTime_lambda3, s_propBI1, s_propBI2, s_propBI3, time_lambda1, 
time_lambda2, time_lambda3, mhProp_theta_var)

##################################################
# number of chains

numReps	    = 2e6
thin	      = 20
burninPerc    = 0.5
path1	      	= "outcome/"
type 		  = "semi-parametric"
model 		    = "Markov"
nChain		      = 2

# fitting Bayesian semi-parametric regression model to semi-competing risks data	
# In """""")))))))])])))])])))])])]))}})})}}}})}}}}})}})})})})})}}})}}}}}})}}}}})}})}}}}})})}})}})}})})})}}}}})}}}}})}}}})}})}}})})}}}})}}}}})}}}}}}}}}}}}}}})}})}}}})}}}}}}}}}}}}}}}}}}}}}})})}}}})}}})})}}})}}})})}}})}}}}"}""}"}}}}"}""}"}}}}}}}}}}}}})}}}}}}}}})}}}}}}"""")}})}}}}