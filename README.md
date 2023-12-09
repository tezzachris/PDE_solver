# Multi-dimensional Black Scholes (BS) PDE Numerical Solution

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Aim: Compute numerical solution of classical BS PDE using monte carlo regression

The numerical scheme corresponds to:

Reference: A MONOTONE SCHEME FOR HIGH-DIMENSIONAL FULLY NONLINEAR PDES 
Journal: The Annals of Applied Probability (2015)
Authors: Wenjie Guo, Jianfeng Zhang and Jia Zhuo
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""



There are two scripts in the R language:

       growthIndianaExample.Rs   and
       ultrasoundExample.Rs

which carry out the illustrative data examples involving streamlined 
variational inference for higher level group specific curve models.

   The script growthIndianaExample.Rs fits two-level group specific curves, 
with a contrast curve extension, to data from a longitudinal study on 
adolescent somatic growth from the Indiana University School of Medicine, 
Indianapolis, Indiana, U.S.A. This analysis corresponds to Figure 2
and surrounding discussion.

   The script ultrasoundExample.Rs fits three-level group specific curves 
to data from an experiment involving quantitative ultrasound technology
from the Bioacoustics Research Laboratory, Department of Electrical 
and Computer Engineering, University of Illinois at Urbana-Champaign, 
Illinois, U.S.A. This analysis corresponds to Figure 3 and surrounding 
discussion.

   The scripts are supported by files containing MATLAB functions.

   The following table summarises all files:

-----------------------------------------------------------------------
   file name                            description 
-----------------------------------------------------------------------

SolveThreeLevelSparseLeastSquares.r     R function for carrying out the
                                        SolveThreeLevelSparseLeastSquares 
                                        algorithm. 
                                        
                                       

SolveTwoLevelSparseLeastSquares.r       R function for carrying out the
                                        SolveTwoLevelSparseLeastSquares     
                                        algorithm.

ZOSull.r                                R function for constructing a design 
                                        matrix consisting of O'Sullivan cubic 
                                        spline functions of an array of abscissae. 
                                        Typicially the array corresponds to either 
                                        observed values of a predictor or an 
                                        abscissa grid for plotting purposes. 
                          
growthIndiana.txt                       Data on adolescent somatic growth 
                                        obtained from a study of the mechanisms 
                                        of human hypertension development conducted 
                                        at the Indiana University School of Medicine, 
                                        Indianapolis, Indiana, U.S.A. The data are 
                                        restricted to a subset of 216 adolescents 
                                        in the original study who had at least 9 
                                        height measurements. There are a total of 
                                        4,123 height measurements taken approximately 
                                        every 6 months. The columns are as follows:
                                
                                         idnum:  identification numbers of the 216 
                                                 adolescents.
                                         height: height in centimeters.
                                         age:    age in years.
                                         male:   indicator of the adolescent being 
                                                 male:
                                                  1 = adolescent is male,
                                                  0 = adolescent is female.
                                         black:  indicator of the adolescent being 
                                                 black:
                                                  1 = adolescent is black,
                                                  0 = adolescent is not black.
            
                                        Fuller details on these data are in:

                                        Pratt, J.H., Jones, J.J., Miller, J.Z., Wagner, M.A. and 
                                        Fineberg, N.S. (1989). Racial differences in aldosterone 
                                        excretion and plasma aldosterone concentrations in children. 
                                        New England Journal of Medicine, 321, 1152-1157.

growthIndianaExample.Rs                 R script which fits two-level group specific 
                                        curves, with a contrast curve extension, to data 
                                        from a longitudinal study on adolescent somatic 
                                        growth from the Indiana University School of 
                                        Medicine, Indianapolis, Indiana, U.S.A. 

matrixSqrt.r                            R function for computing the positive square-root 
                                        of a square positive definite matrix.

threeLevGrpSpecCurvLogML.r              R function for computing the marginal log-likelihood
                                        corresponding to fitting three-level group specific
                                        curve models using streamlined mean field variational
                                        Bayes.
                                        
threeLevGrpSpecCurvMFVB.r               R function for fitting three-level group specific
                                        curve models using streamlined mean field variational
                                        Bayes.

twoLevGrpSpecCurvContrMFVB.r            R function for fitting two-level group specific
                                        curve models, with contrast function extension, 
                                        using streamlined mean field variational Bayes.

ultrasound.txt                          Subset of data from from an experiment involving quantitative 
                                        ultrasound technology. The subset corresponds to a measurements
                                        from the MS400 transducer and 10 tumors, with one tumor in
                                        each of 10 laboratory mice. Each curve corresponds to a 
                                        logarithmically transformed backscatter coefficient, with the 
                                        backscatter coefficient measurements in units  
                                        1/(centimeters*steradian). Each backscatter/frequency curve
                                        corresponds to one of 5 slices of the same tumor, where slices 
                                        correspond to probe locations. The slices are grouped according 
                                        to being from one of 10 tumors.

                                        Fuller details on these data are in:

                                        Wirtzfeld, L.A., Ghoshal, G., Rosado-Mendez, I.M., 
                                        Nam, K., Park, Y.,Pawlicki, A.D., Miller, R.J., 
                                        Simpson, D.G., Zagzebski, J.A., Oelze, M.I., 
                                        Hall, T.J. and O'Brien, W.D. (2015).
                                        Quantitative ultrasound comparison of MAT and 4T1 
                                        mammary tumors in mice and rates across multiple imaging 
                                        systems. Journal of Ultrasound Medicine, 34, 1373-1383.

ultrasoundExample.Rs                    R script which fits three-level group specific 
                                        curves to data from an experiment involving 
                                        quantitative ultrasound technology from the 
                                        Bioacoustics Research Laboratory, Department of 
                                        Electrical and Computer Engineering, University 
                                        of Illinois at Urbana-Champaign, Illinois, U.S.A. 

-------------------------------------------------------------------------------------


                                            ---o---O---o---

Come scegliere basis functions usate per approssimare la conditional expectation tramite regressione polinomiale
CheriditoGersey - COMPUTATION OF CONDITIONAL EXPECTATIONS WITH GUARANTEES





