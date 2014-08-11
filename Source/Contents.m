% SOURCE
% Version 0.4.0 13-Feb-2014
%
% Files
%   AddCompartment                             - AddCompartment Add a compartment to a KroneckerBio model
%   AddInput                                   - AddInput Add an input species to a KroneckerBio model
%   AddLargeSizeReaction                       - AddReaction Add a reaction with an unlimited number of products to a
%   AddOutput                                  - AddOutput Add an output to a KroneckerBio model
%   AddParameter                               - AddParameter Add a rate parameter to a KroneckerBio model
%   AddReaction                                - AddReaction Add a reaction to a KroneckerBio model
%   AddSeed                                    - AddSeed Add a seed parameter to a KroneckerBio model
%   AddState                                   - AddState Add a state species to a KroneckerBio model
%   BestParameterExperiment                    - BestParameterExperiment Determine which experiments will most efficiently
%   BestTopologyExperiment                     - [best data] = BestTopologyExperiment(m, con, obj, objPrior, conPos, objPos, target, opts, pmyStart)
%   chi2cutoff                                 - Chi-square cut-off
%   chi2pvalue                                 - Chi-square p-value
%   completeInputFunction                      - This is a helper function that turns a bunch of input functions from
%   dealfield                                  - From a vector of objects or structs, deal the contents of a
%   dealfield2cell                             - Move all contents of a field to a cell vector
%   emptystruct                                - 
%   entropy                                    - Entropy of a discrete probability density.
%   Experiment                                 - Experiment constructs a KroneckerBio experimental conditions structure
%   FinalizeModel                              - FinalizeModel Update the mathematical components of the model to reflect
%   FiniteObjectiveGradient                    - FiniteObjectiveGradient Evaluate the gradient of a set of objective
%   FiniteObjectiveHessian                     - FiniteObjectiveHessian Compute the hessian of an objective function by
%   FitObjective                               - FitObjective Optimize the parameters of a model to minimize a set of
%   fixInputValueFunction                      - Standardize the input function as @(t,q)
%   goalParameterSpaceBelowTenPercent          - 
%   goalParameterSpaceErrorSumOfSqrt           - 
%   goalParameterSpaceSumOfLog                 - 
%   goalParameterSpaceSumOfLogAboveTenPercent  - 
%   GoodAbsTol                                 - GoodAbsTol Make a reasonable estimate as what the absolute tolerance on
%   Gzero                                      - Gzero Default structure for Kronecker Bio objective functions
%   inf2big                                    - inf2big converts infinities to just big numbers
%   infoeig                                    - Perform the eigendecomposition of an information matrix that is
%   infoinv                                    - Desparsify
%   InitializeModel                            - InitializeModel Create an empty mass action Kronecker model structure
%   intrnd                                     - Draw random integers
%   isoctave                                   - 
%   kroncol                                    - Fast kronecker product of two column vectors
%   linernd                                    - choses a random point along the line segment connecting two points
%   LoadModel                                  - LoadModel Load a model from files using various methods
%   LoadModelMassAction                        - LoadModelMassAction Load a model from the standard Kronecker model file
%   LoadModelSbmlAnalytic                      - LoadModelSbmlAnalytic Convert an SBML model into a pseudo-kronecker model
%   LoadModelSbmlMassAction                    - LoadModelSbmlMassAction Load a mass action model from an SBML file, a
%   lookup                                     - Map a vector onto a dictionary
%   MapParameterSpace                          - 
%   mergestruct                                - MERGESRUCT Merge one structure into another.
%   mvnbndrndgibbs                             - Sample from the bounded multivariate normal distribution
%   mvnpvalue                                  - Multivariant normal p-value
%   nan2zero                                   - nan2zero turns all the nans in matrix A into zeros
%   NetworkDistance                            - NetworkDistance Compute the distance that any two species in a network are
%   normbndrnd                                 - Draw random samples from the bounded normal distribution.
%   ObjectiveGradient                          - ObjectiveGradient Evaluate the gradient of a set of objective functions
%   ObjectiveHessian                           - 
%   ObjectiveInformation                       - ObjectiveInformation Compute the Fisher information matrix of a set of
%   ObjectiveLogLikelihood                     - ObjectiveLogLikelihood Evaluate the log likelihood of a set of 
%   objectiveLogNormalPriorOnKineticParameters - obj = objectiveNormalPriorOnKineticParameters(kbar, Vkbar, name)
%   objectiveLogNormalPriorOnSeedParameters    - obj = objectiveLogNormalPriorOnSeedParameters(m, sbar, Vsbar, name)
%   objectiveNormalPriorOnKineticParameters    - obj = objectiveNormalPriorOnKineticParameters(kbar, Vkbar, name)
%   objectiveOutputTrace                       - obj = objectiveTargetFunction(m, yInd, yTarget, ySd, name)
%   objectivePosteriorLikelihood               - 
%   ObjectiveProbability                       - ObjectiveProbability Evaluate the likelihood of a set of 
%   ObjectiveValue                             - ObjectiveValue Evaluate a set of objective functions
%   objectiveWeightedSumOfSquares              - obj = objectiveWeightedSumOfSquares(outputlist, timelist, sd, measurements, name)
%   objectiveWeightedSumOfSquaresFixedToData   - obj = objectiveWeightedSumOfSquaresFixedToData(m, outputlist, timelist, sd, measurements, name)
%   objectiveWeightedSumOfSquaresNonNeg        - obj = objectiveWeightedSumOfSquaresNonNeg(outputlist, timelist, sd, measurements, name)
%   OptimalAbsTol                              - 
%   ParameterExperimentElimination             - ParameterExperimentElimination removes experiments from consideration in
%   pastestruct                                - Paste an new structure over an old, tossing new fields
%   piecewiselinear                            - interpolates within a series of tvalues and
%   piecewisestep                              - returns the appropriate pvalues to a stepwise-defined
%   plotCurvatureExperiment                    - 
%   plotExperiment                             - 
%   plotSensitivityExperiment                  - 
%   pointvectorintercept                       - takes two points and two vectors eminating from them 
%   posdef                                     - Make symmetric matrix positive definite
%   RemoveCompartment                          - Find added instances of this name
%   RemoveInput                                - Find inputs by this name
%   RemoveOutput                               - Find added instances of this name
%   RemoveParameter                            - Find added instances of this name
%   RemoveReaction                             - 
%   RemoveSeed                                 - Find added instances of this name
%   RemoveState                                - Find states by this name
%   RenameModel                                - 
%   row                                        - Row vectorize a matrix
%   SampleParameterSpace                       - SampleParameterSpace Sample the parameter space according to the
%   SaveModel                                  - SaveModel writes the Kronecker mass action model to a Kronecker mass
%   SaveModelSbml                              - 
%   sdLinear                                   - 
%   simbio2Symbolic                            - simbio2Symbolic converts a Simbiology model into a structure of symbolic
%   Simulate                                   - Simulate Integrate the concentration of every species over a time
%   SimulateCurvature                          - SimulateCurvature Integrate the second-order sensitivities 
%   SimulateCurvatureSelect                    - SimulateCurvatureSelect Integrate the second-order sensitivities 
%   SimulateLna                                - SimulateLna Integrate the concentration of every species and its variance
%   SimulateMfk                                - Work-up
%   SimulateSelect                             - Simulate Integrate the concentration of every species over a time
%   SimulateSensitivity                        - SimulateSensitivity Integrate the sensitivities of every species with
%   SimulateSensitivitySelect                  - SimulateSensitivitySelect Integrate the sensitivities of every species with
%   spermute132                                - Permute the second and third axis of a three-dimensional array
%   spermute213                                - Permute the first and second axis of a three-dimensional array
%   subindex                                   - Subscript indexes
%   symbolic2MassAction                        - 
%   symbolic2PseudoKronecker                   - symbolic2PseudoKronecker converts a symbolic model into a pseudo-kronecker
%   symmat                                     - Symmetrize a matrix.
%   Tmatrix                                    - 
%   TopologyProbability                        - TopologyProbability Compute the relative probability that each member of
%   uncertaintySpectra                         - F is a Fisher information matrix or cell array of them; labels is the
%   uncertaintyVectors                         - 
%   Update                                     - 
%   Uzero                                      - Default structure for Kronecker Bio experimental conditions
%   vec                                        - Vectorize a matrix
%   weightedmean                               - Weighted mean
%   weightedvar                                - Weighted variance
%   zpvalue                                    - Upper-tail univariate normal p-value
