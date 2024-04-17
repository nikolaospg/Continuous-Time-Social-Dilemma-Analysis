# Continuous-Time-Social-Dilemma-Analysis

In this project, I work with Continuous Time Social Dilemmas, as described [in](https://github.com/thanasiskehagias/MyPapers/blob/main/books/2022Kehagias.pdf). In this analysis, we use parametric non linear dynamical system models (action update models) to describe the behaviour of social dilemmas. The equations are polynomials and the state values express the degree of cooperation between players.
Through mathematical analysis, we can work with a game model (supposed to be known to us) and find the dynamical system equilibria, as opposed to Nash equilibria which is more frequent on Game Theory. In more complex dynamics we could also have a more complex steady state behaviour (attractors in general), which is something that is also studied through mathematical analysis, along with notions like the attraction basin.
However, in games with a higher complexity, the mathematical study of the steady state behaviour can be oftentimes very demanding (if not practically impossible). In this project, I propose a methodology which uses simulation and ML techniques to offer a computational way to study the steady state characteristics of a game. This can help us gain some valuable intuition and understanding of how the game works, before resorting (if necessary) to analytical methods. An abstract description of the methodology would be:

1. Compute ẋ for some defined states, create a phase plot (for Num. Players ≤ 3)
2. Simulate for some initial conditions, plot the results.
3. Simulate for a large number of initial conditions, carefully chosen - the final states are the Steady State Dataset
4. Use a Density Based Clustering algorithm to separate the points belonging in different attractors/equilibria.
5. Use a Classification algorithm to predict the attractor (output) of the steady state for some given initial state (input). This is an estimator of the basin of attraction for each initial condition.
6. Estimate the ”Dimensionality” of each attractor, using some proper technique (e.g. PCA).
7. Try to model, if possible, the attractors (regression analysis, curve modeling etc.)

In the report.pdf file I offer a more detailed analysis of the project. I also provide applications of the method through 2 experiments, one on a 2 Player and another on a 3 Player game (parameters chosen randomly).  
