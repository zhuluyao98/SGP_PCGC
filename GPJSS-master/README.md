# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* This is a package for algorithms for automatic rule design in Job Shop Scheduling (JSS) using Genetic Programming (GP). Written by Yi Mei.
* The package is based on the Java ECJ package, which is available from https://cs.gmu.edu/~eclab/projects/ecj/.
* Version 1.0.0

### How do I get set up? ###

1. Download the source files in the `src/` folder and the dependencies in the `libraries/` folder.
2. Create a Java Project using `src/` and `libraries/` (the repository is a IDEA IntelliJ project by default).
3. The ECJ functions are located at `src/ec/`, and the functions for JSS and GP are in `src/yimei/jss/`. Now you are ready to run a number of different algorithms.
4. Before starting, it is highly recommended to get yourself familiar with the ECJ package, especially the GP part. You can start from the four tutorials located at `src/ec/app/tutorialx` (x = 1, ..., 4). Turorial 4 is about GP for symbolic regression, which is very useful for understanding this project. A more thorough manual can be found in https://cs.gmu.edu/~eclab/projects/ecj/docs/manual/manual.pdf.

### Project structure ###

The main project is located in `/src/yimei/jss/`. It contains the following packages:

* `algorithm/` contains a number of algorithms to run. They are good entry points to start with.
* `gp/` contains the core classes for GP for evolving dispatching rules.
* `jobshop/` contains the core classes for representing a job shop.
* `niching/` contains classes for niching-based GP, i.e. using the Clearing method.
* `rule/` contains the core classes for representing dispatching rules.
* `ruleanalysis/` contains the classes for rule analysis, e.g. reading rules from the ECJ result file, testing rules, calculating the program length, depth, number of unique terminals, etc.
* `ruleevalulation/` contains the evaluation models for dispatching rules, such as discrete event simulation, static job shop instances, etc.
* `ruleoptimisation/` contains the classes for optimisation dispatching rules, e.g. RuleOptimisationProblem.
* `simulation/` contains the classes for discrete event simulation for dynamic job shop scheduling.
* `surrogate/` contains the classes for surrogate models and evaluators.

### Running experiments ###

*Example (SGP-PCGP):**

L. Zhu, F. Zhang, X. Zhu, K. Chen*, and M. Zhang, “Phenotype and Genotype Based Sample Aware Surrogate-Assisted Genetic Programming in Dynamic Flexible Job Shop Scheduling,” IEEE Computational Intelligence Magazine.  (under review)

 Run the params file `src/yimei/algorithms/multitreeModelSurrogateSample/multitreegp-dynamic.params`

### Who do I talk to? ###

* Email: zhuluyao58@163.com
