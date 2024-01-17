package yimei.jss.niching;

import ec.EvolutionState;
import ec.Individual;
import ec.Population;
import ec.coevolve.GroupedProblemForm;
import ec.coevolve.MultiPopCoevolutionaryEvaluator;
import ec.util.Parameter;
import org.apache.commons.lang3.ArrayUtils;
import smile.math.distance.EuclideanDistance;
import smile.math.distance.Metric;
import smile.regression.RBFNetwork;
import yimei.jss.ruleevaluation.MultipleRuleEvaluationModel;
import yimei.jss.ruleoptimisation.RuleOptimizationProblem;

import java.util.ArrayList;

import static yimei.jss.niching.surrogateClearing.indsCharLists;

/**
 * Created by dyska on 26/09/17.
 */
//a class can not extend from more than one class
public class surrogateClearingMultiPopCoevolutionaryEvaluatorV3 extends MultiPopCoevolutionaryEvaluator {
    public static final String P_RADIUS = "radius";
    public static final String P_CAPACITY = "capacity";

    //fzhang 2018.10.9 to get the pre-generation value
    public static final String P_PRE_GENERATIONS = "pre-generations";
    ArrayList<Double[]> fitnessesForModel = new ArrayList<>();

    int maxNumInstances  = 10000;

    /*ArrayList<Double> fitnessesForModelSubPop0 = new ArrayList<>();
    ArrayList<Double> fitnessesForModelSubPop1 = new ArrayList<>();*/

    ArrayList<double[][]> indsCharListsForModel = new ArrayList<>();
  /*  int[][] indsCharListsForModelSubPop0 = null; //the dictionary for looking for fitness for pop0
    int[][] indsCharListsForModelSubPop1 = null; //the dictionary for looking for fitness for pop1*/

    protected boolean clear = true;
    protected boolean nonIntermediatePop = true;

    protected double radius;
    protected int capacity;

    protected PhenoCharacterisation[] phenoCharacterisation;

    public double getRadius() {
        return radius;
    }

    public int getCapacity() {
        return capacity;
    }

    public PhenoCharacterisation[] getPhenoCharacterisation() {
        return phenoCharacterisation;
    }

    public PhenoCharacterisation getPhenoCharacterisation(int index) {
        return phenoCharacterisation[index];
    }

    public void setup(final EvolutionState state, final Parameter base) {
        super.setup(state, base);

        radius = state.parameters.getDoubleWithDefault(
                base.push(P_RADIUS), null, 0.0);
        capacity = state.parameters.getIntWithDefault(
                base.push(P_CAPACITY), null, 1);
        String filePath = state.parameters.getString(new Parameter("filePath"), null);
        //It's a little tricky to know whether we have 1 or 2 populations here, so we will assume
        //2 for the purpose of the phenoCharacterisation, and ignore the second object if only
        //1 is used
        phenoCharacterisation = new PhenoCharacterisation[2];
        if (filePath == null) {
            //dynamic simulation
            phenoCharacterisation[0] =
                    SequencingPhenoCharacterisation.defaultPhenoCharacterisation();
            phenoCharacterisation[1] =
                    RoutingPhenoCharacterisation.defaultPhenoCharacterisation();
        } else {
            //static simulation
            phenoCharacterisation[0] =
                    SequencingPhenoCharacterisation.defaultPhenoCharacterisation(filePath);
            phenoCharacterisation[1] =
                    RoutingPhenoCharacterisation.defaultPhenoCharacterisation(filePath);
        }
    }

    @Override
    public void evaluatePopulation(final EvolutionState state) {

     /*   if (clear) { //always true

            //1. remove the duplicated individuals
            Clearing.clearPopulation(state, radius, capacity,
                    phenoCharacterisation);
            //this is only to get the phenotypic characteristic for all individuals
            indsCharLists = surrogateClearing.phenotypicPopulation(state, phenoCharacterisation);
        }*/
         //super.evaluatePopulation(state, indsCharLists); //fzhang  all the evolution process is the same. for one individual, get a fitness value
        //2. full time evaluation, but do not need to evaluate the individuals that already evaluated

        //if it is a normal evaluation
        if (nonIntermediatePop) {
            //fzhang do not clearing individuals before full evaluation. In Jurgen's paper, the initialization population is cleared once and the intermediate is cleared every time.
          /*  ClearingFirst.clearPopulation(state, radius, capacity, phenoCharacterisation); // 1. clear it
            super.evaluatePopulation(state, clear); //2. evaluate it --- do not evaluate the cleared ones*/

         /* if(state.generation == 0){
              ClearingFirst.clearPopulation(state, radius, capacity, phenoCharacterisation);
              super.evaluatePopulation(state, clear); //2. evaluate it --- do not evaluate the cleared ones
          }
          else{
              super.evaluatePopulation(state);
          }*/

        //do not clear
            super.evaluatePopulation(state);

            indsCharLists = surrogateClearing.phenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

            ArrayList<Double[]> tempFitnessesForModel = ((MultipleRuleEvaluationModel)(((RuleOptimizationProblem)this.p_problem).getEvaluationModel())).getFitnessFixedSimulation();
            //ArrayList<Double[]> tempFitnessesForModel = multipleRuleEvaluationModel.getFitnessFixedSimulation();
            //4. update the surrogate model -- KNN: update the table for looking after evaluation on each generation
            //on one hand, filter new generated phenotypic characterisation; one the other hand, put all the PC from previous generations together
          /*  int[][] tempIndsCharListsForModelSubPop0 = null;
            int[][] tempIndsCharListsForModelSubPop1 = null;*/
           ArrayList<double[][]> tempIndsCharListsForModel = new ArrayList<>();
            for (int subpop = 0; subpop < state.population.subpops.length; subpop++) {
                double[][] indsCharListsForModelSubPop = indsCharLists[subpop];
                int removeIdx = 0;

                for (int i = 0; i < tempFitnessesForModel.get(subpop).length; i++) {
                    if (tempFitnessesForModel.get(subpop)[i] == Double.MAX_VALUE) {
                        indsCharListsForModelSubPop = ArrayUtils.remove(indsCharListsForModelSubPop, i - removeIdx);
                        Double[] temp = ArrayUtils.remove(tempFitnessesForModel.get(subpop), i - removeIdx);
                        tempFitnessesForModel.set(subpop, temp);
                        removeIdx++;
                    }
                }

                if(fitnessesForModel.size() == 0 || fitnessesForModel.size() == 1){
                    fitnessesForModel.add(subpop, tempFitnessesForModel.get(subpop));
                }
                else{
                    Double[] combineFitness = ArrayUtils.addAll(fitnessesForModel.get(subpop), tempFitnessesForModel.get(subpop));
                    fitnessesForModel.set(subpop, combineFitness);
                }

                if(subpop == 1){
                    ((MultipleRuleEvaluationModel)(((RuleOptimizationProblem)this.p_problem).getEvaluationModel())).clearingFitnessFixedSimulation();
                }

                //fzhang 2019.9.15 solve the problem that the arraylist is too large
                if(fitnessesForModel.get(subpop).length > maxNumInstances)
                {
                    int copyElement = fitnessesForModel.get(subpop).length - maxNumInstances;
                    Double[] tempForDeleteElement = new Double[maxNumInstances];
                    System.arraycopy(fitnessesForModel.get(subpop), copyElement, tempForDeleteElement, 0, maxNumInstances);
                    fitnessesForModel.set(subpop, tempForDeleteElement);
                }

                tempIndsCharListsForModel.add(subpop, indsCharListsForModelSubPop);//the remaining phenotypic characterisation after removing
                if(indsCharListsForModel.size() == 0 || indsCharListsForModel.size() == 1){
                    indsCharListsForModel.add(subpop, tempIndsCharListsForModel.get(subpop));
                }
                else{
                    double[][] tempCombined = ArrayUtils.addAll(indsCharListsForModel.get(subpop), tempIndsCharListsForModel.get(subpop));
                    indsCharListsForModel.set(subpop, tempCombined);
                }

                //fzhang 2019.9.15 solve the problem that the arraylist is too large
                if(indsCharListsForModel.get(subpop).length > maxNumInstances)
                {
                    int copyElement = indsCharListsForModel.get(subpop).length - maxNumInstances;
                    double[][] tempForDeleteElement = new double[maxNumInstances][];
                    System.arraycopy(indsCharListsForModel.get(subpop), copyElement, tempForDeleteElement, 0, maxNumInstances);
                    indsCharListsForModel.set(subpop, tempForDeleteElement);
                }
            }

            nonIntermediatePop = false;
        } else {
            ClearingFirst.clearPopulation(state, radius, capacity, phenoCharacterisation);//the individuals do not have fitness or they are evaluated yet
                                                                                          //but has a fitness there
            super.evaluatePopulation(state, indsCharListsForModel);
            nonIntermediatePop = true;
            }
        }


    // individuals to evaluate together
    Individual[] inds = null;
    // which individual should have its fitness updated as a result
    boolean[] updates = null;

    public void performCoevolutionaryEvaluation( final EvolutionState state,
                                                 final Population population,
                                                 final GroupedProblemForm prob,
                                                 final Boolean clear) {
        int evaluations = 0;

        inds = new Individual[population.subpops.length];
        updates = new boolean[population.subpops.length];//check if we need to evaluate the individuals, because it is so expensive, should be careful.

        //==========================here skip step 1: load the selectionMethod for current generation==================================
        //numCurrent = 0
        if (numCurrent > 0) {
            for (int i = 0; i < selectionMethodCurrent.length; i++) {
                selectionMethodCurrent[i].prepareToProduce(state, i, 0);  //A default version of prepareToProduce which does nothing.
            }
        }

        //==========================here skip: step 2: load the selectionMethod for previous generation==================================
        if (numPrev > 0) {
            for (int i = 0; i < selectionMethodPrev.length; i++) {
                // do a hack here
                Population currentPopulation = state.population; //state.population is current population.
                state.population = previousPopulation;
                selectionMethodPrev[i].prepareToProduce(state, i, 0);  //A default version of prepareToProduce which does nothing.  here state is important, it is the state
                //for previous population
                state.population = currentPopulation; //currentPopulaiton is a temp parameter to save current population
            }
        }

        //step 3: build subPopulaiton, subpops[0] and subpops[1]
        // build subpopulation array to pass in each time
        int[] subpops = new int[state.population.subpops.length];
        //System.out.println(subpops.length);
        for(int j = 0; j < subpops.length; j++) {
            subpops[j] = j;
        }

        //System.out.println(prob);  //yimei.jss.ruleoptimisation.RuleCoevolutionProblem@5ae63ade
        //here skip: step 3: setup the shuffle: here num-shuffled = 0  here, we do not use it.
        if (numShuffled > 0) {
            int[/*numShuffled*/][/*subpop*/][/*shuffledIndividualIndexes*/] ordering = null;
            // build shuffled orderings
            ordering = new int[numShuffled][state.population.subpops.length][state.population.subpops[0].individuals.length]; // if num-shuffled =1  [1][2][512]
            for(int c = 0; c < numShuffled; c++)
                for(int m = 0; m < state.population.subpops.length; m++)
                {
                    for(int i = 0; i < state.population.subpops[0].individuals.length; i++)
                        ordering[c][m][i] = i; //ordering[0][0][0] = 0, ordering[0][0][1] = 1, ordering[0][0][2] = 2, ordering[0][0][3] = 3
                    if (m != 0)
                        shuffle(state, ordering[c][m]); //ordering = new int[numShuffled][state.population.subpops.length]
                }

            // for each individual
            for (int i = 0; i < state.population.subpops[0].individuals.length; i++) {
                for (int k = 0; k < numShuffled; k++) {
                    for (int ind = 0; ind < inds.length; ind++) {
                        inds[ind] = state.population.subpops[ind].individuals[ordering[k][ind][i]];
                        updates[ind] = true;
                    }
                    prob.evaluate(state, inds, updates, false, subpops, 0);  ////yimei.jss.ruleoptimisation.RuleCoevolutionProblem@5ae63ade
                    evaluations++;
                }
            }
        }

        //==========================useful and important part=======================================
        //step 4: for each subpopulation, j means subPopulation   2*512*4*2= 8192  cost
        for (int j = 0; j < state.population.subpops.length; j++) //0,1
        {
            if (!shouldEvaluateSubpop(state, j, 0)) continue;  // don't evaluate this subpopulation

            // for each individual
            for (int i = 0; i < state.population.subpops[j].individuals.length; i++) //512
            {
                Individual individual = state.population.subpops[j].individuals[i];

                if(individual.evaluated == false){
                    // Test against all the elites
                    for (int k = 0; k < eliteIndividuals[j].length; k++) { //2
                        for (int ind = 0; ind < inds.length; ind++) { //2
                            if (ind == j) {   //j = 0, 1  (ind j) ---> (0 0) or (1 1) that is to say, this is the subpopulation1
                                inds[ind] = individual; //inds[0] = individual = state.population.subpops[0].individuals[0];
                                //inds[1] = individual = state.population.subpops[1].individuals[1];
                                //the individuals to evaluate together
                                updates[ind] = true;   // updates[0] = true    updates[1] = true   evaluate
                            }
                            else {  // this is subpopulation2
                                inds[ind] = eliteIndividuals[ind][k];   // (ind j) ---> (0 1) or (1 0)
                                //inds[1] = eliteIndividuals[1][*]   inds[0] = eliteIndividuals[0][*]
                                updates[ind] = false;  // do not evaluate
                            }
                        }

                        prob.evaluate(state,inds,updates, false, subpops, 0);
                        evaluations++;
                    }
                }else
                    {
                        individual.fitness.trials = new ArrayList();
                        individual.fitness.trials.add(individual.fitness.fitness());
                }

                //here, skip this part: test against random selected individuals of the current population
                for(int k = 0; k < numCurrent; k++) //0  skip this part
                {
                    for(int ind = 0; ind < inds.length; ind++) //2
                    {
                        if (ind == j) { inds[ind] = individual; updates[ind] = true; }
                        else { inds[ind] = produceCurrent(ind, state, 0); updates[ind] = true; }
                    }
                    prob.evaluate(state,inds,updates, false, subpops, 0);
                    evaluations++;
                }

                // here, skip this part. Test against random individuals of previous population
                for(int k = 0; k < numPrev; k++)  // 0  skip this part
                {
                    for(int ind = 0; ind < inds.length; ind++)
                    {
                        if (ind == j) { inds[ind] = individual; updates[ind] = true; }
                        else { inds[ind] = producePrevious(ind, state, 0); updates[ind] = false; }
                    }
                    prob.evaluate(state,inds,updates, false, subpops, 0);
                    evaluations++;
                }
            }
        }
        //============================================================================================================================

        //here, skip this part
        // now shut down the selection methods
        if (numCurrent > 0)
            for( int i = 0 ; i < selectionMethodCurrent.length; i++)
                selectionMethodCurrent[i].finishProducing( state, i, 0 );  //A default version of finishProducing, which does nothing.

        if (numPrev > 0)
            for( int i = 0 ; i < selectionMethodPrev.length ; i++ )
            {
                // do a hack here
                Population currentPopulation = state.population;
                state.population = previousPopulation;
                selectionMethodPrev[i].finishProducing( state, i, 0 );
                state.population = currentPopulation;
            }

        state.output.message("Evaluations: " + evaluations);
    }

//==========================================KNN========================================
/*    public void performCoevolutionaryEvaluation( final EvolutionState state,
                                                 final Population population,
                                                 final GroupedProblemForm prob,
                                                 ArrayList<int[][]> indsCharListsForModelSubPop) {
        int[][][] indsCharListsIntermediatePop = surrogateClearing.phenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
            {
                Individual individual = state.population.subpops[sub].individuals[i];
                if(indsCharListsForModelSubPop.get(sub).length != 0){
                    //KNN
                    //===============================start==============================================
                    double dMin  = Double.MAX_VALUE;
                    int index = 0;
                    if(individual.evaluated == false){
                        int[] pcIntermediate = indsCharListsIntermediatePop[sub][i];
                        //calculate the fitness based on surrogate model
                        for(int pc = 0; pc < indsCharListsForModelSubPop.get(sub).length; pc++){
                            int[] pcModel = indsCharListsForModelSubPop.get(sub)[pc];
                            double d  = PhenoCharacterisation.distance(pcIntermediate, pcModel);
                            if(d == 0){
                                index = pc;
                                break;
                            }

                            if (d < dMin){
                                dMin = d;
                                index = pc;
                            }
                        }
                        ((Clearable)individual.fitness).surrogateFitness(fitnessesForModel.get(sub)[index]);
                        //individual.evaluated = true;
                        //individual.fitness.trials = new ArrayList();
                        //individual.fitness.trials.add(individual.fitness.fitness());
                    }
                }
                else{
                    ((Clearable)individual.fitness).surrogateFitness(Double.MAX_VALUE);
                }
            }
        }
    }*/

//===================SVM=======================
   /* public void performCoevolutionaryEvaluation( final EvolutionState state,
                                                 final Population population,
                                                 final GroupedProblemForm prob,
                                                 ArrayList<double[][]> indsCharListsForModel) {

        double[][][] indsCharListsIntermediatePop = surrogateClearing.phenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            //if there is training data, train the model. Otherwise, do not need to train model, just assign them the same fitness (Double.MAX_VALUE)
            if (indsCharListsForModel.get(sub).length != 0) {
                //===============================start==============================================
            *//*    MercerKernel<double[]> mercerKernel = new MercerKernel<double[]>() {
                    @Override
                    public double k(double[] x, double[] y) {
                        return 0;
                    }
                };
*//*
              *//*  GaussianKernel<double[]> mercerKernel = new GaussianKernel<double[]>() {
                    @Override
                    public double k(double[] x, double[] y) {
                        return 0;
                    }
                };*//*

                GaussianKernel mercerKernel = new GaussianKernel(100);

                SVR.Trainer<double[]> trainer = new SVR.Trainer<double[]>(mercerKernel, 0.001, 1); //change the C from 1 to 10, no effect on fitness
                double[] fitnessForModel = ArrayUtils.toPrimitive(fitnessesForModel.get(sub));
                SVR<double[]> network = trainer.train(indsCharListsForModel.get(sub), fitnessForModel);

                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if (individual.evaluated == false) {
                        double[] pcIntermediate = indsCharListsIntermediatePop[sub][i];
                        //calculate the fitness based on surrogate model
                        double estimatedFitness = network.predict(pcIntermediate);
                        ((Clearable) individual.fitness).surrogateFitness(estimatedFitness);
                    }
                }
            } else {// do not use the surrogate, in this way, equals to use the original way to get offsprings
                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if(individual.evaluated == false){
                        ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE - 1000); //just remove the duplicated ones
                    }
                }
            }
        }
    }*/


  //=====================================RBF==========================================2019.9.15
    public void performCoevolutionaryEvaluation( final EvolutionState state,
                                                 final Population population,
                                                 final GroupedProblemForm prob,
                                                 ArrayList<double[][]> indsCharListsForModel) {

        double[][][] indsCharListsIntermediatePop = surrogateClearing.phenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            //if there is training data, train the model. Otherwise, do not need to train model, just assign them the same fitness (Double.MAX_VALUE)
            int maxNeighbour = (int) Math.round(0.01 * 10 * (indsCharListsForModel.get(sub).length - 10)); //insure the RBF can work
            if (indsCharListsForModel.get(sub).length != 0 && maxNeighbour >=1) {
                //===============================start==============================================
                Metric<double[]> metric = new EuclideanDistance();
                RBFNetwork.Trainer<double[]> trainer = new RBFNetwork.Trainer<double[]>(metric);
                double[] fitnessForModel = ArrayUtils.toPrimitive(fitnessesForModel.get(sub));
                RBFNetwork<double[]> network = trainer.train(indsCharListsForModel.get(sub), fitnessForModel);

                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if (individual.evaluated == false) {
                        double[] pcIntermediate = indsCharListsIntermediatePop[sub][i];
                        //calculate the fitness based on surrogate model
                        double estimatedFitness = network.predict(pcIntermediate);
                        if(estimatedFitness <= 0){
                            ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE);
                        }
                        else{
                            ((Clearable) individual.fitness).surrogateFitness(estimatedFitness);
                        }
                    }
                }
            } else {
                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if(individual.evaluated == false){
                        ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE - 1000); //just remove the duplicated ones
                    }
                }
            }
        }
    }

//2019.9.16
    //=======================================Gaussian Process Regression==============================================
/*    public void performCoevolutionaryEvaluation( final EvolutionState state,
                                                 final Population population,
                                                 final GroupedProblemForm prob,
                                                 ArrayList<double[][]> indsCharListsForModel) {

        double[][][] indsCharListsIntermediatePop = surrogateClearing.phenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            //if there is training data, train the model. Otherwise, do not need to train model, just assign them the same fitness (Double.MAX_VALUE)
            if (indsCharListsForModel.get(sub).length != 0) {
                //===============================start==============================================
                MercerKernel<double[]> mercerKernel = new MercerKernel<double[]>() {
                    @Override
                    public double k(double[] x, double[] y) {
                        return 0;
                    }
                };
                GaussianProcessRegression.Trainer<double[]> trainer = new GaussianProcessRegression.Trainer<double[]>(mercerKernel,0.1);
                double[] fitnessForModel = ArrayUtils.toPrimitive(fitnessesForModel.get(sub));
                GaussianProcessRegression<double[]> network = trainer.train(indsCharListsForModel.get(sub), fitnessForModel);

                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if (individual.evaluated == false) {
                        double[] pcIntermediate = indsCharListsIntermediatePop[sub][i];
                        //calculate the fitness based on surrogate model
                        double estimatedFitness = network.predict(pcIntermediate);
                        ((Clearable) individual.fitness).surrogateFitness(estimatedFitness);
                    }
                }
            } else {
                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if(individual.evaluated == false){
                        ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE - 1000); //just remove the duplicated ones
                    }
                }
            }
        }
    }*/


//LASSO
 /*   public void performCoevolutionaryEvaluation( final EvolutionState state,
                                                 final Population population,
                                                 final GroupedProblemForm prob,
                                                 ArrayList<double[][]> indsCharListsForModel) {

        double[][][] indsCharListsIntermediatePop = surrogateClearing.phenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            //if there is training data, train the model. Otherwise, do not need to train model, just assign them the same fitness (Double.MAX_VALUE)
            if (indsCharListsForModel.get(sub).length != 0) {
                //===============================start==============================================
                LASSO.Trainer trainer = new LASSO.Trainer(1);
                double[] fitnessForModel = ArrayUtils.toPrimitive(fitnessesForModel.get(sub));
                LASSO network = trainer.train(indsCharListsForModel.get(sub), fitnessForModel);

                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if (individual.evaluated == false) {
                        double[] pcIntermediate = indsCharListsIntermediatePop[sub][i];
                        //calculate the fitness based on surrogate model
                        double estimatedFitness = network.predict(pcIntermediate);
                        if(estimatedFitness <= 0){
                            ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE);
                        }
                        else{
                            ((Clearable) individual.fitness).surrogateFitness(estimatedFitness);
                        }
                    }
                }
            } else {
                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if(individual.evaluated == false){
                        ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE - 1000); //just remove the duplicated ones
                    }
                }
            }
        }
    }*/

    //NN
   /* public void performCoevolutionaryEvaluation( final EvolutionState state,
                                                 final Population population,
                                                 final GroupedProblemForm prob,
                                                 ArrayList<double[][]> indsCharListsForModel) {

        double[][][] indsCharListsIntermediatePop = surrogateClearing.phenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            //if there is training data, train the model. Otherwise, do not need to train model, just assign them the same fitness (Double.MAX_VALUE)
            if (indsCharListsForModel.get(sub).length != 0) {
                //===============================start==============================================

                int numInput = indsCharListsForModel.get(sub)[0].length;
                NeuralNetwork.Trainer trainer = new NeuralNetwork.Trainer(numInput,200, 200, 200,1);//the number of input, ...the number of nodes in each layer..., the nunmber of output
                double[] fitnessForModel = ArrayUtils.toPrimitive(fitnessesForModel.get(sub));
                NeuralNetwork network = trainer.train(indsCharListsForModel.get(sub), fitnessForModel);

                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if (individual.evaluated == false) {
                        double[] pcIntermediate = indsCharListsIntermediatePop[sub][i];
                        //calculate the fitness based on surrogate model
                        double estimatedFitness = network.predict(pcIntermediate);
                        ((Clearable) individual.fitness).surrogateFitness(estimatedFitness);
                    }
                }
            } else {
                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if(individual.evaluated == false){
                        ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE - 1000); //just remove the duplicated ones
                    }
                }
            }
        }
    }*/

    public void setClear(boolean clear) {
        this.clear = clear;
    }
}
