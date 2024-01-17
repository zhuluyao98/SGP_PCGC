package yimei.jss.algorithm.multitreeModelSurrogateSample;

import ec.Individual;
import ec.Population;
import ec.util.Checkpoint;
import ec.util.Parameter;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.Pair;
import yimei.jss.gp.GPRuleEvolutionState;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.ruleoptimisation.RuleOptimizationProblem;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static yimei.jss.gp.GPRun.out_dir;


//fzhang 2019.9.11 change the steps of evolutionary process --- intermediate population
public class PhenotypicGPRuleEvolutionState extends GPRuleEvolutionState {
    public static final String P_REPLICATIONS = "num-Rep";
    public int numRep;

    public static final String P_ELITES = "num-elites";
    public static int numElites;

    ArrayList[] fitnessDiffGens = new ArrayList[6];
    public static Boolean realEvaluationIntermediate = false;
    public static ArrayList<Integer> selectedRealUnselectedEstimated = new ArrayList<>();
    public static ArrayList<Integer> selectedRealUnselectedEstimatedTop30Percent = new ArrayList<>();

    ArrayList<Pair<Double, Double>> distanceFitness;
    ArrayList<Double> criticalBetween = new ArrayList<>();
    public static ArrayList<Integer> indexEvaluated = new ArrayList<>();

    public int evolve() {
        if (generation > 0)
            output.message("Generation " + generation);

        statistics.preEvaluationStatistics(this);

        evaluator.evaluatePopulation(this);  //// here, after this we evaluate the population
        //((MultiPopCoevolutionaryClearingEvaluator)evaluator).setClear(false);

        //2021.6.17 keep the fitness in the whole pop to investigate the distribution of the inds
        //we save the fitness of inds at generation 0, 10, 20, 30, 40, 50
/*        if (this.generation == 0 || this.generation == 10 || this.generation == 20 || this.generation == 30
                || this.generation == 40 || this.generation == 50) {
            ArrayList tempFitDiffGens = new ArrayList();

            ArrayList<int[]> pcCurrentGens = new ArrayList<>();
            int[][] pcCurrentTemp = phenotypicForSurrogate.muchBetterPhenotypicPopulation(this, surrogateClearingMultitreeEvaluatorV1.phenoCharacterisation);

            for (int ind = 0; ind < this.population.subpops[0].individuals.length; ind++) {
                if (this.population.subpops[0].individuals[ind].fitness.fitness() != Double.POSITIVE_INFINITY &
                        this.population.subpops[0].individuals[ind].fitness.fitness() != Double.MAX_VALUE) {
                    tempFitDiffGens.add(this.population.subpops[0].individuals[ind].fitness.fitness());

                    pcCurrentGens.add(pcCurrentTemp[ind]);
                }
            }
            fitnessDiffGens[this.generation / 10] = tempFitDiffGens;


//            pcDiffGens[this.generation/10] = tempPCDiffGens;
            writePCDiffGens(pcCurrentGens);
        }*/
        //================================end==================================

        statistics.postEvaluationStatistics(this);

        // SHOULD WE QUIT?
        if (evaluator.runComplete(this) && quitOnRunComplete) {
            output.message("Found Ideal Individual");
            return R_SUCCESS;
        }
        // SHOULD WE QUIT?
/*        if (generation == numGenerations - 1) {

            generation++; // in this way, the last generation value will be printed properly.  fzhang 28.3.2018
//            writeFitnessDiffGens(fitnessDiffGens);
//            writeCorrealtionSurrogateAccuracyGens(surrogateClearingMultitreeEvaluatorV1.correlationSurrogateAccuracy);
//            writeDistanceRankSurrogateAccuracyGens(surrogateClearingMultitreeEvaluatorV1.distanceRankSurrogateAccuracy);
//            writeCorrealtionSurrogateAccuracyGensWhole(surrogateClearingMultitreeEvaluatorV1.correlationSurrogateAccuracyWhole);
//            writeDistanceRankSurrogateAccuracyGensWhole(surrogateClearingMultitreeEvaluatorV1.distanceRankSurrogateAccuracyWhole);

//            writeSelectedRealUnselectedEstimated(selectedRealUnselectedEstimated);
//            writeSelectedRealUnselectedEstimatedtop30Percent(selectedRealUnselectedEstimatedTop30Percent);
//            writeAverageFitnessGensToFile(averageFitnessGens);
            return R_FAILURE;
        }*/

        if (totalTime >= endTime) {
            generation++; // in this way, the last generation value will be printed properly.  fzhang 28.3.2018
            //writeSumTrainTime(sumTrainTime);
            return R_FAILURE;
        }

        // PRE-BREEDING EXCHANGING
        statistics.prePreBreedingExchangeStatistics(this);
        population = exchanger.preBreedingExchangePopulation(this);  /** Simply returns state.population. */
        statistics.postPreBreedingExchangeStatistics(this);

        String exchangerWantsToShutdown = exchanger.runComplete(this);  /** Always returns null */
        if (exchangerWantsToShutdown != null) {
            output.message(exchangerWantsToShutdown);
            return R_SUCCESS;
        }

        // BREEDING
        statistics.preBreedingStatistics(this);
        population = preselection();//actually, breed intermediate population and then resize the population to the normal size.  pre-selection

        //population = breeder.breedPopulation(this); //!!!!!!   return newpop;  if it is NSGA-II, the population here is 2N

        // POST-BREEDING EXCHANGING
        statistics.postBreedingStatistics(this);   //position 1  here, a new pop has been generated.

        // POST-BREEDING EXCHANGING
        statistics.prePostBreedingExchangeStatistics(this);
        population = exchanger.postBreedingExchangePopulation(this);   /** Simply returns state.population. */
        statistics.postPostBreedingExchangeStatistics(this);  //position 2

        // Generate new instances if needed
        RuleOptimizationProblem problem = (RuleOptimizationProblem) evaluator.p_problem;
        if (problem.getEvaluationModel().isRotatable()) {
            problem.rotateEvaluationModel();
        }

        // INCREMENT GENERATION AND CHECKPOINT
        generation++;
        if (checkpoint && generation % checkpointModulo == 0) {
            output.message("Checkpointing");
            statistics.preCheckpointStatistics(this);
            Checkpoint.setCheckpoint(this);
            statistics.postCheckpointStatistics(this);
        }

        return R_NOTDONE;
    }

    public Population preselection() {

        //ensure the elites can be saved in next generation
        ArrayList<Individual[]> elites = new ArrayList<Individual[]>(this.population.subpops.length);
        PopulationUtils.sort(this.population);

        numElites = this.parameters.getIntWithDefault(new Parameter(P_ELITES), null, 1);

        for (int pop = 0; pop < this.population.subpops.length; pop++) {
            List<Individual> tempElites = new ArrayList<>();
            for (int e = 0; e < 500; e++) {
                if (this.population.subpops[pop].individuals[e].evaluated == true)
                    tempElites.add(this.population.subpops[pop].individuals[e]);
                if (tempElites.size() == numElites)
                    break;
            }
            if (elites.size() == 0 || elites.size() == 1) {
                elites.add(pop, tempElites.toArray(new Individual[tempElites.size()]));
            } else {
                Individual[] combineElites = ArrayUtils.addAll(elites.get(pop), tempElites.toArray(new Individual[tempElites.size()]));
                elites.set(pop, combineElites);
            }
        }

        Population newPop = (Population) this.population.emptyClone();//save the population with k*populationsize individuals
        Population tempNewPop; //save the population with populationsize individuals for combining them together to newPop
        numRep = this.parameters.getIntWithDefault(new Parameter(P_REPLICATIONS), null, 1);
        for (int i = 0; i < numRep; i++) {
            tempNewPop = breeder.breedPopulation(this);
            for (int sub = 0; sub < this.population.subpops.length; sub++) {
                //combinedInds = new Individual[subpopsLength];
                //System.arraycopy(tempNewPop.subpops[sub].individuals, 0, combinedInds, 0, tempNewPop.subpops[sub].individuals.length);
                if (i == 0) {
                    newPop.subpops[sub].individuals = tempNewPop.subpops[sub].individuals;
                } else {
                    newPop.subpops[sub].individuals = ArrayUtils.addAll(newPop.subpops[sub].individuals, tempNewPop.subpops[sub].individuals);
                }
            }
        }

        population = newPop;

        //evaluate the population based on surrogate model
        evaluator.evaluatePopulation(this);
        //2021.7.29 when sorting the inds in the intermediate pop, mapping with the distance
//        if (this.generation == 0 || this.generation == 10 || this.generation == 20 || this.generation == 30
//                || this.generation == 40 || this.generation == 49) {
//            distanceFitness = new ArrayList<>();
//            for (int ind = 0; ind < this.population.subpops[0].individuals.length; ind++){
//                distanceFitness.add(new ImmutablePair<>(pcDistance.get(ind),this.population.subpops[0].individuals[ind].fitness.fitness()));
//            }
//            distanceFitness.sort(Comparator.comparingDouble(Pair::getValue)); // in a increase order, the smaller the better
//            pcDistance = new ArrayList<>();
//        }
        //================================end===========================

        PopulationUtils.sort(population); //sort the inds in the intermediate pop by estimated fitness


        //2021.7.26 this is for investigating the effectiveness of the surrogate. To see whether if there are good inds are ignored based on surrogate
        //=====================begin==================
       /* ArrayList<Pair<Integer, Double>> indexEstimatedFitness;
        indexEstimatedFitness = new ArrayList<>();
        for (int ind = 0; ind < this.population.subpops[0].individuals.length; ind++) {
            indexEstimatedFitness.add(new ImmutablePair<>(ind, this.population.subpops[0].individuals[ind].fitness.fitness()));
        }
        indexEstimatedFitness.sort(Comparator.comparingDouble(Pair::getValue)); // in a increase order, the smaller the better
        int[] indexSelectedIndsEstimated = new int[population.subpops[0].individuals.length / numRep];
        for (int i = 0; i < population.subpops[0].individuals.length / numRep; i++) {
            indexSelectedIndsEstimated[i] = indexEstimatedFitness.get(i).getKey();
        }
        //2021.7.29 only consider the top 30% inds to investigate whether the selected inds by real evaluation are selected by surrogate or not.
        int[] indexSelectedIndsEstimatedTop30Percent = new int[(int)((population.subpops[0].individuals.length / numRep) * 0.3)];
        for (int i = 0; i < indexSelectedIndsEstimatedTop30Percent.length; i++) {
            indexSelectedIndsEstimatedTop30Percent[i] = indexEstimatedFitness.get(i).getKey();
        }

        //real evaluation
        realEvaluationIntermediate = true;
        evaluator.evaluatePopulation(this);
        realEvaluationIntermediate = false;
        ArrayList<Pair<Integer, Double>> indexRealFitness;
        indexRealFitness = new ArrayList<>();
        for (int ind = 0; ind < this.population.subpops[0].individuals.length; ind++) {
            indexRealFitness.add(new ImmutablePair<>(ind, this.population.subpops[0].individuals[ind].fitness.fitness()));//this is real fitness
        }

        //2021.7.29 get the distance VS fitnessGap
        if (this.generation == 0 || this.generation == 10 || this.generation == 20 || this.generation == 30
                || this.generation == 40 || this.generation == 49) {
            ArrayList<Double> distanceFinal  = new ArrayList<>();
            ArrayList<Double> fitnessGap = new ArrayList<>();
            ArrayList<Double> realFitness = new ArrayList<>();
            ArrayList<Double> surrogateFitness = new ArrayList<>();
            ArrayList<Integer> overOrUnderEstimated = new ArrayList<>();
            for (int i = 0; i < distanceFitness.size() / numRep; i++) {
                if (distanceFitness.get(i).getKey() != Double.MAX_VALUE) { //means this inds are evaluated, 1. cleared inds 2. generate by reproductions
                    distanceFinal.add(distanceFitness.get(i).getKey());

//                    System.out.println(i + "," + distanceFitness.get(i).getKey() + "," + Math.abs(this.population.subpops[0].individuals[i].fitness.fitness() - distanceFitness.get(i).getValue()));
                    fitnessGap.add(Math.abs(this.population.subpops[0].individuals[i].fitness.fitness() - distanceFitness.get(i).getValue()));
                    realFitness.add(this.population.subpops[0].individuals[i].fitness.fitness());
                    surrogateFitness.add(distanceFitness.get(i).getValue());
                    if (this.population.subpops[0].individuals[i].fitness.fitness() - distanceFitness.get(i).getValue() > 0) {
                        overOrUnderEstimated.add(-1); //estimated a smaller value --- bad becomes good
                    } else if (this.population.subpops[0].individuals[i].fitness.fitness() - distanceFitness.get(i).getValue() < 0) {
                        overOrUnderEstimated.add(1); //estimated a larger value --- treat good as bad (serious)
                    } else {
                        overOrUnderEstimated.add(0); //exactly the same
                    }
                }
            }
            writeDistanceFitnessGap(distanceFinal, realFitness, surrogateFitness, fitnessGap, overOrUnderEstimated);
        }

        //2021.7.27 the inds in the intermediate population is sorted based on the estimated fitness. We also get the real fitness of inds---ordered
        //by the estimated fitness
        if (this.generation == 0 || this.generation == 10 || this.generation == 20 || this.generation == 30
                || this.generation == 40 || this.generation == 49) {
            writeEstimatedRealFitnessIndsIntermediatePop(indexEstimatedFitness, indexRealFitness);
        }

        indexRealFitness.sort(Comparator.comparingDouble(Pair::getValue)); // in a increase order, the smaller the better

        int[] indexSelectedIndsReal = new int[population.subpops[0].individuals.length / numRep];
        for (int i = 0; i < population.subpops[0].individuals.length / numRep; i++) {
            indexSelectedIndsReal[i] = indexRealFitness.get(i).getKey();
        }

        //2021.7.29 only consider the top 30% inds to investigate whether the selected inds by real evaluation are selected by surrogate or not.
        int[] indexSelectedIndsReal30Percent = new int[(int)((population.subpops[0].individuals.length / numRep) * 0.3)];
        for (int i = 0; i < indexSelectedIndsReal30Percent.length; i++) {
            indexSelectedIndsReal30Percent[i] = indexRealFitness.get(i).getKey();
        }

        int tempSelectedRealUnselectedEstimated = 0;
        for (int selectedReal = 0; selectedReal < indexSelectedIndsReal.length; selectedReal++) {
            if (ArrayUtils.contains(indexSelectedIndsReal, indexSelectedIndsEstimated[selectedReal]) == false) {//if the selected inds
                //by real evaluation is not selected by estimated fitness (surrogate)
                tempSelectedRealUnselectedEstimated++;
            }
        }
        selectedRealUnselectedEstimated.add(tempSelectedRealUnselectedEstimated);

        //2021.7.29
        int tempSelectedRealUnselectedEstimatedtop30Percent = 0;
        for (int selectedReal = 0; selectedReal < indexSelectedIndsReal30Percent.length; selectedReal++) {
            if (ArrayUtils.contains(indexSelectedIndsReal30Percent, indexSelectedIndsEstimatedTop30Percent[selectedReal]) == false) {//if the selected inds
                //by real evaluation is not selected by estimated fitness (surrogate)
                tempSelectedRealUnselectedEstimatedtop30Percent++;
            }
        }
        selectedRealUnselectedEstimatedTop30Percent.add(tempSelectedRealUnselectedEstimatedtop30Percent);
//        System.out.println(selectedRealUnselectedEstimated);

        if (this.generation == 0 || this.generation == 10 || this.generation == 20 || this.generation == 30
                || this.generation == 40 || this.generation == 49) {
            writeRealEstimatedFitness(indexEstimatedFitness, indexRealFitness); //ordered by index
        }*/
        //=====================================end=====================================

        for (int sub = 0; sub < this.population.subpops.length; sub++) {
            population.subpops[sub].resize(population.subpops[sub].individuals.length / numRep);
        }


        //------begin------  记录下有多少PC的fitness在pre-selection最后一个个体的左右两侧
/*        int num = 0;
        for (int a=0; a<surrogateClearingMultitreeEvaluatorVV9.duplicatedPCMinMaxFit.length; a++){
            if(surrogateClearingMultitreeEvaluatorVV9.duplicatedPCMinMaxFit[a][0] <= population.subpops[0].individuals[499].fitness.fitness()&&
                    population.subpops[0].individuals[499].fitness.fitness() <= surrogateClearingMultitreeEvaluatorVV9.duplicatedPCMinMaxFit[a][1])
                num++;
        }
        criticalBetween.add((double)num/surrogateClearingMultitreeEvaluatorVV9.duplicatedPCMinMaxFit.length);

        if(generation == numGenerations-2)
            writeCriticalBetweenRatio(criticalBetween);*/
        //------end------

        //System.out.println("最后一个个体的fitness" + population.subpops[0].individuals[499].fitness.fitness());

        for (int sub = 0; sub < this.population.subpops.length; sub++) {
            int e = 0;
            for (int replace = population.subpops[sub].individuals.length - 1; replace >= population.subpops[sub].individuals.length - numElites; replace--) {
                population.subpops[sub].individuals[replace] = elites.get(sub)[e];
                e++;
            }
        }


        return population;
    }

    //2021.6.17 save the pc distribution of individuals
    public void writePCDiffGens(ArrayList[] pcDiffGens) {
        //find the maximum list length
        int maxLength = -1;
        for (ArrayList list : pcDiffGens) {
            if (list.size() > maxLength)
                maxLength = list.size();
        }

        File pcDiffGensFile = new File(out_dir + "/job." + jobSeed + ".pcDiffGens.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(pcDiffGensFile));
//            writer.write("gen0,gen10,gen20,gen30,gen40,gen50");
            writer.newLine();
            for (int len = 0; len < maxLength; len++) {
                String toWrite = "";
                for (int index = 0; index < pcDiffGens.length - 1; index++) {
                    ArrayList<double[]> list = pcDiffGens[index];
                    if (list.size() > len) {
                        for (int ele = 0; ele < list.size(); ele++) {
                            toWrite += Arrays.toString(list.get(ele)).replaceAll("\\[", "").replaceAll(" ", "").replaceAll("]", "\n");
                        }
                    }
                    toWrite += ",";
                }
                ArrayList<double[]> list = pcDiffGens[pcDiffGens.length - 1];
                if (list.size() > len) {
                    for (int ele = 0; ele < list.size(); ele++) {
                        toWrite += Arrays.toString(list.get(ele)).replaceAll("\\[", "").replaceAll(" ", "").replaceAll("]", "\n");
                    }
                }
                writer.write(toWrite);
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writePCDiffGens(ArrayList<int[]> pcCurrentGens) {
        File pcCurrentGensFile = new File(out_dir + "/job." + jobSeed + "." + this.generation + ".pcCurrentGens.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(pcCurrentGensFile));
            writer.write("D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15,D16,D17,D18,D19,D20,D21,D22,D23,D24,D25,D26,D27,D28,D29,D30,D31," +
                    "D32,D33,D34,D35,D36,D37,D38,D39,D40");
            writer.newLine();
            String toWrite = "";
            for (int index = 0; index < pcCurrentGens.size(); index++) {
                int[] list = pcCurrentGens.get(index);
                for (int ele = 0; ele < list.length - 1; ele++) {
                    toWrite += list[ele] + ",";
                }
                toWrite += list[list.length - 1]; //delete the "," after the last column
                toWrite += "\n";
            }

            writer.write(toWrite);
            writer.newLine();

            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //2021.6.17 save the fitness distribution of individuals
    public void writeFitnessDiffGens(ArrayList[] fitnessDiffGens) {
        //find the maximum list length
        int maxLength = -1;
        for (ArrayList list : fitnessDiffGens) {
            if (list.size() > maxLength)
                maxLength = list.size();
        }

        File fitnessDiffGensFile = new File(out_dir + "/job." + jobSeed + ".fitnessDiffGens.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(fitnessDiffGensFile));
            writer.write("gen0,gen10,gen20,gen30,gen40,gen50");
            writer.newLine();
            for (int len = 0; len < maxLength; len++) {
                String toWrite = "";
                for (int index = 0; index < fitnessDiffGens.length - 1; index++) {
                    ArrayList list = fitnessDiffGens[index];
                    if (list.size() > len) {
                        toWrite += list.get(len);
                    }
                    toWrite += ",";
                }
                ArrayList list = fitnessDiffGens[fitnessDiffGens.length - 1];
                if (list.size() > len) {
                    toWrite += list.get(len);
                }
                writer.write(toWrite);
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeCorrealtionSurrogateAccuracyGens(ArrayList<Double> correlationSurrogateAccuracy) {
        File correlationSurrogateAccuracyFile = new File(out_dir + "/job." + jobSeed + ".correlationSurrogateAccuracy.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(correlationSurrogateAccuracyFile));
            writer.write("Gen,correlationSurrogateAccuracy");
            writer.newLine();
            for (int gen = 0; gen < correlationSurrogateAccuracy.size(); gen++) {
                writer.write(gen + "," + correlationSurrogateAccuracy.get(gen));
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeCorrealtionSurrogateAccuracyGensWhole(ArrayList<Double> correlationSurrogateAccuracyWhole) {
        File correlationSurrogateAccuracyWholeFile = new File(out_dir + "/job." + jobSeed + ".correlationSurrogateAccuracyWhole.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(correlationSurrogateAccuracyWholeFile));
            writer.write("Gen,correlationSurrogateAccuracyWhole");
            writer.newLine();
            for (int gen = 0; gen < correlationSurrogateAccuracyWhole.size(); gen++) {
                writer.write(gen + "," + correlationSurrogateAccuracyWhole.get(gen));
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeDistanceRankSurrogateAccuracyGens(ArrayList<Double> distanceRankSurrogateAccuracy) {
        File distanceRankSurrogateAccuracyFile = new File(out_dir + "/job." + jobSeed + ".distanceRankSurrogateAccuracy.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(distanceRankSurrogateAccuracyFile));
            writer.write("Gen,distanceRankSurrogateAccuracy");
            writer.newLine();
            for (int gen = 0; gen < distanceRankSurrogateAccuracy.size(); gen++) {
                writer.write(gen + "," + distanceRankSurrogateAccuracy.get(gen));
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeDistanceRankSurrogateAccuracyGensWhole(ArrayList<Double> distanceRankSurrogateAccuracyWhole) {
        File distanceRankSurrogateAccuracyWholeFile = new File(out_dir + "/job." + jobSeed + ".distanceRankSurrogateAccuracyWhole.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(distanceRankSurrogateAccuracyWholeFile));
            writer.write("Gen,distanceRankSurrogateAccuracyWhole");
            writer.newLine();
            for (int gen = 0; gen < distanceRankSurrogateAccuracyWhole.size(); gen++) {
                writer.write(gen + "," + distanceRankSurrogateAccuracyWhole.get(gen));
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public void writeRealEstimatedFitness(ArrayList<Pair<Integer, Double>> indexEstimatedFitness, ArrayList<Pair<Integer, Double>> indexRealFitness) {
        File indexEstimatedRealFitnessFile = new File(out_dir + "/job." + jobSeed + "." + this.generation + ".indexEstimatedRealFitness.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(indexEstimatedRealFitnessFile));
            writer.write("indID,indexEstimatedFitness,estimatedFitness,indexRealFitness,realFitness");
            writer.newLine();
            for (int gen = 0; gen < indexEstimatedFitness.size(); gen++) {
                writer.write(gen + "," + indexEstimatedFitness.get(gen).getKey() + "," + indexEstimatedFitness.get(gen).getValue()
                        + "," + indexRealFitness.get(gen).getKey() + "," + indexRealFitness.get(gen).getValue());
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //2021.7.27 write the estimated fitness and real fitness of the inds in the intermediate population
    public void writeEstimatedRealFitnessIndsIntermediatePop(ArrayList<Pair<Integer, Double>> indexEstimatedFitness, ArrayList<Pair<Integer, Double>> indexRealFitness) {
        File EstimatedRealFitnessIndsIntermediatePopFile = new File(out_dir + "/job." + jobSeed + "." + this.generation + ".estimatedRealFitnessIndsIntermediatePop.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(EstimatedRealFitnessIndsIntermediatePopFile));
            writer.write("ID,estimatedFitness,realFitness");
            writer.newLine();
            for (int gen = 0; gen < indexEstimatedFitness.size(); gen++) {
                writer.write(gen + "," + indexEstimatedFitness.get(gen).getValue() + "," + indexRealFitness.get(gen).getValue());
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeSelectedRealUnselectedEstimated(ArrayList<Integer> selectedRealUnselectedEstimated) {
        File selectedRealUnselectedEstimatedFile = new File(out_dir + "/job." + jobSeed + ".selectedRealUnselectedEstimatedWholePop.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(selectedRealUnselectedEstimatedFile));
            writer.write("gen,selectedRealUnselectedEstimated");
            writer.newLine();
            for (int gen = 0; gen < selectedRealUnselectedEstimated.size(); gen++) {
                writer.write(gen + "," + selectedRealUnselectedEstimated.get(gen));
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //2021.7.29 whether the selected top 30% inds by real evaluation are selected by surrogate
    public void writeSelectedRealUnselectedEstimatedtop30Percent(ArrayList<Integer> selectedRealUnselectedEstimatedTop30Percent) {
        File selectedRealUnselectedEstimatedTop30PercentFile = new File(out_dir + "/job." + jobSeed + ".selectedRealUnselectedEstimatedtop30Percent.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(selectedRealUnselectedEstimatedTop30PercentFile));
            writer.write("gen,selectedRealUnselectedEstimatedTop30Percent");
            writer.newLine();
            for (int gen = 0; gen < selectedRealUnselectedEstimatedTop30Percent.size(); gen++) {
                writer.write(gen + "," + selectedRealUnselectedEstimatedTop30Percent.get(gen));
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeDistanceFitnessGap(ArrayList<Double> distanceFinal, ArrayList<Double> realFitness, ArrayList<Double> surrogateFitness, ArrayList<Double> fitnessGap, ArrayList<Integer> overOrUnderEstimated) {
        File distanceFitnessGapFile = new File(out_dir + "/job." + jobSeed + "." + this.generation + ".distanceFitnessGap.csv"); // jobSeed = 0
//        File distanceFitnessGapFile = new File("job." + jobSeed + "." + this.generation + ".distanceFitnessGap.xlsx"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(distanceFitnessGapFile));
            writer.write("Gen,distance,realFitness,surrogateFitness,fitnessGap,overOrUnderEstimated");
            writer.newLine();
            for (int gen = 0; gen < distanceFinal.size(); gen++) {
                writer.write(gen + "," + distanceFinal.get(gen) + "," + realFitness.get(gen) + "," + surrogateFitness.get(gen) + "," + fitnessGap.get(gen) + "," + overOrUnderEstimated.get(gen));
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //2021.8.3 write the average fitness value at each generation
    public static void writeAverageFitnessGensToFile(ArrayList<Double> averageFitnessGens) {
        File averageFitnessGensFile = new File(out_dir + "/job." + jobSeed + ".averageFitnessGens.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(averageFitnessGensFile));
            writer.write("Gen,averageFitness");
            writer.newLine();
            for (int i = 0; i < averageFitnessGens.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + "," + averageFitnessGens.get(i) + "\n");
            }
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void writeListSize(ArrayList<Integer> listPCSize, ArrayList<Integer> mostPCDuplicatesNum) {
//        Parameter p;
        // Get the job seed.
//        p = new Parameter("seed").push("" + 0);
//        long jobSeed = parameters.getLongWithDefault(p, null, 0);

        File listSize = new File(out_dir + "/job." + GPRuleEvolutionState.jobSeed + ".listSizeTask.csv");

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(listSize));
            writer.write("generation,listSizeTask0,mostPCDuplicatesNum");
            writer.newLine();
            for (int i = 0; i < listPCSize.size(); i++) {
                writer.write(i + "," + listPCSize.get(i) + "," + mostPCDuplicatesNum.get(i));
                writer.newLine();
            }
            listPCSize.clear();
            mostPCDuplicatesNum.clear();
//            writer.write(state.generation-1 + "," + 0 + "," + 0);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    //2021.4.16 fzhang save the diversity value to csv
    public static void writeDiversityToFile(ArrayList<Double> entropyDiversity) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + GPRuleEvolutionState.jobSeed + ".diversity.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Gen,diversitySubpop0");
            writer.newLine();
            for (int i = 0; i < entropyDiversity.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + entropyDiversity.get(i) + "\n");
            }
            entropyDiversity.clear();
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void writeCriticalBetweenRatio(ArrayList<Double> criticalBetween) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + GPRuleEvolutionState.jobSeed + ".criticalBetweenRatio.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("generation,criticalBetweenRatio");
            writer.newLine();
            for (int i = 0; i < criticalBetween.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + criticalBetween.get(i) + "\n");
            }
            criticalBetween.clear();
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void writeSumTrainTime(ArrayList<Double> sumTrainTime) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + GPRuleEvolutionState.jobSeed + ".sumTrainTime.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("generation,sumTrainTime");
            writer.newLine();
            for (int i = 0; i < sumTrainTime.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + sumTrainTime.get(i) + "\n");
            }
            sumTrainTime.clear();
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
