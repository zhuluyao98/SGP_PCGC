package yimei.jss.algorithm.multitreeModelSurrogateSample;

import ec.EvolutionState;
import ec.Individual;
import ec.Subpopulation;
import ec.simple.SimpleEvaluator;
import ec.simple.SimpleProblemForm;
import ec.util.Parameter;
import org.apache.commons.lang3.ArrayUtils;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.niching.*;

import java.util.ArrayList;
import java.util.Map;

import static yimei.jss.algorithm.multitreeModelSurrogate.PhenotypicGPRuleEvolutionState.realEvaluationIntermediate;
import static yimei.jss.algorithm.multitreeModelSurrogate.PhenotypicGPRuleEvolutionState.writeDiversityToFile;

//用上一代的surrogate评估下一代中最佳的250个个体用于真实评估
public class pastSurrogateMultitreeEvaluator extends SimpleEvaluator {


    public ArrayList<Map.Entry<int[], Integer>> listPcClusterFit;

    public static ArrayList<Integer> indexEvaluated = new ArrayList<>();
    public ArrayList<int[][]> genPCModel = new ArrayList<>();
    public ArrayList<double[]> genFitModel = new ArrayList<>();
    public static ArrayList<Integer> listPCSize = new ArrayList<>();
    public static ArrayList<Integer> mostPCDuplicatesNum = new ArrayList<>();
    public static final String P_RADIUS = "radius";
    public static final String P_CAPACITY = "capacity";

//    public static final String P_PERCENTAGESURRACC = "percentage-surrogateAcc";
//    protected double percentageSurrogateAcc;

    //fzhang 2018.10.9 to get the pre-generation value
    public static final String P_PRE_GENERATIONS = "pre-generations";
    //ArrayList<Double[]> fitnessesForModel = new ArrayList<>();
    double[] tempfitnessesForModel; //for transfer the fitness into the surrogate model to estimate the fitness
    int[][] tempindsCharListsMultiTree;

    protected boolean clear = true;
    protected boolean nonIntermediatePop = true;

    protected double radius;
    protected int capacity;

    public static ArrayList<Double> correlationSurrogateAccuracy = new ArrayList<>();
    public static ArrayList<Double> correlationSurrogateAccuracyWhole = new ArrayList<>();
    public static ArrayList<Double> distanceRankSurrogateAccuracy = new ArrayList<>();
    public static ArrayList<Double> distanceRankSurrogateAccuracyWhole = new ArrayList<>();

    //    protected PhenoCharacterisation[] phenoCharacterisation;
    public static PhenoCharacterisation[] phenoCharacterisation;

    public static ArrayList<Double> entropyDiversity = new ArrayList<>();

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

    public static ArrayList<Double> pcDistance = new ArrayList<>();

    public static ArrayList<Double> averageFitnessGens = new ArrayList<>();

    public void setup(final EvolutionState state, final Parameter base) {
        super.setup(state, base);

        radius = state.parameters.getDoubleWithDefault(
                base.push(P_RADIUS), null, 0.0);
        capacity = state.parameters.getIntWithDefault(
                base.push(P_CAPACITY), null, 1);
//        percentageSurrogateAcc = state.parameters.getDoubleWithDefault(
//                base.push(P_PERCENTAGESURRACC), null, 0.0);

        String filePath = state.parameters.getString(new Parameter("filePath"), null);
        //It's a little tricky to know whether we have 1 or 2 populations here, so we will assume
        //2 for the purpose of the phenoCharacterisation, and ignore the second object if only
        //1 is used
        phenoCharacterisation = new PhenoCharacterisation[2];

        //this is the baseline PhenoCharacterisation with baseline rule "SPT" "WIQ", it will be set again by the best rule, so it is useful here
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

        //if it is a normal evaluation
        if (realEvaluationIntermediate == true) {
//            multitreeClearingFirst.clearPopulation(state, radius, capacity, phenoCharacterisation);
            super.evaluatePopulation(state);
        } else {
            if (nonIntermediatePop) {

                int[][] indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic
                double diversityValue = PopulationUtils.entropy(indsCharListsMultiTree);
                entropyDiversity.add(diversityValue);


                //先用上代的surrogate评估出最好的250个个体，然后对齐进行real-evaluation
                if (state.generation == 0) {

                    for (int a = 0; a < (state.population.subpops[0].individuals.length); a++)
                        indexEvaluated.add(a);

                    super.evaluatePopulation(state);
                    //get the fitness for training model
                    double[] fitnessesForModel = new double[state.population.subpops[0].individuals.length];
                    for (int ind = 0; ind < state.population.subpops[0].individuals.length; ind++) {
                        fitnessesForModel[ind] = state.population.subpops[0].individuals[ind].fitness.fitness();
                    }

                    int[][] indsCharListsMultiTreeCopy = indsCharListsMultiTree.clone();

                    // int removeIdx = 0;
                    for (int i = 0; i < fitnessesForModel.length; i++) {
                        if (fitnessesForModel[i] == Double.MAX_VALUE) {
                   /* indsCharListsMultiTree = ArrayUtils.remove(indsCharListsMultiTree, i - removeIdx);
                    fitnessesForModel = ArrayUtils.remove(fitnessesForModel, i - removeIdx);*/
                            indsCharListsMultiTreeCopy = ArrayUtils.remove(indsCharListsMultiTreeCopy, i);
                            fitnessesForModel = ArrayUtils.remove(fitnessesForModel, i);
                            i--;
                            // removeIdx++;
                        }
                    }

                    tempfitnessesForModel = fitnessesForModel;
                    tempindsCharListsMultiTree = indsCharListsMultiTreeCopy;
                    nonIntermediatePop = false;
                } else {
                    for (int sub = 0; sub < state.population.subpops.length; sub++) {
                        for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //
                        {
                            Individual individual = state.population.subpops[sub].individuals[i];
                            //KNN
                            //===============================start==============================================
                            double dMin = Double.MAX_VALUE;
                            int index = 0;
                            if (individual.evaluated == false) {
                                int[] pcIntermediate = indsCharListsMultiTree[i];
                                //calculate the fitness based on surrogate model
                                for (int pc = 0; pc < genPCModel.get(state.generation - 1).length; pc++) {
                                    int[] pcModel = genPCModel.get(state.generation - 1)[pc];
                                    double d = PhenoCharacterisation.distance(pcIntermediate, pcModel);
                                    if (d == 0) {
                                        dMin = d;
                                        index = pc;
                                        break;
                                    }

                                    if (d < dMin) {
                                        dMin = d;
                                        index = pc;
                                    }
                                }
                                ((Clearable) individual.fitness).surrogateFitness(genFitModel.get(state.generation - 1)[index]);
                                //individual.evaluated = true;
                                //individual.fitness.trials = new ArrayList();
                                //individual.fitness.trials.add(individual.fitness.fitness());
                            }
                        }
                    }


                    PopulationUtils.sort(state.population);

                    super.evaluatePopulation(state);

                    for (int a = 0; a < (state.population.subpops[0].individuals.length * 0.5); a++)
                        indexEvaluated.add(a);

                    //get the fitness for training model
                    double[] fitnessesForModel = new double[indexEvaluated.size()];
                    int[][] indsCharListsMultiTreeEvaluated = new int[indexEvaluated.size()][40];
                    for (int ind = 0; ind < indexEvaluated.size(); ind++) {
                        fitnessesForModel[ind] = state.population.subpops[0].individuals[indexEvaluated.get(ind)].fitness.fitness();
                        indsCharListsMultiTreeEvaluated[ind] = indsCharListsMultiTree[indexEvaluated.get(ind)];
                    }

                    // int removeIdx = 0;
                    for (int i = 0; i < fitnessesForModel.length; i++) {
                        if (fitnessesForModel[i] == Double.MAX_VALUE) {
                   /* indsCharListsMultiTree = ArrayUtils.remove(indsCharListsMultiTree, i - removeIdx);
                    fitnessesForModel = ArrayUtils.remove(fitnessesForModel, i - removeIdx);*/
                            indsCharListsMultiTreeEvaluated = ArrayUtils.remove(indsCharListsMultiTreeEvaluated, i);
                            fitnessesForModel = ArrayUtils.remove(fitnessesForModel, i);
                            i--;
                            // removeIdx++;
                        }
                    }
                    tempfitnessesForModel = fitnessesForModel;
                    tempindsCharListsMultiTree = indsCharListsMultiTreeEvaluated;
                    nonIntermediatePop = false;
                }

                //将每一代的surrogate保存
                genPCModel.add(tempindsCharListsMultiTree);
                genFitModel.add(tempfitnessesForModel);

                //-----begin-----
                //luyao 8/7/2022 用hashmap记录下每个重复PC的个数
                /*ConcurrentHashMap<int[], Integer> duplicates = new ConcurrentHashMap<>();
                for(int i = 0; i < (indsCharListsMultiTree).length; i++){
                    int judge = 0;
                    for(int[] k : duplicates.keySet()){
                        if(Arrays.equals(k,indsCharListsMultiTree[i])){
                            judge = 1;
                            int tempValue = duplicates.get(k);
                            tempValue++;
                            duplicates.replace(k,tempValue);
                            break;
                        }
                    }
                    if(judge == 0)
                        duplicates.put(indsCharListsMultiTree[i],1 );
                }
                listPcClusterFit = new ArrayList<>(duplicates.entrySet());
                //通过list对hashmap排序
                Collections.sort(listPcClusterFit, new Comparator<ConcurrentHashMap.Entry<int[], Integer>>()
                {
                    @Override
                    public int compare(ConcurrentHashMap.Entry<int[], Integer> o1, ConcurrentHashMap.Entry<int[], Integer> o2)
                    {
                        //按照value值，从大到小排序
                        return o2.getValue() - o1.getValue();
                        //按照value值，用compareTo()方法默认是从小到大排序
                        //return o1.getValue().compareTo(o2.getValue());
                    }
                });

                listPCSize.add(listPcClusterFit.size());
                mostPCDuplicatesNum.add(listPcClusterFit.get(0).getValue());*/
                //------end------
                if (state.generation == state.numGenerations - 1) {
                    //writeListSize(this.listPCSize, this.mostPCDuplicatesNum);
                    writeDiversityToFile(entropyDiversity);
                }


            } else {
                //multitreeClearingFirst.clearPopulation(state, radius, capacity, phenoCharacterisation);//the individuals do not have fitness or they are evaluated yet
                //but has a fitness there
                this.evaluatePopulation(state, tempindsCharListsMultiTree, tempfitnessesForModel);
                nonIntermediatePop = true;
            }
        }

    }

    //==========================================KNN========================================
    public void evaluatePopulation(final EvolutionState state, int[][] indsCharListsMultiTree, double[] fitnessesForModel) {

        int[][] indsCharListsIntermediatePop = phenotypicForSurrogate.muchBetterPhenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //
            {
                if (!indexEvaluated.contains(i)) {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if (indsCharListsMultiTree.length != 0) {
                        //KNN
                        //===============================start==============================================
                        double dMin = Double.MAX_VALUE;
                        int index = 0;
                        if (individual.evaluated == false) {
                            int[] pcIntermediate = indsCharListsIntermediatePop[i];
                            //calculate the fitness based on surrogate model
                            for (int pc = 0; pc < indsCharListsMultiTree.length; pc++) {
                                int[] pcModel = indsCharListsMultiTree[pc];
                                double d = PhenoCharacterisation.distance(pcIntermediate, pcModel);
                                if (d == 0) {
                                    dMin = d;
                                    index = pc;
                                    break;
                                }

                                if (d < dMin) {
                                    dMin = d;
                                    index = pc;
                                }
                            }
                            ((Clearable) individual.fitness).surrogateFitness(fitnessesForModel[index]);
                            //individual.evaluated = true;
                            //individual.fitness.trials = new ArrayList();
                            //individual.fitness.trials.add(individual.fitness.fitness());
                        }
                    } else {
                        ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE);

                    }
                }
            }
        }

    }

    protected void evalPopChunk(EvolutionState state, int[] numinds, int[] from,
                                int threadnum, SimpleProblemForm p) {
        ((ec.Problem) p).prepareToEvaluate(state, threadnum);

        Subpopulation[] subpops = state.population.subpops;
        int len = subpops.length;

        for (int pop = 0; pop < len; pop++) {
            // start evaluatin'!
            int fp = from[pop];
            int upperbound = fp + numinds[pop];
            Individual[] inds = subpops[pop].individuals;
                for (int x = fp; x < indexEvaluated.size(); x++)
                    p.evaluate(state, inds[indexEvaluated.get(x)], pop, threadnum);

        }

        ((ec.Problem) p).finishEvaluating(state, threadnum);
    }

}
