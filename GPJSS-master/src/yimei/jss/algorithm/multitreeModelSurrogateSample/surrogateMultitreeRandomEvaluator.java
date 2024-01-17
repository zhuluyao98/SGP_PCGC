package yimei.jss.algorithm.multitreeModelSurrogateSample;

import ec.EvolutionState;
import ec.Individual;
import ec.Subpopulation;
import ec.simple.SimpleEvaluator;
import ec.simple.SimpleProblemForm;
import ec.util.Parameter;
import org.apache.commons.lang3.ArrayUtils;
import yimei.jss.gp.GPRuleEvolutionState;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.niching.*;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

import static yimei.jss.algorithm.multitreeModelSurrogate.PhenotypicGPRuleEvolutionState.realEvaluationIntermediate;
import static yimei.jss.algorithm.multitreeModelSurrogateSample.PhenotypicGPRuleEvolutionState.indexEvaluated;
import static yimei.jss.algorithm.multitreeModelSurrogateSample.PhenotypicGPRuleEvolutionState.writeDiversityToFile;


/**
 * Created by luyao on 2022.8.27.
 * 只随机选取一半个体进行评估，然后建立代理模型评估剩余的个体，同时中间种群先设置为 1
 */
//a class can not extend from more than one class
public class surrogateMultitreeRandomEvaluator extends SimpleEvaluator {
    public static final String P_RADIUS = "radius";
    public static final String P_CAPACITY = "capacity";

//    public static final String P_PERCENTAGESURRACC = "percentage-surrogateAcc";
//    protected double percentageSurrogateAcc;

    //fzhang 2018.10.9 to get the pre-generation valu

    public static ArrayList<Integer> listPCSize = new ArrayList<>();
    public static ArrayList<Integer> mostPCDuplicatesNum = new ArrayList<>();

    public static final String P_PRE_GENERATIONS = "pre-generations";
    //ArrayList<Double[]> fitnessesForModel = new ArrayList<>();
    double[] tempfitnessesForModel; //for transfer the fitness into the surrogate model to estimate the fitness
    int[][] tempindsCharListsMultiTree;

    protected boolean clear = true;
    protected boolean nonIntermediatePop = true;

    protected double radius;
    protected int capacity;

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

                if(state.generation > 0)
                    PopulationUtils.sort(state.population);

                int[][] indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic
                double diversityValue = PopulationUtils.entropy(indsCharListsMultiTree);
                entropyDiversity.add(diversityValue);

                //---------test---------------
                //luyao 8/7/2022 用hashmap记录下每个重复PC的个数
                HashMap<int[], Integer> duplicates = new HashMap<>();
                for (int i = 0; i < (indsCharListsMultiTree).length; i++) {
                    int judge = 0;
                    for (int[] k : duplicates.keySet()) {
                        if (Arrays.equals(k, indsCharListsMultiTree[i])) {
                            judge = 1;
                            int tempValue = duplicates.get(k);
                            tempValue++;
                            duplicates.replace(k, tempValue);
                            break;
                        }
                    }
                    if (judge == 0)
                        duplicates.put(indsCharListsMultiTree[i], 1);
                }
                ArrayList<Map.Entry<int[], Integer>> listPcClusterFit = new ArrayList<>(duplicates.entrySet());
                //通过list对hashmap排序
                Collections.sort(listPcClusterFit, new Comparator<ConcurrentHashMap.Entry<int[], Integer>>() {
                    @Override
                    public int compare(ConcurrentHashMap.Entry<int[], Integer> o1, ConcurrentHashMap.Entry<int[], Integer> o2) {
                        //按照value值，从大到小排序
                        return o2.getValue() - o1.getValue();
                        //按照value值，用compareTo()方法默认是从小到大排序
                        //return o1.getValue().compareTo(o2.getValue());
                    }
                });

                listPCSize.add(listPcClusterFit.size());
                mostPCDuplicatesNum.add(listPcClusterFit.get(0).getValue());

//                int evaluatedNum = (int) (state.population.subpops[0].individuals.length * 0.5);
                int evaluatedNum = listPcClusterFit.size();
                for (int a = 0; a <state.population.subpops[0].individuals.length ; a++) {
                    if (state.generation == 0) {
                        indexEvaluated.add(a);
                    } else {
                        int random = state.random[0].nextInt(state.population.subpops[0].individuals.length);
                        if (indexEvaluated.contains(random))
                            continue;
                        else
                            indexEvaluated.add(random);
//                        indexEvaluated.add(a);
                        if (indexEvaluated.size() == evaluatedNum)
                            break;
                    }
                }

                super.evaluatePopulation(state);

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


                //用surrogate评估剩余的个体
                for (int sub = 0; sub < state.population.subpops.length; sub++) {
                    for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //
                    {
                        Individual individual = state.population.subpops[sub].individuals[i];
                        if (indsCharListsMultiTree.length != 0) {
                            //KNN
                            //===============================start==============================================
                            double dMin = Double.MAX_VALUE;
                            int index = 0;
                            if (!(indexEvaluated.contains(i))) {
                                int[] pcIntermediate = indsCharListsMultiTree[i];
                                //calculate the fitness based on surrogate model
                                for (int pc = 0; pc < indsCharListsMultiTreeEvaluated.length; pc++) {
                                    int[] pcModel = indsCharListsMultiTreeEvaluated[pc];
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




//                if (state.generation == state.numGenerations - 1) {
                if (((GPRuleEvolutionState)state).totalTime >= ((GPRuleEvolutionState)state).endTime) {
                    PhenotypicGPRuleEvolutionState.writeListSize(listPCSize , mostPCDuplicatesNum);
                    writeDiversityToFile(entropyDiversity);
                }
                //-------将list的PC记录在temp中
//                ArrayList<int[]> temp = new ArrayList<>();
//                for (int i=0; i<listPcClusterFit.size() ; i++){
//                    temp.add(listPcClusterFit.get(i).getKey());
//                }
//                tempList.add(temp);
//                System.out.println( );
//                if(tempList.size() == 50){
//                    int index = 48;
//                int num = 0;
//                for (int a=0 ; a<tempList.get(index).size() ; a++){
//                    for (int b = 0; b<tempList.get(index+1).size(); b++ ){
//                        if(Arrays.equals(tempList.get(index).get(a) , tempList.get(index+1).get(b)))
//                            num++;
//                    }
//                }
//                    System.out.println(tempList.get(index).size() + tempList.get(index+1).size() - num);
//                }
 /*               for (int i = 0; i < listPcClusterFit.size(); i++) {

                    if (listPcClusterFit.get(i).getValue() > 1) {
                        ArrayList<Integer> samePCIndex = new ArrayList<>();
                        for (int j = 0; j < indsCharListsMultiTree.length; j++) {
                            if (Arrays.equals(listPcClusterFit.get(i).getKey(), indsCharListsMultiTree[j]))
                                samePCIndex.add(j);
                        }
                        double[] tempDateSet = new double[samePCIndex.size()];
                        for (int a = 0; a < samePCIndex.size(); a++) {
                            tempDateSet[a] = state.population.subpops[0].individuals[samePCIndex.get(a)].fitness.fitness();
                            System.out.print(tempDateSet[a] + " , ");
                        }
                        System.out.println(" ");
//                        System.out.println(Arrays.stream(tempDateSet).min().getAsDouble() + " , " +  Arrays.stream(tempDateSet).max().getAsDouble());
                        samePCIndex.clear();
                    }
                }*/

                //---------end---------------


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

        //-----begin----- 输出、看中间种群的PC值对应的fitness数
/*        HashMap<int[], Integer> duplicates = new HashMap<>();
        for (int i = 0; i < (indsCharListsIntermediatePop).length; i++) {
            int judge = 0;
            for (int[] k : duplicates.keySet()) {
                if (Arrays.equals(k, indsCharListsIntermediatePop[i])) {
                    judge = 1;
                    int tempValue = duplicates.get(k);
                    tempValue++;
                    duplicates.replace(k, tempValue);
                    break;
                }
            }
            if (judge == 0)
                duplicates.put(indsCharListsIntermediatePop[i], 1);
        }
        ArrayList<Map.Entry<int[], Integer>> listPcClusterFitMidPop = new ArrayList<>(duplicates.entrySet());
        //通过list对hashmap排序
        Collections.sort(listPcClusterFitMidPop, new Comparator<HashMap.Entry<int[], Integer>>() {
            @Override
            public int compare(HashMap.Entry<int[], Integer> o1, HashMap.Entry<int[], Integer> o2) {
                //按照value值，从大到小排序
                return o2.getValue() - o1.getValue();
                //按照value值，用compareTo()方法默认是从小到大排序
                //return o1.getValue().compareTo(o2.getValue());
            }
        });*/
//        int Num0 = 0; int Num1 = 0;
//        for (int a = 0; a < listPcClusterFitMidPop.size() ; a++){
//            if(listPcClusterFitMidPop.get(a).getValue() > 2)
//                Num0 += listPcClusterFitMidPop.get(a).getValue();
//            if(listPcClusterFitMidPop.get(a).getValue() > 1)
//                Num1 += listPcClusterFitMidPop.get(a).getValue();
//        }
//        System.out.println(Num0 + Num1);
        //-----end----- 输出、看中间种群的PC值对应的fitness数

        int num = 0;
        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //
            {
                Individual individual = state.population.subpops[sub].individuals[i];
                if (indsCharListsMultiTree.length != 0) {
                    //KNN
                    //===============================start==============================================
                    double dMin = Double.MAX_VALUE;
                    int index = 0;
                    if (individual.evaluated == false) {
                        num++;
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
                    if (state.generation == 0 || state.generation == 10 || state.generation == 20 || state.generation == 30
                            || state.generation == 40 || state.generation == 49) {
                        pcDistance.add(dMin);
                    }
                } else {
                    ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE);
                }
            }
            //System.out.println(num);
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
            for (int a = fp; a < indexEvaluated.size(); a++) {
                int x = indexEvaluated.get(a);
                if (realEvaluationIntermediate == true) { //2021.7.27 limit this to intermediate population with clearing
                    if (inds[x].evaluated == false) {
                        p.evaluate(state, inds[x], pop, threadnum);
                    } else {
                        continue;
                    }
                } else {
                    p.evaluate(state, inds[x], pop, threadnum);
                }
            }
        }

        ((ec.Problem) p).finishEvaluating(state, threadnum);
    }


    public void setClear(boolean clear) {
        this.clear = clear;
    }
}

