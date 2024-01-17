package yimei.jss.algorithm.multitreeModelSurrogateSample;

import ec.EvolutionState;
import ec.Individual;
import ec.Subpopulation;
import ec.gp.GPIndividual;
import ec.gp.GPNode;
import ec.simple.SimpleEvaluator;
import ec.simple.SimpleProblemForm;
import ec.util.Parameter;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import yimei.jss.gp.GPRuleEvolutionState;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.niching.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import static yimei.jss.algorithm.multitreeModelSurrogate.PhenotypicGPRuleEvolutionState.realEvaluationIntermediate;
import static yimei.jss.algorithm.multitreeModelSurrogateSample.PhenotypicGPRuleEvolutionState.indexEvaluated;
import static yimei.jss.gp.GPRun.out_dir;


/**
 * Created by luyao on 2022.8.27.
 * 对同一PC进行niching, 然后每个niche中只对最小rule size的个体
 */
//a class can not extend from more than one class
public class surrogateMultitreePCNichingEvaluator extends SimpleEvaluator {
    public static final String P_RADIUS = "radius";
    public static final String P_CAPACITY = "capacity";

//    public static final String P_PERCENTAGESURRACC = "percentage-surrogateAcc";
//    protected double percentageSurrogateAcc;

    //fzhang 2018.10.9 to get the pre-generation value
    public static final String P_PRE_GENERATIONS = "pre-generations";
    //ArrayList<Double[]> fitnessesForModel = new ArrayList<>();
    double[][] tempfitnessesForModel; //for transfer the fitness into the surrogate model to estimate the fitness
    int[][] tempindsCharListsMultiTree;

    int[][] listPCIndex;

    double correlationCoeffcient;

    double[][][] listPCTermianl;

    public ArrayList<Integer> listPCSize = new ArrayList<>();
    public ArrayList<Integer> duplicatePCNum = new ArrayList<>();

    public ArrayList<Integer> topPCNum = new ArrayList<>();
    public ArrayList<Integer> topDuplicatePCNum = new ArrayList<>();
    public ArrayList<Integer> topDuplicatePCNichingNum = new ArrayList<>();       //top30%重复PC中niche>2的数量
    public ArrayList<Integer> topExtraEvaluatedNum = new ArrayList<>();           //额外评估的个体数
    //public ArrayList<Integer> topFitDistance = new ArrayList<>();              //niche中的fitness之差

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

                if (state.generation > 0)
                    PopulationUtils.sort(state.population);

                int[][] indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic
                double diversityValue = PopulationUtils.entropy(indsCharListsMultiTree);
                entropyDiversity.add(diversityValue);

                //将新种群的所有个体都inNichied = false,避免对后续niching造成干扰！
                for (int a = 0; a < state.population.subpops[0].individuals.length; a++) {
                    Individual inds = state.population.subpops[0].individuals[a];
                    ((ClearingMultiObjectiveFitness) inds.fitness).notInNiche();
                }

                correlationCoeffcient = state.parameters.getDoubleWithDefault(new Parameter("correlation-coeffcient"),null,0.5);
//                correlationCoeffcient = 0.8 - 0.3*Math.pow(1-(((GPRuleEvolutionState)state).totalTime/9100),2);
//                if(state.generation <= 30)
//                    correlationCoeffcient = 0.5;
//                else if (state.generation <= 150)
//                    correlationCoeffcient = 0.7;
//                else
//                    correlationCoeffcient = 0.85;

                //---------test---------------
                //luyao 8/7/2022 用hashmap记录下每个重复PC的个数
                LinkedHashMap<int[], Integer> duplicates = new LinkedHashMap<>();
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
                ArrayList<Map.Entry<int[], Integer>> listPCNum = new ArrayList<>(duplicates.entrySet());
                //通过list对hashmap排序
/*                Collections.sort(listPcClusterFit, new Comparator<HashMap.Entry<int[], Integer>>() {
                    @Override
                    public int compare(HashMap.Entry<int[], Integer> o1, HashMap.Entry<int[], Integer> o2) {
                        //按照value值，从大到小排序
                        return o2.getValue() - o1.getValue();
                        //按照value值，用compareTo()方法默认是从小到大排序
                        //return o1.getValue().compareTo(o2.getValue());
                    }
                });*/

                //记录PC的种类个数，独一无二PC的个数
                int duplicatePcNum = 0;
                for (int i = 0; i < listPCNum.size(); i++) {
                    if (listPCNum.get(i).getValue() > 1) {
                        duplicatePcNum++;
                    }
                }
                listPCSize.add(listPCNum.size());
                duplicatePCNum.add(duplicatePcNum);

                //将独一无二PC的个体挑选出来，建立代理模型
                if (state.generation == 0) {
                    listPCIndex = new int[500][1];
                    for (int a = 0; a < 500; a++) {
                        indexEvaluated.add(a);
                        listPCIndex[a][0] = a;
                    }
                } else {
                    //只有在第1代往后才记录Top30%的个体
                    int pcLastTopIndex = 0;
                    //记录下Top30%的PC最后一个在list中的索引
                    int topPcNum = (int) (state.population.subpops[0].individuals.length * 0.3);
                    int topLastIndex = topPcNum - 1;
                    for (int a = 0; a < listPCNum.size(); a++) {
                        if (Arrays.equals(indsCharListsMultiTree[topLastIndex], listPCNum.get(a).getKey())) {
                            pcLastTopIndex = a;
                            break;
                        }
                    }

                    topPCNum.add(pcLastTopIndex + 1);
                    //记录下Top30%的PC有多少做了niching
                    int topDuplicatePcNum = 0;
                    int nicheNumOver1 = 0;
                    int topExtraEvaluatedIndNum = 0;
                    listPCIndex = new int[listPCNum.size()][];
                    listPCTermianl = new double[listPCNum.size()][][];
                    for (int i = 0; i < listPCNum.size(); i++) {
                        if (listPCNum.get(i).getValue() > 1) {

                            ArrayList<Integer> samePCIndex = new ArrayList<>();
                            for (int j = 0; j < indsCharListsMultiTree.length; j++) {
                                if (Arrays.equals(listPCNum.get(i).getKey(), indsCharListsMultiTree[j]))
                                    samePCIndex.add(j);
                                if (samePCIndex.size() == listPCNum.get(i).getValue())
                                    break;
                            }

                            if (i < 500) {
//                                if (i < pcLastTopIndex) {        //说明是Top30%的PC,需要niching
                                topDuplicatePcNum++;
                                ArrayList<ArrayList<Individual>> pcNiches = new ArrayList<>();
                                ArrayList<Integer> pcNicheIndsIndex = new ArrayList<>();     //记录该PC所分各个niching中真实评估的个体索引
                                int[] pcRuleSize = new int[samePCIndex.size()];
                                int[] pcIndex = new int[samePCIndex.size()];
                                for (int a = 0; a < samePCIndex.size(); a++) {
                                    Individual inds = state.population.subpops[0].individuals[samePCIndex.get(a)];
                                    pcRuleSize[a] = ((GPIndividual) inds).trees[0].child.numNodes(GPNode.NODESEARCH_ALL) + ((GPIndividual) inds).trees[1].child.numNodes(GPNode.NODESEARCH_ALL);
                                    pcIndex[a] = samePCIndex.get(a);
                                }
                                int[][] sortedPCIndexRule = new int[pcRuleSize.length][2];
                                for (int a = 0; a < pcRuleSize.length; a++) {
                                    sortedPCIndexRule[a][0] = pcRuleSize[a];
                                    sortedPCIndexRule[a][1] = pcIndex[a];
                                }
                                //在Arrays.sort()重载方法中，new Comparator对象
                                Arrays.sort(sortedPCIndexRule, new Comparator<int[]>() {
                                    @Override
                                    public int compare(int[] o1, int[] o2) {
                                        return o1[0] - o2[0];
                                    }
                                });
                                for (int a = 0; a < pcRuleSize.length; a++) {
                                    pcRuleSize[a] = sortedPCIndexRule[a][0];
                                    pcIndex[a] = sortedPCIndexRule[a][1];
                                }
                                ArrayList<double[]> pcTerminalFrequency = new ArrayList<>();
                                for (int a = 0; a < pcIndex.length; a++) {
                                    ArrayList<Individual> pcOneNiches = new ArrayList<>();
                                    Individual inds = state.population.subpops[0].individuals[pcIndex[a]];
                                    double[] indsTerminalFrequencySet = getNumTerminalSet(inds);
                                    if (((ClearingMultiObjectiveFitness) inds.fitness).isInNiche()) {
                                        continue;
                                    } else {
                                        pcOneNiches.add(inds);
                                        indexEvaluated.add(pcIndex[a]);
                                        pcNicheIndsIndex.add(pcIndex[a]);
                                        ((ClearingMultiObjectiveFitness) inds.fitness).inNiche();
                                    }
                                    for (int b = 0; b < pcIndex.length; b++) {
                                        Individual inds2 = state.population.subpops[0].individuals[pcIndex[b]];
                                        if (((ClearingMultiObjectiveFitness) inds2.fitness).isInNiche()) {
                                            continue;
                                        }
                                        double[] indsTerminalFrequencySet2 = getNumTerminalSet(inds2);
                                        double correlation = new SpearmansCorrelation().correlation(indsTerminalFrequencySet, indsTerminalFrequencySet2);
                                        if (correlation >= correlationCoeffcient) {
                                            pcOneNiches.add(inds2);
                                            ((ClearingMultiObjectiveFitness) inds2.fitness).inNiche();
                                        }

                                    }
                                    pcNiches.add(pcOneNiches);
                                }
                                listPCIndex[i] = new int[pcNicheIndsIndex.size()];
                                listPCTermianl[i] = new double[pcNicheIndsIndex.size()][20];
                                for (int a = 0; a < pcNicheIndsIndex.size(); a++) {
                                    listPCIndex[i][a] = pcNicheIndsIndex.get(a);
                                    listPCTermianl[i][a] = getNumTerminalSet(state.population.subpops[0].individuals[listPCIndex[i][a]]);
                                }

                                if (pcNicheIndsIndex.size() > 1) {
                                    nicheNumOver1++;
                                    topExtraEvaluatedIndNum += (pcNicheIndsIndex.size() - 1);
                                }


                            } else {                         //后70%的PC,只需要取出rule size最小的Index
                                int minRuleSize = 1000;
                                int selectedIndex = 0;
                                for (int a : samePCIndex) {
                                    Individual inds = state.population.subpops[0].individuals[a];
                                    int ruleSize = ((GPIndividual) inds).trees[0].child.numNodes(GPNode.NODESEARCH_ALL) + ((GPIndividual) inds).trees[1].child.numNodes(GPNode.NODESEARCH_ALL);
                                    if (ruleSize < minRuleSize) {
                                        selectedIndex = a;
                                        minRuleSize = ruleSize;
                                    }
                                }
                                indexEvaluated.add(selectedIndex);
                                listPCIndex[i] = new int[1];
                                listPCTermianl[i] = new double[1][20];
                                listPCIndex[i][0] = selectedIndex;
                            }

                        } else {
                            listPCIndex[i] = new int[1];
                            listPCTermianl[i] = new double[1][20];
                            for (int j = 0; j < indsCharListsMultiTree.length; j++) {
                                if (Arrays.equals(listPCNum.get(i).getKey(), indsCharListsMultiTree[j])) {
                                    indexEvaluated.add(j);
                                    listPCIndex[i][0] = j;
                                    break;
                                }
                            }
                        }
                    }
                    topDuplicatePCNum.add(topDuplicatePcNum);
                    topDuplicatePCNichingNum.add(nicheNumOver1);
                    topExtraEvaluatedNum.add(topExtraEvaluatedIndNum);
                }

                super.evaluatePopulation(state);

                double[][] fitnessesForModel = new double[listPCIndex.length][];
                int[][] indsCharListsMultiTreeEvaluated = new int[listPCIndex.length][40];
                for (int a = 0; a < listPCIndex.length; a++) {
                    indsCharListsMultiTreeEvaluated[a] = indsCharListsMultiTree[listPCIndex[a][0]];
                    if (listPCIndex[a].length == 1) {
                        fitnessesForModel[a] = new double[1];
                        fitnessesForModel[a][0] = state.population.subpops[0].individuals[listPCIndex[a][0]].fitness.fitness();
                    } else {
                        fitnessesForModel[a] = new double[listPCIndex[a].length];
                        for (int b = 0; b < listPCIndex[a].length; b++) {
                            fitnessesForModel[a][b] = state.population.subpops[0].individuals[listPCIndex[a][b]].fitness.fitness();
                        }
                    }
                }

                for (int i = 0; i < fitnessesForModel.length; i++) {
                    if (fitnessesForModel[i][0] == Double.MAX_VALUE) {
                        indsCharListsMultiTreeEvaluated = ArrayUtils.remove(indsCharListsMultiTreeEvaluated, i);
                        fitnessesForModel = ArrayUtils.remove(fitnessesForModel, i);
                        if (state.generation > 0)
                            listPCTermianl = ArrayUtils.remove(listPCTermianl, i);
                        i--;

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
                                if (fitnessesForModel[index].length == 1)
                                    ((Clearable) individual.fitness).surrogateFitness(fitnessesForModel[index][0]);
                                else {
                                    double[] indsTerminalFrequency = getNumTerminalSet(individual);
                                    int selectedIndex = 0;
                                    double dMax = -1;
                                    for (int a = 0; a < listPCTermianl[index].length; a++) {
                                        double correlation = new SpearmansCorrelation().correlation(indsTerminalFrequency, listPCTermianl[index][a]);
                                        if (correlation == 1) {
                                            selectedIndex = a;
                                            break;
                                        }
                                        if (correlation > dMax) {
                                            dMax = correlation;
                                            selectedIndex = a;
                                        }
                                        ((Clearable) individual.fitness).surrogateFitness(fitnessesForModel[index][selectedIndex]);
                                    }
                                }
                                individual.evaluated = false;
                                //individual.fitness.trials = new ArrayList();
                                //individual.fitness.trials.add(individual.fitness.fitness());
                            }

                        } else {
                            ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE);
                        }
                    }
                }

//                        if (state.generation == state.numGenerations - 1) {
                if (((GPRuleEvolutionState)state).totalTime >= ((GPRuleEvolutionState)state).endTime) {
                    writeListSize(listPCSize, duplicatePCNum);
                    writeDiversityToFile(entropyDiversity);
                    writeTopListSize(topPCNum,topDuplicatePCNum,topDuplicatePCNichingNum,topExtraEvaluatedNum);
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
    public void evaluatePopulation(final EvolutionState state, int[][] indsCharListsMultiTree, double[][] fitnessesForModel) {

        int[][] indsCharListsIntermediatePop = phenotypicForSurrogate.muchBetterPhenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //
            {
                Individual individual = state.population.subpops[sub].individuals[i];
                if (indsCharListsMultiTree.length != 0) {
                    //KNN
                    //===============================start==============================================
                    double dMin = Double.MAX_VALUE;
                    int index = 0;
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
                    if (fitnessesForModel[index].length == 1)
                        ((Clearable) individual.fitness).surrogateFitness(fitnessesForModel[index][0]);
                    else {
                        double[] indsTerminalFrequency = getNumTerminalSet(individual);
                        int selectedIndex = 0;
                        double dMax = -1;
                        for (int a = 0; a < listPCTermianl[index].length; a++) {
                            double correlation = new SpearmansCorrelation().correlation(indsTerminalFrequency, listPCTermianl[index][a]);
                            if (correlation == 1) {
                                selectedIndex = a;
                                break;
                            }
                            if (correlation > dMax) {
                                dMax = correlation;
                                selectedIndex = a;
                            }
                            ((Clearable) individual.fitness).surrogateFitness(fitnessesForModel[index][selectedIndex]);
                        }
                    }
                    individual.evaluated = false;

                } else {
                    ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE);
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


    public void setClear(boolean clear) {
        this.clear = clear;
    }

    public static int countSubString(String string, String subString) {
        // 定义一个count来存放字符串出现的次数
        int count = 0;
        // 调用String类的indexOf(String str)方法，返回第一个相同字符串出现的下标
        while (string.indexOf(subString) != -1) {
            count++;
            int index = string.indexOf(subString);
            //将长的那个字符串做截取
            string = string.substring(index + 1);
        }
        return count;
    }

    public static double[] getNumTerminalSet(Individual inds) {
        double[] numTerminalSet = new double[20];
        int countWIQ1, countTIS1, countW1, countNIQ1, countNOR1, countWKR1, countOWT1, countNPT1, countPT1, countMWT1;
        int sumTerminalSet1 = 0;
        int countWIQ2, countTIS2, countW2, countNIQ2, countNOR2, countWKR2, countOWT2, countNPT2, countPT2, countMWT2;
        int sumTerminalSet2 = 0;

        String T1 = StringUtils.substringBetween(inds.toString(), "T0:", "T1:");
        String T2 = StringUtils.substringAfter(inds.toString(), "T1:");

        countWIQ1 = countSubString(T1, "WIQ");
        countTIS1 = countSubString(T1, "TIS");
        countMWT1 = countSubString(T1, "MWT");
        countNIQ1 = countSubString(T1, "NIQ");
        countNOR1 = countSubString(T1, "NOR");
        countWKR1 = countSubString(T1, "WKR");
        countOWT1 = countSubString(T1, "OWT");
        countNPT1 = countSubString(T1, "NPT");
        countPT1 = countSubString(T1, "PT") - countNPT1;
        countW1 = countSubString(T1, "W") - countWIQ1 - countWKR1 - countOWT1 - countMWT1;
        sumTerminalSet1 = countWIQ1 + countTIS1 + countW1 + countNIQ1 + countNOR1 + countWKR1 + countOWT1 + countNPT1 + countPT1 + countMWT1;

        countWIQ2 = countSubString(T2, "WIQ");
        countTIS2 = countSubString(T2, "TIS");
        countMWT2 = countSubString(T2, "MWT");
        countNIQ2 = countSubString(T2, "NIQ");
        countNOR2 = countSubString(T2, "NOR");
        countWKR2 = countSubString(T2, "WKR");
        countOWT2 = countSubString(T2, "OWT");
        countNPT2 = countSubString(T2, "NPT");
        countPT2 = countSubString(T2, "PT") - countNPT2;
        countW2 = countSubString(T2, "W") - countWIQ2 - countWKR2 - countOWT2 - countMWT2;
        sumTerminalSet2 = countWIQ2 + countTIS2 + countW2 + countNIQ2 + countNOR2 + countWKR2 + countOWT2 + countNPT2 + countPT2 + countMWT2;
            /*System.out.println("WIQ有"+countWIQ+"个");
            System.out.println("TIS有"+countTIS+"个");
            System.out.println("W有"+countW+"个");
            System.out.println("NIQ有"+countNIQ+"个");
            System.out.println("NOR有"+countNOR+"个");
            System.out.println("WR有"+countWR+"个");
            System.out.println("OWT有"+countOWT+"个");
            System.out.println("NPT有"+countNPT+"个");
            System.out.println("PT有"+countPT+"个");
            System.out.println("MWT有"+countMWT+"个");
            System.out.println("TermibalSet共有"+sumTerminalSet+"个");*/
        numTerminalSet[0] = (float) countWIQ1 / sumTerminalSet1;
        numTerminalSet[1] = (float) countTIS1 / sumTerminalSet1;
        numTerminalSet[2] = (float) countMWT1 / sumTerminalSet1;
        numTerminalSet[3] = (float) countNIQ1 / sumTerminalSet1;
        numTerminalSet[4] = (float) countNOR1 / sumTerminalSet1;
        numTerminalSet[5] = (float) countWKR1 / sumTerminalSet1;
        numTerminalSet[6] = (float) countOWT1 / sumTerminalSet1;
        numTerminalSet[7] = (float) countNPT1 / sumTerminalSet1;
        numTerminalSet[8] = (float) countPT1 / sumTerminalSet1;
        numTerminalSet[9] = (float) countW1 / sumTerminalSet1;
        numTerminalSet[10] = (float) countWIQ2 / sumTerminalSet2;
        numTerminalSet[11] = (float) countTIS2 / sumTerminalSet2;
        numTerminalSet[12] = (float) countMWT2 / sumTerminalSet2;
        numTerminalSet[13] = (float) countNIQ2 / sumTerminalSet2;
        numTerminalSet[14] = (float) countNOR2 / sumTerminalSet2;
        numTerminalSet[15] = (float) countWKR2 / sumTerminalSet2;
        numTerminalSet[16] = (float) countOWT2 / sumTerminalSet2;
        numTerminalSet[17] = (float) countNPT2 / sumTerminalSet2;
        numTerminalSet[18] = (float) countPT2 / sumTerminalSet2;
        numTerminalSet[19] = (float) countW2 / sumTerminalSet2;

        return numTerminalSet;
    }

    public static void writeListSize(ArrayList<Integer> listPCSize, ArrayList<Integer> duplicatePCNum) {
//        Parameter p;
        // Get the job seed.
//        p = new Parameter("seed").push("" + 0);
//        long jobSeed = parameters.getLongWithDefault(p, null, 0);

        File listSize = new File(out_dir + "/job." + GPRuleEvolutionState.jobSeed + ".listSizeTask.csv");

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(listSize));
            writer.write("generation,listSizeTask0,duplicatePCNum,uniquePCNum");
            writer.newLine();
            for (int i = 0; i < listPCSize.size(); i++) {
                writer.write(i + "," + listPCSize.get(i) + "," + duplicatePCNum.get(i) + "," + (listPCSize.get(i) - duplicatePCNum.get(i)));
                writer.newLine();
            }
            listPCSize.clear();
            duplicatePCNum.clear();
//            writer.write(state.generation-1 + "," + 0 + "," + 0);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeTopListSize(ArrayList<Integer> topPCNum, ArrayList<Integer> topDuplicatePCNum,ArrayList<Integer> topDuplicatePCNichingNum, ArrayList<Integer> topExtraEvaluatedNum) {
//        Parameter p;
        // Get the job seed.
//        p = new Parameter("seed").push("" + 0);
//        long jobSeed = parameters.getLongWithDefault(p, null, 0);

        File listSize = new File(out_dir + "/job." + GPRuleEvolutionState.jobSeed + ".topListPCSizeTask.csv");

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(listSize));
            writer.write("generation,listSizeTask0,duplicatePCNum,uniquePCNum,topDuplicatePCNichingNum,topExtraEvaluatedNum");
            writer.newLine();
            for (int i = 0; i < topPCNum.size(); i++) {
                writer.write(i + "," + topPCNum.get(i) + "," + topDuplicatePCNum.get(i) + "," + (topPCNum.get(i) - topDuplicatePCNum.get(i)) + "," + topDuplicatePCNichingNum.get(i) + "," + topExtraEvaluatedNum.get(i));
                writer.newLine();
            }
            topPCNum.clear();
            topDuplicatePCNum.clear();
            topDuplicatePCNichingNum.clear();
            topExtraEvaluatedNum.clear();
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


}

