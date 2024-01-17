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
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import yimei.jss.algorithm.elbowSelectedFeatures.IndexPoint;
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
import static yimei.jss.algorithm.multitreeModelSurrogateSample.PhenotypicGPRuleEvolutionState.writeListSize;
import static yimei.jss.gp.GPRun.out_dir;


/**
 * Created by luyao on 2022.8.27.
 * 对同一PC进行niching, 然后每个niche中只对最小rule size的个体,此时GC是将两个RULE在一起
 * Top PC是通过2次knee-point选择的，pre-selection,top PC,非全部个体
 */
//a class can not extend from more than one class
public class surrogateMultitreePCNichingEvaluatorV5 extends SimpleEvaluator {
    public static final String P_RADIUS = "radius";
    public static final String P_CAPACITY = "capacity";

//    public static final String P_PERCENTAGESURRACC = "percentage-surrogateAcc";
//    protected double percentageSurrogateAcc;

    //fzhang 2018.10.9 to get the pre-generation value
    public static final String P_PRE_GENERATIONS = "pre-generations";
    //ArrayList<Double[]> fitnessesForModel = new ArrayList<>();
    double[] tempfitnessesForModel; //for transfer the fitness into the surrogate model to estimate the fitness
    int[][] tempindsCharListsMultiTree;

    int[][] listPCIndex;

    public ArrayList<Double> percentage = new ArrayList<>();

    double correlationCoeffcient;

    double[][][] listPCTermianl;

    public ArrayList<Integer> listPCSize = new ArrayList<>();
    public ArrayList<Integer> duplicatePCNum = new ArrayList<>();

    public ArrayList<Integer> topPCNum = new ArrayList<>();
    public ArrayList<Integer> topDuplicatePCNum = new ArrayList<>();
    public ArrayList<Integer> topDuplicatePCNichingNum = new ArrayList<>();       //top30%重复PC中niche>2的数量
    public ArrayList<Integer> topExtraEvaluatedNum = new ArrayList<>();           //额外评估的个体数
    //public ArrayList<Integer> topFitDistance = new ArrayList<>();              //niche中的fitness之差

    //记录此PC每个niche的索引，第一个永远是被真实评估个体的索引(第几个小生境的个体索引)
    ArrayList<ArrayList<int[]>> pcNicheIndex;

    protected boolean clear = true;
    protected boolean nonIntermediatePop = true;

    protected double radius;
    protected int capacity;

    public static PhenoCharacterisation[] phenoCharacterisation;


    public ArrayList<Double> entropyDiversity = new ArrayList<>();

    public static ArrayList<Double> entropyGCDiversity = new ArrayList<>();

    public static ArrayList<Double> correlationValue = new ArrayList<>();

    public static ArrayList<Double> co = new ArrayList<>();

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

                    //-------------计算相关性阈值-----------------------
                   /* ArrayList<ArrayList<Integer>> DuplicatedPCIndex = new ArrayList<>();
                    for (int i = 0; i < listPCNum.size(); i++) {
                        if (listPCNum.get(i).getValue() > 1) {
                            ArrayList<Integer> OnePCIndex = new ArrayList<>();
                            for (int j = 0; j < indsCharListsMultiTree.length; j++) {
                                if (Arrays.equals(listPCNum.get(i).getKey(), indsCharListsMultiTree[j]))
                                    OnePCIndex.add(j);
                                if (OnePCIndex.size() == listPCNum.get(i).getValue())
                                    break;
                            }
                            DuplicatedPCIndex.add(OnePCIndex);
                        }
                    }
                    int[][] pcIndex1 = new int[DuplicatedPCIndex.size()][];
                    for (int a=0; a<DuplicatedPCIndex.size(); a++){
                        pcIndex1[a] = new int[DuplicatedPCIndex.get(a).size()];
                        for (int b=0; b<DuplicatedPCIndex.get(a).size(); b++){
                            pcIndex1[a][b] = DuplicatedPCIndex.get(a).get(b);
                        }
                    }
                    double sumCorrelation = 0;
                    int num = 0;
                    double BaseSequencingRuleSize = 0;
                    double BaseRoutingRuleSize = 0;
                    double CounterSequencingRuleSize = 0;
                    double CounterRoutingRuleSize = 0;
                    for (int a=0; a<pcIndex1.length; a++){
                        Individual indsBase = state.population.subpops[0].individuals[pcIndex1[a][0]];
                        BaseSequencingRuleSize = ((GPIndividual) indsBase).trees[0].child.numNodes(GPNode.NODESEARCH_ALL);
                        BaseRoutingRuleSize = ((GPIndividual) indsBase).trees[1].child.numNodes(GPNode.NODESEARCH_ALL);
                        double[] indsTerminalFrequencySetRouBase = getNumTerminalSet(indsBase,"rou");
                        double[] indsTerminalFrequencySetSeqBase = getNumTerminalSet(indsBase,"seq");
//                                    NaturalRanking temp = new NaturalRanking();
//                                    indsTerminalFrequencySetSeqBase = temp.rank(indsTerminalFrequencySetSeqBase);
//                                    indsTerminalFrequencySetRouBase = temp.rank(indsTerminalFrequencySetRouBase);
                        for(int b=1; b<pcIndex1[a].length; b++){
                            Individual indsCounter = state.population.subpops[0].individuals[pcIndex1[a][b]];
                            CounterSequencingRuleSize = ((GPIndividual) indsCounter).trees[0].child.numNodes(GPNode.NODESEARCH_ALL);
                            CounterRoutingRuleSize = ((GPIndividual) indsCounter).trees[1].child.numNodes(GPNode.NODESEARCH_ALL);
                            double[] indsTerminalFrequencySetRouCounter = getNumTerminalSet(indsCounter,"rou");
                            double[] indsTerminalFrequencySetSeqCounter = getNumTerminalSet(indsCounter,"seq");
                            double correlationRou = new SpearmansCorrelation().correlation(indsTerminalFrequencySetRouBase, indsTerminalFrequencySetRouCounter);
                            double correlationSeq = new SpearmansCorrelation().correlation(indsTerminalFrequencySetSeqBase, indsTerminalFrequencySetSeqCounter);
                            sumCorrelation += ((BaseRoutingRuleSize+CounterRoutingRuleSize)/(BaseRoutingRuleSize+CounterRoutingRuleSize+BaseSequencingRuleSize+CounterSequencingRuleSize))*correlationRou +
                                    ((BaseSequencingRuleSize+CounterSequencingRuleSize)/(BaseRoutingRuleSize+CounterRoutingRuleSize+BaseSequencingRuleSize+CounterSequencingRuleSize))*correlationSeq;
                            num++;
                        }
                    }

//                    if(state.generation%50 ==0)
//                        System.out.println("sumRou = " + sumRoutingRuleSize + " , sumSeq = " + sumSequencingRuleSize);

                    correlationCoeffcient = sumCorrelation/num;
                    correlationValue.add(correlationCoeffcient);*/
//                    System.out.println(correlationCoeffcient);

                    //-----------------end-------------------------

                    //--------------------begin------------------通过熵的计算求得基因型的多样性
                    double[][] pcRuleSizeGC = new double[state.population.subpops[0].individuals.length][20]; //每个PC取rule size最小的个体的基因型数组
                    NaturalRanking rank = new NaturalRanking();

                    for (int i = 0; i < state.population.subpops[0].individuals.length; i++){
                        pcRuleSizeGC[i] = rank.rank(getNumTerminalSet(state.population.subpops[0].individuals[i]));
                    }
                    double gcEntropy = PopulationUtils.entropy(pcRuleSizeGC);
                    entropyGCDiversity.add(gcEntropy);

/*                    for (int i = 0; i < listPCNum.size(); i++) {
                        if (listPCNum.get(i).getValue() > 1) {
                            ArrayList<Integer> samePCIndex = new ArrayList<>();
                            for (int j = 0; j < indsCharListsMultiTree.length; j++) {
                                if (Arrays.equals(listPCNum.get(i).getKey(), indsCharListsMultiTree[j]))
                                    samePCIndex.add(j);
                            }

                            //每个PC随机选一个进行评估
//                            int random = state.random[0].nextInt(samePCIndex.size());
//                            indexEvaluated.add(samePCIndex.get(random));

                            //每个PC选一个rule size最小的个体进行评估
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
                            pcRuleSizeGC[i] = rank.rank(getNumTerminalSet(state.population.subpops[0].individuals[selectedIndex]));
                        } else {
                            for (int j = 0; j < indsCharListsMultiTree.length; j++) {
                                if (Arrays.equals(listPCNum.get(i).getKey(), indsCharListsMultiTree[j])) {
                                    pcRuleSizeGC[i] = rank.rank(getNumTerminalSet(state.population.subpops[0].individuals[j]));
                                    break;
                                }
                            }
                        }
                    }*//*
                    double gcEntropy = PopulationUtils.entropy(pcRuleSizeGC);
                    entropyGCDiversity.add(gcEntropy);
//                    System.out.println(gcEntropy);

                    correlationCoeffcient = 0.9 - 0.4*(gcEntropy - 4.5/1.5);
                    if(state.generation<40){
                        correlationCoeffcient = 0.5;
                    }
                    else if (state.generation < 70)
                        correlationCoeffcient = 0.6;
                    else if (state.generation < 85)
                    //--------------------end-------------------------*/

                    pcNicheIndex = new ArrayList<>();

                    //------------------begin---------------使用knee-point筛选top individuals

                    Individual[] individuals = state.population.subpops[0].individuals;
                    ArrayList<IndexPoint> points = new ArrayList<IndexPoint>();
                    ArrayList<Integer> numSelIndividual = new ArrayList<>();

                    for(int indIndex = 0; indIndex < state.population.subpops[0].individuals.length; indIndex++){
                        //System.out.println(individuals[indIndex].fitness.fitness());
                        if(individuals[indIndex].fitness.fitness() != Double.POSITIVE_INFINITY & individuals[indIndex].fitness.fitness() != Double.MAX_VALUE)
                            points.add(new IndexPoint(new double[]{indIndex,individuals[indIndex].fitness.fitness()}));
                    }

              /* System.out.println("after delete inifinity fitness of individuals \n");
               for(int ind = 0; ind < points.size(); ind++){
                   System.out.println(points.get(ind).position[1]);
               }*/


                    int kneePoint = -1;
                    for (int num=0; num<1; num++) {
                        if(num>0)
                            points.subList(kneePoint, points.size()).clear();
                        IndexPoint p1 = points.get(0);  //the min Point
                        IndexPoint p2 = points.get(points.size() - 1); //the max point

                        //calculate the line factors
                        double aa = p2.position[1] - p1.position[1];
                        double bb = p1.position[0] - p2.position[0];
                        double cc = p2.position[0] * p1.position[1] - p1.position[0] * p2.position[1];
              /* System.out.println(a*p1.position[0]+b*p1.position[1]+c);
               System.out.println(a*p2.position[0]+b*p2.position[1]+c);*/

                        double maxDistance = 0.0;


                        for (int idxPoints = 1; idxPoints < points.size(); idxPoints++) {
                            IndexPoint p = points.get(idxPoints);
                            double distance = Math.abs(aa * p.position[0] + bb * p.position[1] + cc);
                            if (distance > maxDistance) {
                                maxDistance = distance;
                                kneePoint = idxPoints;
                            }
                        }
                    }
                    //System.out.println("index of knee point: " + kneePoint);



                    //----------------------end--------------------------------

                    //只有在第1代往后才记录Top30%的个体
                    int pcLastTopIndex = 0;
                    //记录下Top30%的PC最后一个在list中的索引
                    int topPcNum = (int) (state.population.subpops[0].individuals.length * 0.3);
//                    int topPcNum = kneePoint;
                    int topLastIndex = topPcNum - 1;
                    for (int a = 0; a < listPCNum.size(); a++) {
                        if (Arrays.equals(indsCharListsMultiTree[topLastIndex], listPCNum.get(a).getKey())) {
                            pcLastTopIndex = a;
                            break;
                        }
                    }
//                    System.out.println((double) pcLastTopIndex/listPCNum.size());

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

//                            if (i <= 500) {
                            if (i <= pcLastTopIndex) {        //说明是Top30%的PC,需要niching
                                topDuplicatePcNum++;

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



                                ArrayList<int[]> nicheIndex = new ArrayList<>();
                                for (int a = 0; a < pcIndex.length; a++) {
                                    ArrayList<Integer> nicheIndexTemp = new ArrayList<>();
                                    Individual inds = state.population.subpops[0].individuals[pcIndex[a]];
                                    double[] indsTerminalFrequencySet = getNumTerminalSet(inds);
                                    if (((ClearingMultiObjectiveFitness) inds.fitness).isInNiche()) {
                                        continue;
                                    } else {

                                        nicheIndexTemp.add(pcIndex[a]);
                                        indexEvaluated.add(pcIndex[a]);
                                        ((ClearingMultiObjectiveFitness) inds.fitness).inNiche();
                                    }
                                    for (int b = 0; b < pcIndex.length; b++) {
                                        Individual inds2 = state.population.subpops[0].individuals[pcIndex[b]];
                                        if (((ClearingMultiObjectiveFitness) inds2.fitness).isInNiche()) {
                                            continue;
                                        }
                                        double[] indsTerminalFrequencySet2 = getNumTerminalSet(inds2);
                                        double correlation = new SpearmansCorrelation().correlation(indsTerminalFrequencySet, indsTerminalFrequencySet2);
                                        //---------------------


                                        if (correlation >= correlationCoeffcient) {
                                            nicheIndexTemp.add(pcIndex[b]);
                                            ((ClearingMultiObjectiveFitness) inds2.fitness).inNiche();
                                        }
                                    }
                                    int[] oneNicheIndex = new int[nicheIndexTemp.size()];     //记录一个生境内的个体索引，第一个是真实评估的个体索引
                                    for (int index=0; index<nicheIndexTemp.size(); index++){
                                        oneNicheIndex[index] = nicheIndexTemp.get(index);
                                    }
                                    nicheIndex.add(oneNicheIndex);
                                }

                                pcNicheIndex.add(nicheIndex);

                                if (nicheIndex.size() > 1) {
                                    nicheNumOver1++;
                                    topExtraEvaluatedIndNum += (nicheIndex.size() - 1);
                                }


                            } else {                         //后70%的PC,只需要取出rule size最小的Index
                                ArrayList<int[]> nicheIndex = new ArrayList<>();
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
                                int[] oneNicheIncex = new int[samePCIndex.size()];
                                oneNicheIncex[0] = selectedIndex;
                                samePCIndex.remove(Integer.valueOf(selectedIndex));
                                for (int a=0; a<samePCIndex.size(); a++){
                                    oneNicheIncex[a+1] = samePCIndex.get(a);
                                }
                                nicheIndex.add(oneNicheIncex);
                                pcNicheIndex.add(nicheIndex);
                                /*listPCIndex[i] = new int[1];
                                listPCTermianl[i] = new double[1][20];
                                listPCIndex[i][0] = selectedIndex;*/
                            }

                        } else {              //一个PC对应一个fitness
                            /*listPCIndex[i] = new int[1];
                            listPCTermianl[i] = new double[1][20];*/
                            for (int j = 0; j < indsCharListsMultiTree.length; j++) {
                                if (Arrays.equals(listPCNum.get(i).getKey(), indsCharListsMultiTree[j])) {
                                    indexEvaluated.add(j);
//                                    listPCIndex[i][0] = j;
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

                //-----------------begin将同niche的fitness设置成相同的fitness---------------------
                if (state.generation > 0) {
                    for (int a=0; a<pcNicheIndex.size(); a++){
                        for (int b=0; b<pcNicheIndex.get(a).size(); b++){
                            Individual indRealEvaluation = state.population.subpops[0].individuals[pcNicheIndex.get(a).get(b)[0]];
                            for (int c=1; c<pcNicheIndex.get(a).get(b).length; c++){
                                Individual ind = state.population.subpops[0].individuals[pcNicheIndex.get(a).get(b)[c]];
                                ((Clearable) ind.fitness).surrogateFitness(indRealEvaluation.fitness.fitness());
                                ind.evaluated = false;
                            }
                        }
                    }
                }

                //------------------end---------------------------------------
                double[] fitnessesForModel = new double[state.population.subpops[0].individuals.length];
                for(int ind = 0; ind < state.population.subpops[0].individuals.length; ind++){
                    fitnessesForModel[ind] = state.population.subpops[0].individuals[ind].fitness.fitness();
                }

                int [][] indsCharListsMultiTreeCopy =  indsCharListsMultiTree.clone();

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

                //-----------------------------begin*(输出同PC各个niche中心的fitness)--------------------------------


                /*int sum = 0;
                int sumbase = 0;
                for (int a=0; a<tempfitnessesForModel.length; a++){
                    if(tempfitnessesForModel[a].length > 1){
                        sumbase ++;

                        ArrayList<Double> fitDis = new ArrayList<>();
                        ArrayList<Double> corDis = new ArrayList<>();
                        for (int b=0; b<tempfitnessesForModel[a].length; b++){
                            //System.out.print(tempfitnessesForModel[a][b] + " , ");
                            double correlation = new SpearmansCorrelation().correlation(listPCTermianl[a][0], listPCTermianl[a][b]);
                            corDis.add(correlation);
                            fitDis.add(Math.abs(tempfitnessesForModel[a][b]-tempfitnessesForModel[a][0]));
                        }
                        for (int b=0; b<tempfitnessesForModel[a].length; b++){
                            if(tempfitnessesForModel[a][b] < tempfitnessesForModel[a][0]){
                                sum++;
                                break;}
                        }
                        System.out.print("第"+ a +"个PC的"+"相关性是  ");
                        for (int index=0; index<corDis.size(); index++){
                            System.out.print(corDis.get(index)+ ",");
                        }
                        System.out.print(" ---- "+"第"+ a +"个PC的"+"适应度之差是  " );
                        for (int index=0; index<fitDis.size(); index++){
                            System.out.print(fitDis.get(index)+ ",");
                        }
                        System.out.println();
                    }
                }
                double p = 0;
                if(sumbase == 0){
                    p = -1;
                }
                else {
                    p = (double) sum/sumbase;}
//                if(state.generation > 0)
//                    System.out.println(p);
                percentage.add(p);*/


                //-----------------------------end--------------------------------
                if (((GPRuleEvolutionState)state).totalTime >= ((GPRuleEvolutionState)state).endTime) {
//                if (state.generation == state.numGenerations - 1) {
                    writeListSize(listPCSize, duplicatePCNum);
                    writeDiversityToFile(entropyDiversity);
                    writeGCDiversityToFile(entropyGCDiversity);
                    writeCorrelationToFile(correlationValue);
                    writeTopListSize(topPCNum,topDuplicatePCNum,topDuplicatePCNichingNum,topExtraEvaluatedNum);
//                    writePercentageToFile(percentage);
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
                Individual individual = state.population.subpops[sub].individuals[i];
                if(indsCharListsMultiTree.length != 0){
                    //KNN
                    //===============================start==============================================
                    double dMin  = Double.MAX_VALUE;
                    int index = 0;
                    if(individual.evaluated == false){
                        int[] pcIntermediate = indsCharListsIntermediatePop[i];
                        //calculate the fitness based on surrogate model
                        for(int pc = 0; pc < indsCharListsMultiTree.length; pc++){
                            int[] pcModel = indsCharListsMultiTree[pc];
                            double d  = PhenoCharacterisation.distance(pcIntermediate, pcModel);
                            if(d == 0){
                                dMin = d;
                                index = pc;
                                break;
                            }

                            if (d < dMin){
                                dMin = d;
                                index = pc;
                            }
                        }
                        ((Clearable)individual.fitness).surrogateFitness(fitnessesForModel[index]);
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

 /*   public static double[] getNumTerminalSet(Individual inds, String rule) {
        double[] numTerminalSet = new double[10];
        int countWIQ1, countTIS1, countW1, countNIQ1, countNOR1, countWKR1, countOWT1, countNPT1, countPT1, countMWT1;
        int sumTerminalSet1 = 0;
        String T1 = "";

        if(rule == "seq")
            T1 = StringUtils.substringBetween(inds.toString(), "T0:", "T1:");
        if(rule == "rou")
            T1 = StringUtils.substringAfter(inds.toString(), "T1:");

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

        return numTerminalSet;
    }

    public void writeListSize(ArrayList<Integer> listPCSize, ArrayList<Integer> duplicatePCNum) {
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
    }*/

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

    public static double[] getNumTerminalSet(Individual inds, String rule) {
        double[] numTerminalSet = new double[10];
        int countWIQ1, countTIS1, countW1, countNIQ1, countNOR1, countWKR1, countOWT1, countNPT1, countPT1, countMWT1;
        int sumTerminalSet1 = 0;
        String T1 = "";

        if(rule == "seq")
            T1 = StringUtils.substringBetween(inds.toString(), "T0:", "T1:");
        if(rule == "rou")
            T1 = StringUtils.substringAfter(inds.toString(), "T1:");

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

        return numTerminalSet;
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
    public void writeDiversityToFile(ArrayList<Double> entropyDiversity) {
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

    public void writeGCDiversityToFile(ArrayList<Double> entropyDiversity) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + GPRuleEvolutionState.jobSeed + ".diversityGC.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Gen,diversitySubpop0GC");
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

    public void writeCorrelationToFile(ArrayList<Double> correlationValue) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + GPRuleEvolutionState.jobSeed + ".correlationValue.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Gen,correlationValue");
            writer.newLine();
            for (int i = 0; i < correlationValue.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + correlationValue.get(i) + "\n");
            }
            correlationValue.clear();
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writePercentageToFile(ArrayList<Double> percentage) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + GPRuleEvolutionState.jobSeed + ".percentage.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Gen,percentage");
            writer.newLine();
            for (int i = 0; i < percentage.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + percentage.get(i) + "\n");
            }
            percentage.clear();
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }




}

