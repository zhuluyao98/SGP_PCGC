package yimei.jss.algorithm.multitreeModelSurrogate;

import ec.EvolutionState;
import ec.Individual;
import ec.simple.SimpleEvaluator;
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
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CopyOnWriteArrayList;

import static jdk.nashorn.internal.objects.Global.Infinity;
import static yimei.jss.algorithm.multitreeModelSurrogate.PhenotypicGPRuleEvolutionState.*;
import static yimei.jss.gp.GPRun.out_dir;


/**
 * Created by fzhang on 2019.9.24.
 */
//a class can not extend from more than one class
public class surrogateClearingMultitreeEvaluatorV2 extends SimpleEvaluator {

    int[][] indsCharListsIntermediatePop;
    double[] estimatedFitness;
    double[] distance;
    double minEstimatedFit = 0;
    double maxEstimatedFit = 0;
    double minDistance = 0;
    double maxDistance = 0;

    public ArrayList<Integer> useG2Time = new ArrayList<>();
    public  ArrayList<Integer>  samePCNumArray = new ArrayList<>();

    public ArrayList<Double>  factorDiversity0 = new ArrayList<>();

    public  ArrayList<Integer>  savedExtraInds = new ArrayList<>();

    public ArrayList<Map.Entry<int[], Integer>> listPcClusterFit ;

    public static ArrayList<Integer> listPCSize = new ArrayList<>();
    public static ArrayList<Integer> mostPCDuplicatesNum = new ArrayList<>();

    public ArrayList<double[]> dataSet = new ArrayList<>();   //用于Kmeans的样本数据

    KmeansGroupDuplicates clusterer;

    public int k;  //for kmeans

    double f;

    public ArrayList<Integer> samePCIndex = new ArrayList<>();

    public double[][][] pcIndexKTerminal;
    public double[][][] pcIndexKFitness;

    double[][] pcIndexFitness;
    double[][] pcIndexIndvidualIndex;
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
        }
        else{
            if (nonIntermediatePop) {
                super.evaluatePopulation(state);


                int[][] indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic
                double diversityValue = PopulationUtils.entropy(indsCharListsMultiTree);
                entropyDiversity.add(diversityValue);

                //get the fitness for training model
                double[] fitnessesForModel = new double[state.population.subpops[0].individuals.length];
                for(int ind = 0; ind < state.population.subpops[0].individuals.length; ind++){
                    fitnessesForModel[ind] = state.population.subpops[0].individuals[ind].fitness.fitness();
                }

                double[] fitnessesForModelCopy = fitnessesForModel.clone();
                int[][] indsCharListsMultiTreeCopy = indsCharListsMultiTree.clone();


                // int removeIdx = 0;
                for (int i = 0; i < fitnessesForModelCopy.length; i++) {
                    if (fitnessesForModelCopy[i] == Double.MAX_VALUE) {
                   /* indsCharListsMultiTree = ArrayUtils.remove(indsCharListsMultiTree, i - removeIdx);
                    fitnessesForModel = ArrayUtils.remove(fitnessesForModel, i - removeIdx);*/
                        indsCharListsMultiTreeCopy = ArrayUtils.remove(indsCharListsMultiTreeCopy, i);
                        fitnessesForModelCopy = ArrayUtils.remove(fitnessesForModelCopy, i);
                        i--;
                        // removeIdx++;
                    }
                }
                tempfitnessesForModel = fitnessesForModelCopy;
                tempindsCharListsMultiTree = indsCharListsMultiTreeCopy;
                nonIntermediatePop = false;


                //-------begin------新加的分组引导策略
                //自适应参数计算
                double part1 =0;
                double part2 =0;
                if(state.generation == 0){
                    double maxEntropy = diversityValue;
                    part1 = 1;
                    part2 = 1;
                }
                if(state.generation == 1){
                    part1 = 1;
                    part2 = 1;
                }
                if(state.generation >= 2){
                    double d1 = entropyDiversity.get(state.generation - 2) - entropyDiversity.get(state.generation - 1);
                    double d2 = entropyDiversity.get(state.generation - 1) - entropyDiversity.get(state.generation);
                    part1 = (0.5 + d2/(2*(Math.abs(d1)+Math.abs(d2))));
                    part2 = 1;
                    //part2 = (maxEntropy - entropyDiversity0.get(generation))/maxEntropy;
                }

                f = part1 * part2 * Math.pow(1 - (float) state.generation / state.numGenerations, 2);

                //luyao 8/7/2022 用hashmap记录下每个重复PC的个数
                ConcurrentHashMap<int[], Integer> duplicates = new ConcurrentHashMap<>();
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
                mostPCDuplicatesNum.add(listPcClusterFit.get(0).getValue());

                //luyao 取出k值
                k = state.parameters.getIntWithDefault(new Parameter("num-kmeans"), null, 5);
                int numMax = 0;  //记录下样本中无穷大的数量，避免聚类出现某个簇为空

                pcIndexKTerminal = new double[listPcClusterFit.size()][][];
                pcIndexKFitness = new double[listPcClusterFit.size()][][];

                pcIndexFitness = new double[listPcClusterFit.size()][];
                pcIndexIndvidualIndex = new double[listPcClusterFit.size()][];

                //luyao 4/7/2022 1.根据list的size判断有多少相同PC需要聚类；2.然后聚类
                for (int i=0; i<listPcClusterFit.size() ; i++){

                    if(listPcClusterFit.get(i).getValue() > 1 ){

                        for(int j=0; j<indsCharListsMultiTree.length ; j++){
                            if(Arrays.equals(listPcClusterFit.get(i).getKey(),indsCharListsMultiTree[j]))
                                samePCIndex.add(j);


                        }
                        double[][] dateset = new double[samePCIndex.size()][1];

                        pcIndexFitness[i] = new double[samePCIndex.size()];
                        pcIndexIndvidualIndex[i] = new double[samePCIndex.size()];

                        for(int index=0; index<samePCIndex.size(); index++)
                        {
                            pcIndexFitness[i][index] = state.population.subpops[0].individuals[samePCIndex.get(index)].fitness.fitness();
                            pcIndexIndvidualIndex[i][index] = samePCIndex.get(index);
                        }
                        for (int a=0; a<pcIndexFitness[i].length; a++){
                            dateset[a][0] = pcIndexFitness[i][a];
                            dataSet.add(dateset[a]);
                            if(pcIndexFitness[i][a] == Double.MAX_VALUE)
                                numMax++;}

                        if(dataSet.size() ==  4)
                            k=4;
                        if(dataSet.size() ==  3)
                            k=3;
                        if(dataSet.size() ==  2)
                            k=2;

                        //看有几个独一无二的fitness,如果该值小于K，则将此值设为K
                        double[] tempDateSet = new double[dateset.length];
                        for (int a=0; a<dateset.length; a++)
                            tempDateSet[a] = dateset[a][0];
                        Arrays.sort(tempDateSet);
                        int numUnique = 1;
                        for (int b=1; b<tempDateSet.length; b++){
                            if(tempDateSet[b] != tempDateSet[b-1])
                                numUnique++;
                        }

                        if(dataSet.size() - numMax == 2)
                            k=3;
                        if(dataSet.size() - numMax == 1)
                            k=2;
                        if(dataSet.size() - numMax == 0)
                            k=1;
                        if(dataSet.size() == 2 && numMax == 0)
                            k=2;

                        if(numUnique < k)
                            k = numUnique;

                        clusterer = new KmeansGroupDuplicates(k);
                        clusterer.setDataSet(dataSet);
                        clusterer.execute(state);
                        //可能会变，因为在kmeans中将空的簇删除了
                        k = clusterer.getCluster().size();

                        pcIndexKTerminal[i] = new double[k][20];
                        pcIndexKFitness[i] = new double[k][1];

                        Individual[] clusterPoint =new Individual[k];
                        double pointFit =0;
                        ArrayList<Integer> selectedCenterIndex = new ArrayList<>();
                        double[][] indsTerminalFrequencySet;
                        for (int a=0 ; a<k ; a++){
                            double[] ttt = clusterer.center.get(a);
                            pointFit = getApproaximate(dateset, ttt[0]);
                            pcIndexKFitness[i][a][0] = pointFit;
                            indsTerminalFrequencySet = new double[clusterer.getCluster().get(a).size()][20];
                            for(int b=0; b<clusterer.getCluster().get(a).size(); b++) {
                                //selectedCenterIndex = new int[clusterer.getCluster().get(a).size()];

                                for (int centerIndex = 0; centerIndex < pcIndexFitness[i].length; centerIndex++) {
                                    int judge = 0;
                                    if (pcIndexFitness[i][centerIndex] == clusterer.getCluster().get(a).get(b)[0]) {
                                        for (int c=0; c<selectedCenterIndex.size(); c++){
                                            if(centerIndex == selectedCenterIndex.get(c))
                                            {judge = 1;
                                                break;}
                                            else
                                                continue;}
                                        if(judge == 1)
                                            continue;
                                        selectedCenterIndex.add(centerIndex);
                                        indsTerminalFrequencySet[b] = getNumTerminalSet((Individual) state.population.subpops[0].individuals[(int) pcIndexIndvidualIndex[i][selectedCenterIndex.get(b)]].clone());
                                        break;
                                    }
                                }
                                for(int index=0; index<20; index++){
                                    pcIndexKTerminal[i][a][index] += indsTerminalFrequencySet[b][index]/clusterer.getCluster().get(a).size();
                                }

                            }
                            selectedCenterIndex.clear();}
                        samePCIndex.clear();
                        dataSet.clear();
                        numMax = 0;

                    }

                    //将PC对应一个fitness的信息也放到PCTerminalFit中
                    else {
                        double fit = 0;
                        pcIndexKTerminal[i] = new double[1][20];
                        pcIndexKFitness[i] = new double[1][1];
                        for (int a=0; a<tempindsCharListsMultiTree.length; a++){
                            if(Arrays.equals(tempindsCharListsMultiTree[a],listPcClusterFit.get(i).getKey())){
                                fit = tempfitnessesForModel[a];
                                break;}
                            else
                                fit = Double.MAX_VALUE; //如果在原来构造的surrogate中没找到样本，说明其fitness是Double.Max
                        }
                        pcIndexKFitness[i][0][0] = fit;
                    }
                }

                pcIndexFitness = null;
                pcIndexIndvidualIndex = null;

                System.out.println(listPCSize.size());
                //------end------

                if(state.generation == state.numGenerations-1 ){
                    writeListSize(this.listPCSize,this.mostPCDuplicatesNum);
                    writeDiversityToFile(entropyDiversity);
                }


            } else {
                //multitreeClearingFirst.clearPopulation(state, radius, capacity, phenoCharacterisation);//the individuals do not have fitness or they are evaluated yet
                //but has a fitness there
                this.evaluatePopulation(state, tempindsCharListsMultiTree, tempfitnessesForModel);

                //如果需要fitness加上diversity请取消注释
/*                    estimatedFitness =new double[state.population.subpops[0].individuals.length];
                    for (int i=0; i<state.population.subpops[0].individuals.length; i++){
                        estimatedFitness[i] = state.population.subpops[0].individuals[i].fitness.fitness();}
                    //因为用了分组引导set fitness后，有double.max，所以在归一化时要删掉
                    CopyOnWriteArrayList<Double> tempArray = new CopyOnWriteArrayList<>();
                    for (int i=0; i<estimatedFitness.length; i++){
                        if(estimatedFitness[i]!= Double.MAX_VALUE)
                            tempArray.add(estimatedFitness[i]);
                    }
                    maxEstimatedFit = Collections.max(tempArray);
                    minEstimatedFit = Collections.min(tempArray);
                this.diversityEvaluate(state,indsCharListsIntermediatePop);

                factorDiversity0.add(f);
                for (int i=0; i<state.population.subpops[0].individuals.length; i++){
                    double a = (estimatedFitness[i]-minEstimatedFit)/(maxEstimatedFit-minEstimatedFit);
                    double b = (distance[i]-minDistance)/(maxDistance-minDistance);
                    ((Clearable)state.population.subpops[0].individuals[i].fitness).surrogateFitness((1-f)*a-f*b);
                    //((Clearable)state.population.subpops[0].individuals[i].fitness).surrogateFitness(-distance[i]);
                }*/

                nonIntermediatePop = true;

                if(state.generation == state.numGenerations-2) {
                    writeUseG2Time(state, useG2Time);
                    writeSamePCNum(state,samePCNumArray);
                    writeSavedExtraInds(state, savedExtraInds);
                    //writeFactorDiversityTasks(state,factorDiversity0);
                }


            }
        }

    }

    //==========================================KNN========================================
    public void evaluatePopulation(final EvolutionState state, int[][] indsCharListsMultiTree, double[] fitnessesForModel) {

        indsCharListsIntermediatePop = phenotypicForSurrogate.muchBetterPhenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        //test - begin
       /* HashMap<int[], Integer> duplicatesMidPop = new HashMap<>();
        for(int i = 0; i < (indsCharListsIntermediatePop).length; i++){
            int judge = 0;
            for(int[] k : duplicatesMidPop.keySet()){
                if(Arrays.equals(k,indsCharListsIntermediatePop[i])){
                    judge = 1;
                    int tempValue = duplicatesMidPop.get(k);
                    tempValue++;
                    duplicatesMidPop.replace(k,tempValue);
                    break;
                }
            }
            if(judge == 0)
                duplicatesMidPop.put(indsCharListsIntermediatePop[i],1 );
        }
        ArrayList<HashMap.Entry<int[], Integer>> listPcClusterFitMidPop = new ArrayList<>(duplicatesMidPop.entrySet());
        //通过list对hashmap排序
        Collections.sort(listPcClusterFitMidPop, new Comparator<HashMap.Entry<int[], Integer>>()
        {
            @Override
            public int compare(HashMap.Entry<int[], Integer> o1, HashMap.Entry<int[], Integer> o2)
            {
                //按照value值，从大到小排序
                return o2.getValue() - o1.getValue();
                //按照value值，用compareTo()方法默认是从小到大排序
                //return o1.getValue().compareTo(o2.getValue());
            }
        });*/
        //test - end

        int noUseTime = 0;
        int samePCNum = 0;
        int selectedPCIndex = 0;
        int savedExtraIndsNum = 0;

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //
            {
                Individual individual = state.population.subpops[sub].individuals[i];
                if (indsCharListsMultiTree.length != 0) {
                    //KNN
                    //===============================start==============================================
                    double dMin = Double.MAX_VALUE;
                    if (individual.evaluated == false) {
                        int[] pcIntermediate = indsCharListsIntermediatePop[i];
                        //pcTerminalFit只记录了重复PC的信息
                        int[] nearestPC = new int[40];
                        for (int pcIndex = 0; pcIndex < pcIndexKTerminal.length; pcIndex++) {

                            double d1 = PhenoCharacterisation.distance(pcIntermediate, listPcClusterFit.get(pcIndex).getKey());
                            if (d1 == 0) {
                                selectedPCIndex = pcIndex;
                                samePCNum++;
                                break;
                            }
                            if (d1 < dMin) {
                                dMin = d1;
                                selectedPCIndex = pcIndex;
                            }
                        }

                        double dMax = -Double.MAX_VALUE;
                        double[] indTerminalSet = getNumTerminalSet(individual);
                        int selectedKIndex = 0;
                        double[] tempZero = new double[20];
                        if(Arrays.equals(pcIndexKTerminal[selectedPCIndex][selectedKIndex],tempZero)){
                            ((Clearable)state.population.subpops[0].individuals[i].fitness).surrogateFitness(pcIndexKFitness[selectedPCIndex][selectedKIndex][0]);
                            noUseTime++;
                        }
                        else {
                            for (int kIndex = 0; kIndex < pcIndexKTerminal[selectedPCIndex].length; kIndex++) {

                                double correlation = new SpearmansCorrelation().correlation(indTerminalSet, pcIndexKTerminal[selectedPCIndex][kIndex]);
                                if (correlation == 1) {
                                    selectedKIndex = kIndex;
                                    break;
                                }
                                if (correlation > dMax) {
                                    dMax = correlation;
                                    selectedKIndex = kIndex;
                                }
                            }
                            //如果有多个簇，判断该个体是否在最好的簇内，是的话赋值，不然给无穷大
                            if (pcIndexKTerminal[selectedPCIndex].length > 1) {
                                int bestIndex = 0;
                                for (int index=1; index<pcIndexKTerminal[selectedPCIndex].length; index++){
                                    if(pcIndexKFitness[selectedPCIndex][index][0] < pcIndexKFitness[selectedPCIndex][bestIndex][0])
                                        bestIndex = index;
                                }

                                if (selectedKIndex == bestIndex) {
                                    ((Clearable) state.population.subpops[0].individuals[i].fitness).surrogateFitness(pcIndexKFitness[selectedPCIndex][selectedKIndex][0]);
                                    savedExtraIndsNum ++;
                                } else
                                    ((ClearingMultiObjectiveFitness) state.population.subpops[0].individuals[i].fitness).setFitness();
                            } else {
                                ((Clearable) state.population.subpops[0].individuals[i].fitness).surrogateFitness(pcIndexKFitness[selectedPCIndex][selectedKIndex][0]);
                            }
                        }
                        // judge = 1;
                        //((MultiObjectiveFitness)individual.fitness).surrogateFitness(((PhenotypicGPRuleEvolutionState)state).terminalSetFit.get(((PhenotypicGPRuleEvolutionState)state).clusterCentreTerminal[a][index]));

                    }
                    //如果judge=0,说明此PC是独一无二的，也可能是没有40位一模一样的
/*                        if (judge == 0) {//calculate the fitness based on surrogate model
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
                        }*/
                        /*if (state.generation == 0 || state.generation == 10 || state.generation == 20 || state.generation == 30
                                || state.generation == 40 || state.generation == 49) {
                            pcDistance.add(dMin);
                        }*/
                } else {
                    ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE);
                }
            }
        }


        pcIndexKFitness = null;
        pcIndexKTerminal = null;

        savedExtraInds.add(savedExtraIndsNum);
        useG2Time.add(state.population.subpops[0].individuals.length-noUseTime);
        samePCNumArray.add(samePCNum);
    }

    public void setClear(boolean clear) {
        this.clear = clear;
    }

    public static int countSubString(String string,String subString){
        // 定义一个count来存放字符串出现的次数
        int count = 0 ;
        // 调用String类的indexOf(String str)方法，返回第一个相同字符串出现的下标
        while (string.indexOf(subString) != -1){
            count ++ ;
            int index = string.indexOf(subString);
            //将长的那个字符串做截取
            string = string.substring(index+1);
        }
        return count;
    }

    public static double getApproaximate(double[][] src , double x){
        if(src == null){
            return -1;
        }
        if(src.length == 1){
            return src[0][0];
        }
        if(x == Infinity)
            x = Double.MAX_VALUE;
        double minDifference = Math.abs(src[0][0] - x);
        int minIndex = 0;
        for (int i=0; i<src.length;i++){
            double temp =Math.abs(src[i][0] - x);
            if(temp < minDifference){
                minIndex = i;
                minDifference = temp;
            }
        }
        return src[minIndex][0];
    }

    public static double[] getNumTerminalSet(Individual inds){
        double[] numTerminalSet = new double[20];
        int countWIQ1,countTIS1,countW1,countNIQ1,countNOR1,countWKR1,countOWT1,countNPT1,countPT1,countMWT1;
        int sumTerminalSet1 = 0 ;
        int countWIQ2,countTIS2,countW2,countNIQ2,countNOR2,countWKR2,countOWT2,countNPT2,countPT2,countMWT2;
        int sumTerminalSet2 = 0 ;

        String T1 = StringUtils.substringBetween(inds.toString(),"T0:","T1:");
        String T2 = StringUtils.substringAfter(inds.toString(),"T1:");

        countWIQ1  = countSubString(T1, "WIQ");
        countTIS1  = countSubString(T1, "TIS");
        countMWT1  = countSubString(T1, "MWT");
        countNIQ1  = countSubString(T1, "NIQ");
        countNOR1  = countSubString(T1, "NOR");
        countWKR1  = countSubString(T1, "WKR");
        countOWT1  = countSubString(T1, "OWT");
        countNPT1  = countSubString(T1, "NPT");
        countPT1   = countSubString(T1, "PT")-countNPT1;
        countW1    = countSubString(T1, "W")-countWIQ1-countWKR1-countOWT1-countMWT1;
        sumTerminalSet1 = countWIQ1+countTIS1+countW1+countNIQ1+countNOR1+countWKR1+countOWT1+countNPT1+countPT1+countMWT1;

        countWIQ2  = countSubString(T2, "WIQ");
        countTIS2  = countSubString(T2, "TIS");
        countMWT2  = countSubString(T2, "MWT");
        countNIQ2  = countSubString(T2, "NIQ");
        countNOR2  = countSubString(T2, "NOR");
        countWKR2  = countSubString(T2, "WKR");
        countOWT2  = countSubString(T2, "OWT");
        countNPT2  = countSubString(T2, "NPT");
        countPT2   = countSubString(T2, "PT")-countNPT2;
        countW2    = countSubString(T2, "W")-countWIQ2-countWKR2-countOWT2-countMWT2;
        sumTerminalSet2 = countWIQ2+countTIS2+countW2+countNIQ2+countNOR2+countWKR2+countOWT2+countNPT2+countPT2+countMWT2;
        numTerminalSet[0]=(float)countWIQ1/sumTerminalSet1;
        numTerminalSet[1]=(float)countTIS1/sumTerminalSet1;
        numTerminalSet[2]=(float)countMWT1/sumTerminalSet1;
        numTerminalSet[3]=(float)countNIQ1/sumTerminalSet1;
        numTerminalSet[4]=(float)countNOR1/sumTerminalSet1;
        numTerminalSet[5]=(float)countWKR1/sumTerminalSet1;
        numTerminalSet[6]=(float)countOWT1/sumTerminalSet1;
        numTerminalSet[7]=(float)countNPT1/sumTerminalSet1;
        numTerminalSet[8]=(float)countPT1/sumTerminalSet1;
        numTerminalSet[9]=(float)countW1/sumTerminalSet1;
        numTerminalSet[10]=(float)countWIQ2/sumTerminalSet2;
        numTerminalSet[11]=(float)countTIS2/sumTerminalSet2;
        numTerminalSet[12]=(float)countMWT2/sumTerminalSet2;
        numTerminalSet[13]=(float)countNIQ2/sumTerminalSet2;
        numTerminalSet[14]=(float)countNOR2/sumTerminalSet2;
        numTerminalSet[15]=(float)countWKR2/sumTerminalSet2;
        numTerminalSet[16]=(float)countOWT2/sumTerminalSet2;
        numTerminalSet[17]=(float)countNPT2/sumTerminalSet2;
        numTerminalSet[18]=(float)countPT2/sumTerminalSet2;
        numTerminalSet[19]=(float)countW2/sumTerminalSet2;


        return numTerminalSet;
    }

    public void diversityEvaluate(final EvolutionState state, int[][] indsCharListsMultiTree) {

        distance = new double[state.population.subpops[0].individuals.length];
        for (int x = 0; x < state.population.subpops[0].individuals.length; x++) {

            //计算所有个体之间距离
            for (int otherRows = 0; otherRows < indsCharListsMultiTree.length; otherRows++) {

                distance[x] += PhenoCharacterisation.distance(indsCharListsMultiTree[x], indsCharListsMultiTree[otherRows]);
            }

        }
        maxDistance = Arrays.stream(distance).max().getAsDouble();
        minDistance = Arrays.stream(distance).min().getAsDouble();
    }

    public void writeFactorDiversityTasks(EvolutionState state, CopyOnWriteArrayList<Double> factorDiversity0) {
//        Parameter p;
//        // Get the job seed.
//        p = new Parameter("seed").push("" + 0);
//        long jobSeed = state.parameters.getLongWithDefault(p, null, 0);

        File indSizeForTask = new File(out_dir+"/job." + GPRuleEvolutionState.jobSeed + ".factorDiversityTasks.csv");

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(indSizeForTask));
            writer.write("generation,factorTask0");
            writer.newLine();
            for (int i = 0; i < factorDiversity0.size(); i++) {
                writer.write(i + "," + factorDiversity0.get(i) );
                writer.newLine();
            }
//            writer.write(state.generation-1 + "," + 0 + "," + 0);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public void writeUseG2Time(EvolutionState state, ArrayList<Integer> useG2Time) {
//        Parameter p;
//        // Get the job seed.
//        p = new Parameter("seed").push("" + 0);
//        long jobSeed = state.parameters.getLongWithDefault(p, null, 0);
        File indSizeForTask = new File(out_dir+"/job." + GPRuleEvolutionState.jobSeed + ".useG2Time.csv");

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(indSizeForTask));
            writer.write("generation,useTime0");
            writer.newLine();
            for (int i = 0; i < useG2Time.size(); i++) {
                writer.write(i + "," + useG2Time.get(i) );
                writer.newLine();
            }
//            writer.write(state.generation-1 + "," + 0 + "," + 0);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeSavedExtraInds(EvolutionState state, ArrayList<Integer> savedExtraInds) {
//        Parameter p;
//        // Get the job seed.
//        p = new Parameter("seed").push("" + 0);
//        long jobSeed = state.parameters.getLongWithDefault(p, null, 0);
        File indSizeForTask = new File(out_dir+"/job." + GPRuleEvolutionState.jobSeed + ".useG2Time.csv");

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(indSizeForTask));
            writer.write("generation,savedExtraInds");
            writer.newLine();
            for (int i = 0; i < savedExtraInds.size(); i++) {
                writer.write(i + "," + savedExtraInds.get(i) );
                writer.newLine();
            }
//            writer.write(state.generation-1 + "," + 0 + "," + 0);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeSamePCNum(EvolutionState state, ArrayList<Integer> samePCNumArray) {
//        Parameter p;
//        // Get the job seed.
//        p = new Parameter("seed").push("" + 0);
//        long jobSeed = state.parameters.getLongWithDefault(p, null, 0);

        File indSizeForTask = new File(out_dir+"/job." + GPRuleEvolutionState.jobSeed + ".samePCNum.csv");

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(indSizeForTask));
            writer.write("generation,samePCNum");
            writer.newLine();
            for (int i = 0; i < samePCNumArray.size(); i++) {
                writer.write(i + "," + samePCNumArray.get(i) );
                writer.newLine();
            }
//            writer.write(state.generation-1 + "," + 0 + "," + 0);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


}
