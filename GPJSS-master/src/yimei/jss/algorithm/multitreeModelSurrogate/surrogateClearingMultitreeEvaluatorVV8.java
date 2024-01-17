package yimei.jss.algorithm.multitreeModelSurrogate;

import ec.EvolutionState;
import ec.Individual;
import ec.gp.GPIndividual;
import ec.gp.GPNode;
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

import static yimei.jss.algorithm.multitreeModelSurrogate.PhenotypicGPRuleEvolutionState.*;
import static yimei.jss.gp.GPRun.out_dir;


/**
 * Created by luyao 2022/7/25  pre-selection仍是原surrogate ， 但是排重中根据策略保留最好的那个或者rule size最小的
 */
//a class can not extend from more than one class
public class surrogateClearingMultitreeEvaluatorVV8 extends SimpleEvaluator {

    public ArrayList<ArrayList<Individual>> reservedIndsV8 = new ArrayList<>();

    public static double[][] duplicatedPCMinMaxFit;

    ArrayList<int[]> savedAndDeleteNum = new ArrayList<>();

    public static int[][] indsCharListsIntermediatePop;

    public ArrayList<Integer> useG2Time = new ArrayList<>();
    public  ArrayList<Integer>  samePCNumArray = new ArrayList<>();

    public ArrayList<Map.Entry<int[], Integer>> listPcClusterFit ;

    public static ArrayList<Integer> listPCSize = new ArrayList<>();
    public static ArrayList<Integer> mostPCDuplicatesNum = new ArrayList<>();


    public int k;  //for kmeans

    public ArrayList<Integer> samePCIndex = new ArrayList<>();

    public double[][][] pcIndexKTerminal;
    public double[][][] pcIndexKFitness;

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

                //-----begin----- luyao 8/7/2022 用hashmap记录下每个重复PC的个数
                HashMap<int[], Integer> duplicates = new HashMap<>();
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
                //-----end-----  用hashmap记录下每个重复PC的个数

                //-----begin----- 通过list对hashmap排序
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
                //-----end----- 通过list对hashmap排序

                //-----begin-----
                pcIndexKTerminal = new double[listPcClusterFit.size()][][];
                pcIndexKFitness = new double[listPcClusterFit.size()][][];

                int duplicatedPCNum = 0;
                for (int i=0; i<listPcClusterFit.size() ; i++){
                    if(listPcClusterFit.get(i).getValue() > 1 ){
                        duplicatedPCNum++;
                    }
                }
                duplicatedPCMinMaxFit = new double[duplicatedPCNum][2];

                //luyao 4/7/2022 1.根据list的size判断有多少PC对应多个fitness
                for (int i=0; i<listPcClusterFit.size() ; i++){

                    if(listPcClusterFit.get(i).getValue() > 1 ){

                        for(int j=0; j<indsCharListsMultiTree.length ; j++){
                            if(Arrays.equals(listPcClusterFit.get(i).getKey(),indsCharListsMultiTree[j]))
                                samePCIndex.add(j);
                        }

                        k = samePCIndex.size();

                        pcIndexKTerminal[i] = new double[k][20];
                        pcIndexKFitness[i] = new double[k][1];

                        for(int index=0; index<k; index++) {
                            pcIndexKFitness[i][index][0] = state.population.subpops[0].individuals[samePCIndex.get(index)].fitness.fitness();
                            pcIndexKTerminal[i][index] = getNumTerminalSet((Individual) state.population.subpops[0].individuals[samePCIndex.get(index)].clone());
                        }
                        samePCIndex.clear();

                        //-----begin----- 记录当前PC对应最小与最大的fitness
                        double maxFit = pcIndexKFitness[i][0][0];
                        double minFit = pcIndexKFitness[i][0][0];
                        for (int n=0; n<pcIndexKFitness[i].length; n++){
                            if(pcIndexKFitness[i][n][0] > maxFit)
                                maxFit = pcIndexKFitness[i][n][0];
                            if(pcIndexKFitness[i][n][0] < minFit)
                                minFit = pcIndexKFitness[i][n][0];
                        }
                        duplicatedPCMinMaxFit[i][0] = minFit;
                        duplicatedPCMinMaxFit[i][1] = maxFit;

                        //-----end-----

                    }
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


                if(state.generation == state.numGenerations-1 ){
                    writeListSize(this.listPCSize,this.mostPCDuplicatesNum);
                    writeDiversityToFile(entropyDiversity);

                }
            } else {
                muchBetterClear(state);
                //but has a fitness there
                this.evaluatePopulation(state, tempindsCharListsMultiTree, tempfitnessesForModel);

                nonIntermediatePop = true;

                //输出记录的信息到csv
                if(state.generation == state.numGenerations-2) {
//                    writeUseG2Time(state, useG2Time);
                    writeSamePCNum(state,samePCNumArray);
                    writeClearSaveAndDeleteNum(state,savedAndDeleteNum);
                    writeReservedIndsAverageFitnessAndRuleSize(state,reservedIndsV8);
//                    writeFactorDiversityTasks(state,factorDiversity0);
                }
            }
        }
    }

    private void muchBetterClear(final EvolutionState state) {
        indsCharListsIntermediatePop = phenotypicForSurrogate.muchBetterPhenotypicPopulation(state, phenoCharacterisation);
        HashMap<int[], Integer> duplicatesMidPop = new HashMap<>();
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
        ArrayList<Map.Entry<int[], Integer>> listPcClusterFitMidPop = new ArrayList<>(duplicatesMidPop.entrySet());
        //-----end-----  用hashmap记录下每个重复PC的个数

        //-----begin----- 通过list对hashmap排序
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
        });

        ArrayList<Individual> savedIndsCurrentGen = new ArrayList<>();
        int basedRuleSizeSaveNum = 0; int basedRuleSizeDeleteNum = 0; int basedGuidanceSaveNum = 0; int basedGuidanceDeleteNum = 0;
        for (int i=0; i<listPcClusterFitMidPop.size() ; i++){

            if(listPcClusterFitMidPop.get(i).getValue() > 1 ){

                int selectedPCIndex = 0;
                double Dmin  = Double.MAX_VALUE;

                ArrayList<Integer> samePCIndexMidPop = new ArrayList<>();
                for(int j=0; j<indsCharListsIntermediatePop.length ; j++){
                    if(Arrays.equals(listPcClusterFitMidPop.get(i).getKey(),indsCharListsIntermediatePop[j]))
                        samePCIndexMidPop.add(j);
                }

                for (int pcIndex=0; pcIndex<listPcClusterFit.size(); pcIndex++){
                    double distance = PhenoCharacterisation.distance(listPcClusterFit.get(pcIndex).getKey(),listPcClusterFitMidPop.get(i).getKey());
                    if(distance == 0){
                        selectedPCIndex = pcIndex;
                        break;
                    }
                    if(distance < Dmin){
                        selectedPCIndex = pcIndex;
                        Dmin = distance;
                    }
                }

                //如果该PC只对应一个fitness,则选择rule size最小的那个
                int minRuleSize = 1000;
                int selectedIndex = 0;

                if(listPcClusterFit.get(selectedPCIndex).getValue() == 1){
                    for (int a:samePCIndexMidPop){
                        Individual inds = state.population.subpops[0].individuals[a];
                        int ruleSize = ((GPIndividual)inds).trees[0].child.numNodes(GPNode.NODESEARCH_ALL) + ((GPIndividual) inds).trees[1].child.numNodes(GPNode.NODESEARCH_ALL);
                        if(ruleSize < minRuleSize){
                            selectedIndex = a;
                            minRuleSize = ruleSize;
                        }
                    }
                    for (int index:samePCIndexMidPop){
                        Individual inds = state.population.subpops[0].individuals[index];
                        if(index == selectedIndex){
                            ((Clearable) inds.fitness).surrogateFitness(pcIndexKFitness[selectedPCIndex][0][0]);
                            inds.evaluated = true;
                            basedRuleSizeSaveNum++;
                            savedIndsCurrentGen.add(inds);
                        }
                        else {
                            ((ClearingMultiObjectiveFitness) inds.fitness).setFitness();
                            inds.evaluated = true;
                            basedRuleSizeDeleteNum++;
                        }
                    }
                }
                else {
                    //找出最佳的fitness索引
                    int bestKIndex = 0;
                    int selectedKIndex = 0;
                    for (int k=1; k<pcIndexKFitness[selectedPCIndex].length; k++){
                        if(pcIndexKFitness[selectedPCIndex][k][0]<pcIndexKFitness[selectedPCIndex][bestKIndex][0])
                            bestKIndex = k;
                    }
                    double dMax = -Double.MAX_VALUE;
                    for (int a:samePCIndexMidPop){
                        Individual inds = state.population.subpops[0].individuals[a];
                        double[] indsTerminalFrequencySet = getNumTerminalSet(inds);
                        double correlation = new SpearmansCorrelation().correlation(indsTerminalFrequencySet, pcIndexKTerminal[selectedPCIndex][bestKIndex]);
                        if (correlation == 1) {
                            selectedKIndex = a;
                            break;
                        }
                        if (correlation > dMax) {
                            dMax = correlation;
                            selectedKIndex = a;
                        }
                    }
                    for (int index:samePCIndexMidPop){
                        Individual inds = state.population.subpops[0].individuals[index];
                        if(index == selectedKIndex){
                            ((Clearable) inds.fitness).surrogateFitness(pcIndexKFitness[selectedPCIndex][bestKIndex][0]);
                            inds.evaluated = true;
                            savedIndsCurrentGen.add(inds);
                            basedGuidanceSaveNum++;
                        }
                        else {
                            ((ClearingMultiObjectiveFitness) inds.fitness).setFitness();
                            inds.evaluated = true;
                            basedGuidanceDeleteNum++;
                        }
                    }

                }

                samePCIndexMidPop.clear();
            }

        }
        reservedIndsV8.add(savedIndsCurrentGen);
        int[] saveAndDeleteNum = new int[]{basedRuleSizeSaveNum,basedRuleSizeDeleteNum,basedGuidanceSaveNum,basedGuidanceDeleteNum};
        savedAndDeleteNum.add(saveAndDeleteNum);

    }

    //==========================================KNN========================================
    public void evaluatePopulation(final EvolutionState state, int[][] indsCharListsMultiTree, double[] fitnessesForModel) {

        int samePCNum = 0;

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
                                samePCNum ++;
                                break;

                            }

                            if (d < dMin){
                                dMin = d;
                                index = pc;
                            }
                        }
                        ((Clearable)individual.fitness).surrogateFitness(fitnessesForModel[index]);

                    }
                    if (state.generation == 0 || state.generation == 10 || state.generation == 20 || state.generation == 30
                            || state.generation == 40 || state.generation == 49) {
                        pcDistance.add(dMin);
                    }
                }
                else{
                    ((Clearable)individual.fitness).surrogateFitness(Double.MAX_VALUE);
                }
            }
        }
        samePCNumArray.add(samePCNum);
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

    public static void writeReservedIndsAverageFitnessAndRuleSize(EvolutionState state, ArrayList<ArrayList<Individual>> reservedInds) {
        File reservedIndsAverageFit = new File(out_dir+"/job." + GPRuleEvolutionState.jobSeed + ".reservedIndsAverageFitAndRuleSize.csv");

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(reservedIndsAverageFit));
            writer.write("generation,reservedIndsAverageFit,maxNum,AverageSeqRuleSize,AverageRouRuleSize,AverageRuleSize");
            writer.newLine();
            ArrayList <Double> averageFitness = new ArrayList<>();
            ArrayList <Integer> averageSeqRuleSize = new ArrayList<>();
            ArrayList <Integer> averageRouRuleSize = new ArrayList<>();
            ArrayList <Integer> averageRuleSize = new ArrayList<>();
            ArrayList <Integer> maxNum = new ArrayList<>();
            for (int a=0; a<reservedInds.size(); a++){

                double sumFitness=0;
                int max = 0;
                int SeqSizeTree0 = 0;
                int RouSizeTree1 = 0;

                for (int b=0; b<reservedInds.get(a).size(); b++){
                    if(reservedInds.get(a).get(b).fitness.fitness() == Double.MAX_VALUE){
                        max++;
                        continue;
                    }
                    sumFitness +=reservedInds.get(a).get(b).fitness.fitness();

                    GPIndividual indi = (GPIndividual)reservedInds.get(a).get(b) ;
                    SeqSizeTree0 += indi.trees[0].child.numNodes(GPNode.NODESEARCH_ALL);
                    RouSizeTree1 += indi.trees[1].child.numNodes(GPNode.NODESEARCH_ALL);

                }
                averageFitness.add(sumFitness/(reservedInds.get(a).size()-max));
                maxNum.add(max);

                averageSeqRuleSize.add(SeqSizeTree0 / reservedInds.get(a).size());
                averageRouRuleSize.add(RouSizeTree1 / reservedInds.get(a).size());
                averageRuleSize.add((SeqSizeTree0+RouSizeTree1)/reservedInds.get(a).size());

            }
            for (int i = 0; i < averageFitness.size(); i++) {
                writer.write(i + "," + averageFitness.get(i) + "," + maxNum.get(i) + "," + averageSeqRuleSize.get(i)
                        + "," + averageRouRuleSize.get(i) + "," + averageRuleSize.get(i));
                writer.newLine();
            }
//            writer.write(state.generation-1 + "," + 0 + "," + 0);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeClearSaveAndDeleteNum(EvolutionState state, ArrayList<int[]> savedAndDeleteNum) {

        File ClearSaveAndDeleteNum = new File(out_dir+"/job." + GPRuleEvolutionState.jobSeed + ".clearSaveAndDeleteNum.csv");

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(ClearSaveAndDeleteNum));
            writer.write("generation,saveNumBasedRuleSize,deleteNumBasedRuleSize,saveNumBasedGuidance,deleteNumBasedGuidance");
            writer.newLine();
            for (int i = 0; i < savedAndDeleteNum.size(); i++) {
                writer.write(i + "," + savedAndDeleteNum.get(i)[0] + "," + savedAndDeleteNum.get(i)[1] + "," + savedAndDeleteNum.get(i)[2] + "," + savedAndDeleteNum.get(i)[3] );
                writer.newLine();
            }
//            writer.write(state.generation-1 + "," + 0 + "," + 0);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


}
