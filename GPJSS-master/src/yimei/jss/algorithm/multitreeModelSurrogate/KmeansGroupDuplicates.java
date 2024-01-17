package yimei.jss.algorithm.multitreeModelSurrogate;

import ec.EvolutionState;

import java.util.ArrayList;

public class KmeansGroupDuplicates {
    private int k;// 分成多少簇
    private int m;// 迭代次数
    private int dataSetLength;// 数据集元素个数，即数据集的长度
    private ArrayList<double[]> dataSet;// 数据集链表
    public ArrayList<double[]> center;// 中心链表
    private ArrayList<ArrayList<double[]>> cluster; // 簇
    private ArrayList<Double> jc;// 误差平方和，k越接近dataSetLength，误差越小
//    private Random random;


    /**
     * 设置需分组的原始数据集
     *
     * @param dataSet
     */

    public void setDataSet(ArrayList<double[]> dataSet) {
        this.dataSet = dataSet;
    }

    /**
     * 获取结果分组
     *
     * @return 结果集
     */

    public ArrayList<ArrayList<double[]>> getCluster() {
        return cluster;
    }

    /**
     * 构造函数，传入需要分成的簇数量
     *
     * @param k 簇数量,若k<=0时，设置为1，若k大于数据源的长度时，置为数据源的长度
     */
    public KmeansGroupDuplicates(int k) {
        if (k <= 0) {
            k = 1;
        }
        this.k = k;
    }

    /**
     * 初始化
     */
    private void init(final EvolutionState state) {
        m = 0;
        //random = new Random();
        if (dataSet == null || dataSet.size() == 0) {
            initDataSet();
        }
        dataSetLength = dataSet.size();
        if (k > dataSetLength) {
            k = dataSetLength;
        }
        center = initCenters(state);
        cluster = initCluster();
        jc = new ArrayList<Double>();
    }

    /**
     * 如果调用者未初始化数据集，则采用内部测试数据集
     */
    private void initDataSet() {
        dataSet = new ArrayList<double[]>();
        // 其中{6,3}是一样的，所以长度为15的数据集分成14簇和15簇的误差都为0
        double[][] dataSetArray = new double[][]{{8, 2}, {3, 4}, {2, 5},
                {4, 2}, {7, 3}, {6, 2}, {4, 7}, {6, 3}, {5, 3},
                {6, 3}, {6, 9}, {1, 6}, {3, 9}, {4, 1}, {8, 6}};

        for (int i = 0; i < dataSetArray.length; i++) {
            dataSet.add(dataSetArray[i]);
        }
    }

    /**
     * 初始化中心数据链表，分成多少簇就有多少个中心点,
     *
     * @return 中心点集 ;随机从个体中选出K个作为中心点；
     */
    private ArrayList<double[]> initCenters(final EvolutionState state) {
        ArrayList<double[]> center = new ArrayList<double[]>();
        int[] randoms = new int[k];
        boolean flag;
        int temp = 0;
        //int temp = random.nextInt(dataSetLength);
        for (int i=0; i<dataSetLength; i++){
            if(dataSet.get(i)[0] == Double.MAX_VALUE){
                temp =i;
                break;}
        }
        randoms[0] = temp;
        int alreadyChosenNum = 1;
        for (int i = 1; i < k; i++) {
            flag = true;

            while (flag) {
                int judge = 1;
                //temp = random.nextInt(dataSetLength);
                temp = state.random[0].nextInt(dataSetLength);

                for (int a=0; a<alreadyChosenNum ; a++){
                if(dataSet.get(temp)[0] == dataSet.get(randoms[a])[0]) {  //避免选中重复的中心
                    judge = 0;
                    break;}
                }
                if(judge == 0)
                    continue;

                int j = 0;
                // 不清楚for循环导致j无法加1
                // for(j=0;j<i;++j)
                // {
                // if(temp==randoms[j]);
                // {
                // break;
                // }
                // }
                while (j < i) {
                    if (temp == randoms[j]) {
                        break;
                    }
                    j++;
                }
                if (j == i) {
                    flag = false;
                }
            }
            randoms[i] = temp;
            alreadyChosenNum++;
        }

        // 测试随机数生成情况
        // for(int i=0;i<k;i++)
        // {
        // System.out.println("test1:randoms["+i+"]="+randoms[i]);
        // }

        // System.out.println();
        for (int i = 0; i < k; i++) {
            center.add(dataSet.get(randoms[i]));// 生成初始化中心链表
        }
        return center;
    }

    /**
     * 初始化簇集合
     *
     * @return 一个分为k簇的空数据的簇集合
     */
    private ArrayList<ArrayList<double[]>> initCluster() {
        ArrayList<ArrayList<double[]>> cluster = new ArrayList<ArrayList<double[]>>();
        for (int i = 0; i < k; i++) {
            cluster.add(new ArrayList<double[]>());
        }

        return cluster;
    }

    /**
     * 计算两个点之间的距离
     *
     * @param element 点1
     * @param center  点2
     * @return 距离
     */
    private double distance(double[] element, double[] center) {
        double distance = 0.0;
        double y = 0.0f;
        for (int i = 0; i < element.length; i++) {
            double x = element[i] - center[i];
            y += x * x;
        }
        distance = (double) Math.sqrt(y);

        return distance;
    }

    /**
     * 获取距离集合中最小距离的位置
     *
     * @param distance 距离数组
     * @return 最小距离在距离数组中的位置
     */
    private int minDistance(double[] distance,final EvolutionState state) {
        double minDistance = distance[0];
        int minLocation = 0;
        for (int i = 1; i < distance.length; i++) {
            if (distance[i] < minDistance) {
                minDistance = distance[i];
                minLocation = i;
            } else if (distance[i] == minDistance) // 如果相等，随机返回一个位置
            {
                if (state.random[0].nextInt(10) < 5) {
                    minLocation = i;
                }
            }
        }

        return minLocation;
    }

    /**
     * 核心，将当前元素放到最小距离中心相关的簇中
     */
    private void clusterSet(final EvolutionState state) {
        double[] distance = new double[k];
        for (int i = 0; i < dataSetLength; i++) {
            if(dataSet.get(i)[0] == Double.MAX_VALUE){   //因为初始化的原因，若存在Doouble.Max 第一个簇永远都是Double.Max
                cluster.get(0).add(dataSet.get(i));
                continue;}
            for (int j = 0; j < k; j++) {
                if(center.get(j)[0] == Double.MAX_VALUE){
                    distance[j] = Double.MAX_VALUE;
                    continue;}

                distance[j] = distance(dataSet.get(i), center.get(j));
                // System.out.println("test2:"+"dataSet["+i+"],center["+j+"],distance="+distance[j]);

            }
            int minLocation = minDistance(distance,state);
            // System.out.println("test3:"+"dataSet["+i+"],minLocation="+minLocation);
            // System.out.println();

            cluster.get(minLocation).add(dataSet.get(i));// 核心，将当前元素放到最小距离中心相关的簇中


        }
        for (int i=0; i<k ; i++){
            if(cluster.get(i).size() == 0){
                cluster.remove(i);
                center.remove(i);
                k--;
            }

        }

    }

    /**
     * 求两点误差平方的方法
     *
     * @param element 点1
     * @param center  点2
     * @return 误差平方
     */
    private double errorSquare(double[] element, double[] center) {
        double errSquare = 0.0f;
        for (int i = 0; i < element.length; i++) {
            double x = element[i] - center[i];
            errSquare += Math.abs(x);
        }

        return errSquare;
    }

    /**
     * 计算误差平方和准则函数方法；jc是每个簇中每个个体离中心点距离的平方和的和
     */
    private void countRule() {
        double jcF = 0;
        for (int i = 0; i < cluster.size(); i++) {
/*            if(cluster.get(i).get(0)[0] == Double.MAX_VALUE)
                continue;*/
            for (int j = 0; j < cluster.get(i).size(); j++) {
                jcF += errorSquare(cluster.get(i).get(j), center.get(i));

            }
        }
        jc.add(jcF);
    }

    /**
     * 设置新的簇中心方法
     */
    public void setNewCenter() {
        for (int i = 0; i < k; i++) {
            if(center.get(i)[0] == Double.MAX_VALUE)
                continue;
            int n = cluster.get(i).size();
            if (n != 0) {
                double[] newCenter = {0};//如何表示一个一维40个0元素；
                for (int h = 0; h < 1; h++) {
                    for (int j = 0; j < n; j++) {
                        newCenter[h] += cluster.get(i).get(j)[h];
                    }
                    // 设置一个平均值
                    newCenter[h] = newCenter[h] / n;
                }
                center.set(i, newCenter);
            }
        }
    }

    /*    *//**
     * 打印数据，测试用
     *
     * @param dataArray
     *            数据集
     * @param dataArrayName
     *            数据集名称
     *//*
    public void printDataArray(ArrayList<float[]> dataArray,
                               String dataArrayName) {
        for (int i = 0; i < dataArray.size(); i++) {
            System.out.println("print:" + dataArrayName + "[" + i + "]={"
                    + dataArray.get(i)[0] + "," + dataArray.get(i)[1] + "}");
        }
        System.out.println("===================================");
    }*/

    /**
     * Kmeans算法核心过程方法
     */
    private void kmeans(final EvolutionState state) {
        init(state);
        // printDataArray(dataSet,"initDataSet");
        /*printDataArray(center,"initCenter");*/

        // 循环分组，直到误差不变为止
        while (true) {
            clusterSet(state);
            // for(int i=0;i<cluster.size();i++)
            // {
            // printDataArray(cluster.get(i),"cluster["+i+"]");
            // }

            countRule();

            // System.out.println("count:"+"jc["+m+"]="+jc.get(m));

            // System.out.println();
            // 误差不变了，分组完成
            if (m != 0) {
                if (jc.get(m) - jc.get(m - 1) == 0) {
                    break;
                }
            }

            if (m >=10)
                break;

            setNewCenter();
            /* printDataArray(center,"newCenter");*/
            m++;
            cluster.clear();
            cluster = initCluster();
        }

         //System.out.println("note:the times of repeat:m="+m);//输出迭代次数
    }

    /**
     * 执行算法
     */
    public void execute(final EvolutionState state) {
     /*   long startTime = System.currentTimeMillis();
        System.out.println("kmeans begins");*/
        kmeans(state);
  /*      long endTime = System.currentTimeMillis();
        System.out.println("kmeans running time=" + (endTime - startTime)
                + "ms");
        System.out.println("kmeans ends");
        *//*System.out.println(center[1]);*//*
        System.out.println();*/
        // center  cluster
        //center.sort();

       /* ArrayList<double[]> centerCopy = new ArrayList<>();
        ArrayList<ArrayList<double[]>> clusterCopy = new ArrayList<>();

        for (int a=0; a<center.size() ; a++){
            centerCopy.add(center.get(a));
            clusterCopy.add(cluster.get(a));
        }

        center.clear();
        cluster.clear();

        for (int i = 0; i < centerCopy.size() - 1; ++i) {
            for (int j = 0; j < centerCopy.size() - 1 -i; ++j) {
                if(centerCopy.get(j)[0] > centerCopy.get(j + 1)[0]) {
                    double temp0 = centerCopy.get(j)[0];
                    centerCopy.get(j)[0] = centerCopy.get(j+1)[0];
                    centerCopy.get(j+1)[0] = temp0;

                    ArrayList<double[]> temp2 = (ArrayList<double[]>) clusterCopy.get(j).clone();
                    clusterCopy.get(j).clear();
                    clusterCopy.get(j).addAll(clusterCopy.get(j+1));
                    clusterCopy.get(j+1).clear();
                    clusterCopy.get(j+1).addAll(temp2);
                }
            }
        }

        center = (ArrayList<double[]>) centerCopy.clone();
        cluster = (ArrayList<ArrayList<double[]>>) clusterCopy.clone();*/

    }
}
