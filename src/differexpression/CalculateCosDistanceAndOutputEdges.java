package differexpression;

import utils.MyPrint;
import utils.MySerialization;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * 利用GetDifferExpGenesDataMatrix中构造出来的 差异表达基因数据
 * 计算两两基因之间的余弦相似度，满足阈值范围的 两个基因之间就建立边关系
 *
 * Created by Liyanzhen on 2016/12/9.
 */
public class CalculateCosDistanceAndOutputEdges {
    public static void main(String []args){
        GetDifferExpGenesDataMatrix getDifferExpGenesDataMatrix = new GetDifferExpGenesDataMatrix();

        //先把找到的差异表达基因集合读取进来
        Set<String> differExpGenes = new HashSet<>();
        getDifferExpGenesDataMatrix.readDifferGenes(differExpGenes);

        MyPrint.print("读取到差异表达基因ID集合","共："+differExpGenes.size());

        //开始构造 差异表达基因数据 矩阵。行表示一个基因，列是每一个GSE实验条件的均值(!!!暂时先取所有的replicates作为列！！！！)。
        //读取所有的 归一化之后的数据矩阵

        Map<String,List<Double>> differGenesExpData = new HashMap<>();//用来保存  “基因Id”---表达数据 数据表
        List<String> differExpDataTableHeader = new ArrayList<>();//用来保存数据表的 表头信息，即每一列 对应的是哪个实验条件，方便后续解释实验结果
        getDifferExpGenesDataMatrix.startReadNormalizationData(differExpGenes,differGenesExpData,differExpDataTableHeader);


        CalculateCosDistanceAndOutputEdges obj = new CalculateCosDistanceAndOutputEdges();
        obj.outPutDifferExpData(differGenesExpData,differExpDataTableHeader);
        MyPrint.print("输出差异表达基因数据到文件，完成！！！");

        int similarityType = 1;//similarityType代表选择哪种相似度度量 ：similarityType=0 时选择余弦相似度；similarityType=1时选择皮尔森相关系数
        String [] rowGeneIds = new String[differGenesExpData.size()];
        double[][]cosDistance = obj.calculateSimilarityAndConstructNetwork(differGenesExpData,rowGeneIds,similarityType);

        //释放引用，让GC能回收此内存区域。
        differGenesExpData = null;

        int isWeightedOrNot = 1;//0表示构造无权图 ，1表示加权图
        obj.findAndOutPutEdges(cosDistance,rowGeneIds,isWeightedOrNot);
    }


    /**
     * 计算相似度 矩阵
     * @param differGenesExpData 传进来的是未取多个replications均值的数据。
     * @param rowGeneIds 按序号保存对应的基因名
     * @param similarityType similarityType代表选择哪种相似度度量 ：similarityType=0 时选择余弦相似度；similarityType=1时选择皮尔森相关系数
     *
     * 为了计算的准确性，同时降低维度，应该先做convert，取每个小实验条件的均值作为一个列，  同时由于不同的GSE系列中 相同的replications分布规律不同，convert方法应该针对每一个实验条件特写
     */
    public double[][] calculateSimilarityAndConstructNetwork(Map<String,List<Double>> differGenesExpData, String [] rowGeneIds,int similarityType){
        //先把 Map 结构的数据转换成 二维数组


        double [][] differExpArrayData =this.convertDifferExpMapDataToArray(differGenesExpData,rowGeneIds);
        MyPrint.print("构造的差异表达基因矩阵，基因个数 = "+differExpArrayData.length +" 列数 = "+differExpArrayData[0].length);

        //将构造出的基因表达数据矩阵 differExpArrayData 和 基因ID--序号关系rowGeneIds序列化保存。方便再利用
        MySerialization.serializeObject(differExpArrayData,"D:\\paperdata\\soybean\\differExpressionGenes\\differExpArrayData.obj");
        MySerialization.serializeObject(rowGeneIds,"D:\\paperdata\\soybean\\differExpressionGenes\\rowGeneIds.obj");

        int geneNum = differExpArrayData.length;
        double[][] similarityMatrix = new double[geneNum][geneNum];
        if(similarityType == 0){
            similarityMatrix = getCosSimilarity(differExpArrayData);
        }else if(similarityType == 1){
            similarityMatrix = getPearsonSimilarity(differExpArrayData);
        }


        //输出余弦相似度矩阵，为了看一下数据
//        this.outPutCosDistanceMatrix(similarityMatrix);

        return similarityMatrix;

    }

    public double[][] getCosSimilarity(double[][] differExpArrayData) {
        int rowNum = differExpArrayData.length;
        int columNum = differExpArrayData[0].length;

        double [][] cosDistance = new double[rowNum][rowNum];
        //计算两两基因之间的相似度
        for(int i=0;i < rowNum;i++){
            for(int j=i+1;j < rowNum;j++){
                double mm = 0.0d;
                double modei = 0.0d;
                double modej = 0.0d;
                for(int k=0;k < columNum;k++){
                    mm += differExpArrayData[i][k]*differExpArrayData[j][k];
                    modei +=differExpArrayData[i][k]*differExpArrayData[i][k];
                    modej +=differExpArrayData[j][k]*differExpArrayData[j][k];
                }
                cosDistance[i][j] = mm/(Math.sqrt(modei) * Math.sqrt(modej));
            }

        }
        return cosDistance;
    }

    public double[][] getPearsonSimilarity(double[][] differExpArrayData) {
        int rowNum = differExpArrayData.length;
        int columNum = differExpArrayData[0].length;

        double [][] cosDistance = new double[rowNum][rowNum];
        //计算两两基因之间的相似度
        for(int i=0;i < rowNum;i++){
            for(int j=i+1;j < rowNum;j++){
                double xySum = 0.0d;
                double xSum=0.0d;
                double ySum=0.0d;
                double x2Sum=0.0d;
                double y2Sum=0.0d;
                for(int k=0;k < columNum;k++){
                    xySum+=differExpArrayData[i][k]*differExpArrayData[j][k];
                    xSum+=differExpArrayData[i][k];
                    ySum+=differExpArrayData[j][k];
                    x2Sum+=differExpArrayData[i][k]*differExpArrayData[i][k];
                    y2Sum+=differExpArrayData[j][k]*differExpArrayData[j][k];
                }
                double fenzi = xySum-(xSum*ySum)/columNum;
                double fenmu = Math.sqrt((x2Sum - (xSum*xSum)/columNum)*(y2Sum - (ySum*ySum)/columNum));
                cosDistance[i][j] = fenzi/fenmu;
            }

        }
        return cosDistance;
    }


    /**
     * 把Map<String,List<Double>>格式的差异基因表达数据转换成 数组形式，数组中的列是取均值后的值，便于计算两两基因之间的相似度
     * @param differGenesExpData 原始的差异表达基因数据
     * @param rowGenes  与二维数组的行相对应，每一行的基因ID,保存这个中间结果是为了方便解释 结论用
     */
    public double[][] convertDifferExpMapDataToArray(Map<String,List<Double>> differGenesExpData,String[] rowGenes){

        int rowNum = differGenesExpData.size();//即差异表达的基因数
        int columNum = 2+9+18+24+40+2;//即共有多少个小实验条件，也即取均值后的列数
        double [][]arrayData = new double[rowNum][columNum];

        Iterator iterator = differGenesExpData.entrySet().iterator();
        int i = 0;
        while(iterator.hasNext()){
            Map.Entry<String,List<Double>> entry = (Map.Entry<String,List<Double>>)iterator.next();
            rowGenes[i] = entry.getKey();
            int k = entry.getValue().size();//原始芯片的总个数，注意有2列数据
            List<Double> list = entry.getValue();
            //处理一行数据 计算均值
            int index=0;
            int j=0;
            double sum = 0.0d;
            //计算GSE7108系列
            sum =list.get(0)+list.get(1)+list.get(4);
            arrayData[i][j++] = sum/3;
            sum =list.get(2)+list.get(3)+list.get(5);
            arrayData[i][j++] = sum/3;
            index = 6;

            //计算GSE8432系列的 均值，index属于[6,32]
            while (index <= 14){
                if(index ==11 ||index==12){
                    sum = list.get(index+9) + list.get(index+18);
                    arrayData[i][j++]=sum/2;
                }else {
                    sum = list.get(index) + list.get(index+9) + list.get(index+18);
                    arrayData[i][j++]=sum/3;
                }
                index++;
            }

            //计算 GSE29740、29741系列, index属于[33,158]
            index = 33;
            while(index <= 158){
                sum = list.get(index++) + list.get(index++) + list.get(index++);
                arrayData[i][j++] = sum/3;
            }
            //计算 GSE33410系列, index属于[159,278]
            while(index <= 278){
                sum = list.get(index++) + list.get(index++) + list.get(index++);
                arrayData[i][j++] = sum/3;
            }
            //计算 GSE41724系列, index属于[279,283]
            sum = list.get(index++) + list.get(index++);
            arrayData[i][j++] = sum/2;
            sum = list.get(index++) + list.get(index++) +list.get(index);
            arrayData[i][j++] = sum/3;

            i++;
        }
        return arrayData;
    }

    /**
     * 计算相似度的均值
     * @return
     */
    public double[] getMeanValue(double [][] similarityMatrix){
        double sumPositive=0.0d;
        int countPositive=0;
        double sumNegative = 0.0d;
        int countNegative=0;
        int n=similarityMatrix.length;
        double[] means = new double[2];
        for(int i=0;i < n;i++){
            for(int j=i+1;j < n;j++){
                if(similarityMatrix[i][j] > 0){
                    sumPositive+=similarityMatrix[i][j];
                    countPositive++;
                }
                if(similarityMatrix[i][j] < 0){
                    sumNegative+=similarityMatrix[i][j];
                    countNegative++;
                }
            }
        }
        means[0] = sumPositive/countPositive;
        means[1] = sumNegative/countNegative;
        return means;

    }


    /**
     *根据相似度值来筛选满足阈值的边，并创建想要的网络图（加权图 or 无权图）
     * @param similarityMatrix 基因相似度矩阵
     * @param rowGeneIds 按序号排列的基因名
     * @param isWeightedOrNot 选择要创建的是加权图还是无权图，0表示无权图；1表示加权图
     */
    public void findAndOutPutEdges(double [][] similarityMatrix ,String[] rowGeneIds,int isWeightedOrNot){
        File file1=new File("D:\\paperdata\\soybean\\differExpGenes network\\genesNetworkOfSimilarityThreshold1.txt");
        File file2=new File("D:\\paperdata\\soybean\\differExpGenes network\\genesNetworkOfSimilarityThreshold2.txt");
        try {
            OutputStreamWriter ow1=new OutputStreamWriter(new FileOutputStream(file1), "GBK");
            BufferedWriter br1=new BufferedWriter(ow1);
            OutputStreamWriter ow2=new OutputStreamWriter(new FileOutputStream(file2), "GBK");
            BufferedWriter br2=new BufferedWriter(ow2);
            double[] means = getMeanValue(similarityMatrix);
            MyPrint.print("差异表达基因的相似度均值：正相关的均值="+means[0]+", 负相关的均值 = "+means[1]);
            double threshold1Positive=means[0]*2;//阈值取相似度的倍数
            double threshold2Positive=means[0]*2.78;

            double threshold1Negative=means[1]*2;
            double threshold2Negative=means[1]*3.7;

            int n = similarityMatrix.length;

            List<String> edgesOfThreshold1 = new ArrayList<>();
            List<String> edgesOfThreshold2 = new ArrayList<>();
            int countTh1Positive = 0;
            int countTh1Negative = 0;
            int countTh2Positive = 0;
            int countTh2Negative = 0;

            for(int i=0;i < n;i++){
                for(int j=i+1;j < n;j++){
                    if(similarityMatrix[i][j] > threshold1Positive){
                        //如果基因i、j之间的距离大于mean，在i、j之间创建一条边
                        String edge="";
                        if(isWeightedOrNot == 0) {
                            edge = rowGeneIds[i] + "\t" + rowGeneIds[j] + "\t"+ 1 + "\n";//以对应的基因ID作为 节点
                        }else if(isWeightedOrNot == 1) {
                            edge = rowGeneIds[i] + "\t" + rowGeneIds[j] + "\t"+similarityMatrix[i][j] + "\n";//以对应的基因ID作为 节点

                        }
                        edgesOfThreshold1.add(edge);
                        countTh1Positive++;

                        if(similarityMatrix[i][j] > threshold2Positive){
                            edgesOfThreshold2.add(edge);
                            countTh2Positive++;
                        }
                    }
                    if(similarityMatrix[i][j] < threshold1Negative){
                        String edge="";
                        if(isWeightedOrNot == 0) {
                            edge = rowGeneIds[i] + "\t" + rowGeneIds[j] + "\t"+ 1 + "\n";//以对应的基因ID作为 节点，一行：from to 值为1代表无权图
                        }else if(isWeightedOrNot == 1) {
                            edge = rowGeneIds[i] + "\t" + rowGeneIds[j] + "\t"+similarityMatrix[i][j] + "\n";//以对应的基因ID作为 节点 ，一行：from to 相应权值

                        }
                        edgesOfThreshold1.add(edge);
                        countTh1Negative++;
                        if(similarityMatrix[i][j] < threshold2Negative){
                            edgesOfThreshold2.add(edge);
                            countTh2Negative++;
                        }
                    }
                }
            }

            MyPrint.print("阈值 > threshold1Positive:"+threshold1Positive+"的边共有：",countTh1Positive+"条");
            MyPrint.print("阈值 < threshold1Negative:"+threshold1Negative+"的边共有：",countTh1Negative+"条");
            MyPrint.print("阈值 > threshold2Positive:"+threshold2Positive+"的边共有：",countTh2Positive+"条");
            MyPrint.print("阈值 < threshold2Negative:"+threshold2Negative+"的边共有：",countTh2Negative+"条");

            SimpleDateFormat df1 = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");//设置日期格式
            System.out.println("开始输出满足阈值的边："+df1.format(new Date()));// new Date()为获取当前系统时间
            Iterator iter=edgesOfThreshold1.iterator();
            while(iter.hasNext()){
                String ed=(String)iter.next();
                br1.write(ed);
            }

            br1.close();
            ow1.close();
            MyPrint.print("阈值"+threshold1Positive+"的边输出完毕："+df1.format(new Date()));// new Date()为获取当前系统时间
            Iterator iter2=edgesOfThreshold2.iterator();
            while(iter2.hasNext()){
                String ed=(String)iter2.next();
                br2.write(ed);
            }

            br2.close();
            ow2.close();
            MyPrint.print("阈值"+threshold2Positive+"的边输出完毕："+df1.format(new Date()));// new Date()为获取当前系统时间
        } catch (UnsupportedEncodingException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    //根据相似度的阈值，找距离大于阈值的边
    public void outPutDifferExpData(Map<String,List<Double>> differExpData, List<String> header){
        File file1=new File("D:\\paperdata\\soybean\\differExpressionGenes\\construct_differ_exp_data.txt");
        try {
            OutputStreamWriter ow1=new OutputStreamWriter(new FileOutputStream(file1), "GBK");
            BufferedWriter bw=new BufferedWriter(ow1);
            StringBuffer sb = new StringBuffer();
            for(String str: header){
                sb.append(str);
                sb.append("\t");
            }
            sb.append("\n");
            bw.write(sb.toString());
            for(Map.Entry<String,List<Double>> entry :differExpData.entrySet()){
                sb = new StringBuffer();
                sb.append(entry.getKey());
                sb.append("\t");
                for(int i=0;i <entry.getValue().size()-1;i++){
                    sb.append(entry.getValue().get(i));
                    sb.append("\t");
                }
                sb.append(entry.getValue().get(entry.getValue().size()-1));
                sb.append("\n");
                bw.write(sb.toString());
            }


            bw.close();
            ow1.close();
        } catch (UnsupportedEncodingException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }

    private void outPutCosDistanceMatrix(double[][] cosDistance){
        File file1=new File("D:\\paperdata\\soybean\\differExpressionGenes\\cos_distance_matrix_data.txt");
        try {
            OutputStreamWriter ow1 = new OutputStreamWriter(new FileOutputStream(file1), "GBK");
            BufferedWriter bw = new BufferedWriter(ow1);
            for (int i = 0; i < cosDistance.length; i++) {
                StringBuffer sb = new StringBuffer();
                for(int j=0;j < cosDistance[i].length;j++){
                    sb.append(cosDistance[i][j]);
                    sb.append("\t");
                }
                sb.append("\n");
                bw.write(sb.toString());
            }
            bw.close();
            ow1.close();
        } catch (UnsupportedEncodingException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

}
