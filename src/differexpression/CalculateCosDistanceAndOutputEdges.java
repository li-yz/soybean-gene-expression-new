package differexpression;

import utils.MyPrint;

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
        MyPrint.print("构造差异表达基因数据 完成","行数："+differGenesExpData.size() +" 列数："+differExpDataTableHeader.size());


        CalculateCosDistanceAndOutputEdges obj = new CalculateCosDistanceAndOutputEdges();
        obj.outPutDifferExpData(differGenesExpData,differExpDataTableHeader);
        MyPrint.print("输出差异表达基因数据到文件，完成！！！");
        obj.calculateCosDistanceAndConstructNetwork(differGenesExpData);

    }



    /**
     * 计算余弦相似度 矩阵
     */
    public void calculateCosDistanceAndConstructNetwork(Map<String,List<Double>> differGenesExpData){
        //先把 Map 结构的数据转换成 二维数组
        int columNum = 0;
        int rowNum = differGenesExpData.size();
        for(String key: differGenesExpData.keySet()){
            columNum = differGenesExpData.get(key).size();
            break;
        }
        String [] rowGeneIds = new String[rowNum];
        double[][] differExpArrayData = new double[rowNum][columNum];
        this.convertDifferExpMapDataToArray(differGenesExpData,differExpArrayData,rowGeneIds);

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
        //输出余弦相似度矩阵，为了看一下数据
        this.outPutCosDistanceMatrix(cosDistance);

        this.findAndOutPutEdges(cosDistance,rowGeneIds);

    }

    /**
     * 把Map<String,List<Double>>格式的差异基因表达数据转换成 数组形式，便于计算两两基因之间的相似度
     * @param differGenesExpData 原始的差异表达基因数据
     * @param arrayData 引用传递，达到返回值的目的，二维的差异表达基因数据
     * @param rowGenes  与二维数组的行相对应，每一行的基因ID
     */
    public void convertDifferExpMapDataToArray(Map<String,List<Double>> differGenesExpData,double [][] arrayData,String[] rowGenes){


        Iterator iterator = differGenesExpData.entrySet().iterator();
        int i = 0;
        while(iterator.hasNext()){
            Map.Entry<String,List<Double>> entry = (Map.Entry<String,List<Double>>)iterator.next();
            rowGenes[i] = entry.getKey();
            int columNum = entry.getValue().size();
            for(int j=0;j < columNum;j++){
                arrayData[i][j] = entry.getValue().get(j);
            }
            i++;
        }
    }

    /**
     * 计算相似度的均值
     * @return
     */
    public double getMeanValue(double [][] cosDistance){
        double sum=0.0d;
        int count=0;
        int n=cosDistance.length;
        for(int i=0;i < n;i++){
            for(int j=i+1;j < n;j++){
                if(cosDistance[i][j] > 0){
                    sum+=cosDistance[i][j];
                    count++;
                }
            }
        }
        return sum/count;

    }

    //根据相似度的阈值，找距离大于阈值的边
    public void findAndOutPutEdges(double [][] cosDistance ,String[] rowGeneIds){
        File file1=new File("D:\\paperdata\\soybean\\differExpGenes network\\genesNetworkOfDistanceThreshold3.txt");
        File file2=new File("D:\\paperdata\\soybean\\differExpGenes network\\genesNetworkOfDistanceThreshold5.txt");
        try {
            OutputStreamWriter ow1=new OutputStreamWriter(new FileOutputStream(file1), "GBK");
            BufferedWriter br1=new BufferedWriter(ow1);
            OutputStreamWriter ow2=new OutputStreamWriter(new FileOutputStream(file2), "GBK");
            BufferedWriter br2=new BufferedWriter(ow2);
            double mean = getMeanValue(cosDistance);
            MyPrint.print("差异表达基因的相似度均值：",""+mean);
            double threshold1=mean*3;//阈值取均值的3倍
            double threshold2=mean*5;//阈值取均值的5倍
            int n = cosDistance.length;

            List<String> edgesOf3 = new ArrayList<>();
            List<String> edgesOf5 = new ArrayList<>();

            for(int i=0;i < n;i++){
                for(int j=i+1;j < n;j++){
                    if(cosDistance[i][j] > 1.5*mean){
                        //如果基因i、j之间的距离大于mean，在i、j之间创建一条边
                        String edge=rowGeneIds[i]+"\t"+rowGeneIds[j]+"\n";//以对应的基因ID作为 节点
                        edgesOf3.add(edge);

                        if(cosDistance[i][j] > threshold2){
                            edgesOf5.add(edge);
                        }
                    }
                }
            }

            MyPrint.print("3倍均值的边共有：",edgesOf3.size()+"条");
            MyPrint.print("5倍均值的边共有：",edgesOf5.size()+"条");

            SimpleDateFormat df1 = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");//设置日期格式
            System.out.println("开始输出满足阈值的边："+df1.format(new Date()));// new Date()为获取当前系统时间
            Iterator iter=edgesOf3.iterator();
            while(iter.hasNext()){
                String ed=(String)iter.next();
                br1.write(ed);
            }

            br1.close();
            ow1.close();
            MyPrint.print("阈值为3倍的输出完毕："+df1.format(new Date()));// new Date()为获取当前系统时间
            Iterator iter2=edgesOf5.iterator();
            while(iter2.hasNext()){
                String ed=(String)iter2.next();
                br2.write(ed);
            }

            br2.close();
            ow2.close();
            MyPrint.print("阈值为5倍的输出完毕："+df1.format(new Date()));// new Date()为获取当前系统时间
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
