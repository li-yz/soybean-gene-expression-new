package differexpression;

import utils.MyPrint;
import utils.MySerialization;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * ����GetDifferExpGenesDataMatrix�й�������� �������������
 * ������������֮����������ƶȣ�������ֵ��Χ�� ��������֮��ͽ����߹�ϵ
 *
 * Created by Liyanzhen on 2016/12/9.
 */
public class CalculateCosDistanceAndOutputEdges {
    public static void main(String []args){
        GetDifferExpGenesDataMatrix getDifferExpGenesDataMatrix = new GetDifferExpGenesDataMatrix();

        //�Ȱ��ҵ��Ĳ�������򼯺϶�ȡ����
        Set<String> differExpGenes = new HashSet<>();
        getDifferExpGenesDataMatrix.readDifferGenes(differExpGenes);

        MyPrint.print("��ȡ�����������ID����","����"+differExpGenes.size());

        //��ʼ���� ������������� �����б�ʾһ����������ÿһ��GSEʵ�������ľ�ֵ(!!!��ʱ��ȡ���е�replicates��Ϊ�У�������)��
        //��ȡ���е� ��һ��֮������ݾ���

        Map<String,List<Double>> differGenesExpData = new HashMap<>();//��������  ������Id��---������� ���ݱ�
        List<String> differExpDataTableHeader = new ArrayList<>();//�����������ݱ�� ��ͷ��Ϣ����ÿһ�� ��Ӧ�����ĸ�ʵ�������������������ʵ����
        getDifferExpGenesDataMatrix.startReadNormalizationData(differExpGenes,differGenesExpData,differExpDataTableHeader);


        CalculateCosDistanceAndOutputEdges obj = new CalculateCosDistanceAndOutputEdges();
        obj.outPutDifferExpData(differGenesExpData,differExpDataTableHeader);
        MyPrint.print("���������������ݵ��ļ�����ɣ�����");

        int similarityType = 1;//similarityType����ѡ���������ƶȶ��� ��similarityType=0 ʱѡ���������ƶȣ�similarityType=1ʱѡ��Ƥ��ɭ���ϵ��
        String [] rowGeneIds = new String[differGenesExpData.size()];
        double[][]cosDistance = obj.calculateSimilarityAndConstructNetwork(differGenesExpData,rowGeneIds,similarityType);

        //�ͷ����ã���GC�ܻ��մ��ڴ�����
        differGenesExpData = null;

        int isWeightedOrNot = 1;//0��ʾ������Ȩͼ ��1��ʾ��Ȩͼ
        obj.findAndOutPutEdges(cosDistance,rowGeneIds,isWeightedOrNot);
    }


    /**
     * �������ƶ� ����
     * @param differGenesExpData ����������δȡ���replications��ֵ�����ݡ�
     * @param rowGeneIds ����ű����Ӧ�Ļ�����
     * @param similarityType similarityType����ѡ���������ƶȶ��� ��similarityType=0 ʱѡ���������ƶȣ�similarityType=1ʱѡ��Ƥ��ɭ���ϵ��
     *
     * Ϊ�˼����׼ȷ�ԣ�ͬʱ����ά�ȣ�Ӧ������convert��ȡÿ��Сʵ�������ľ�ֵ��Ϊһ���У�  ͬʱ���ڲ�ͬ��GSEϵ���� ��ͬ��replications�ֲ����ɲ�ͬ��convert����Ӧ�����ÿһ��ʵ��������д
     */
    public double[][] calculateSimilarityAndConstructNetwork(Map<String,List<Double>> differGenesExpData, String [] rowGeneIds,int similarityType){
        //�Ȱ� Map �ṹ������ת���� ��ά����


        double [][] differExpArrayData =this.convertDifferExpMapDataToArray(differGenesExpData,rowGeneIds);
        MyPrint.print("����Ĳ����������󣬻������ = "+differExpArrayData.length +" ���� = "+differExpArrayData[0].length);

        //��������Ļ��������ݾ��� differExpArrayData �� ����ID--��Ź�ϵrowGeneIds���л����档����������
        MySerialization.serializeObject(differExpArrayData,"D:\\paperdata\\soybean\\differExpressionGenes\\differExpArrayData.obj");
        MySerialization.serializeObject(rowGeneIds,"D:\\paperdata\\soybean\\differExpressionGenes\\rowGeneIds.obj");

        int geneNum = differExpArrayData.length;
        double[][] similarityMatrix = new double[geneNum][geneNum];
        if(similarityType == 0){
            similarityMatrix = getCosSimilarity(differExpArrayData);
        }else if(similarityType == 1){
            similarityMatrix = getPearsonSimilarity(differExpArrayData);
        }


        //����������ƶȾ���Ϊ�˿�һ������
//        this.outPutCosDistanceMatrix(similarityMatrix);

        return similarityMatrix;

    }

    public double[][] getCosSimilarity(double[][] differExpArrayData) {
        int rowNum = differExpArrayData.length;
        int columNum = differExpArrayData[0].length;

        double [][] cosDistance = new double[rowNum][rowNum];
        //������������֮������ƶ�
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
        //������������֮������ƶ�
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
     * ��Map<String,List<Double>>��ʽ�Ĳ������������ת���� ������ʽ�������е�����ȡ��ֵ���ֵ�����ڼ�����������֮������ƶ�
     * @param differGenesExpData ԭʼ�Ĳ������������
     * @param rowGenes  ���ά����������Ӧ��ÿһ�еĻ���ID,��������м�����Ϊ�˷������ ������
     */
    public double[][] convertDifferExpMapDataToArray(Map<String,List<Double>> differGenesExpData,String[] rowGenes){

        int rowNum = differGenesExpData.size();//��������Ļ�����
        int columNum = 2+9+18+24+40+2;//�����ж��ٸ�Сʵ��������Ҳ��ȡ��ֵ�������
        double [][]arrayData = new double[rowNum][columNum];

        Iterator iterator = differGenesExpData.entrySet().iterator();
        int i = 0;
        while(iterator.hasNext()){
            Map.Entry<String,List<Double>> entry = (Map.Entry<String,List<Double>>)iterator.next();
            rowGenes[i] = entry.getKey();
            int k = entry.getValue().size();//ԭʼоƬ���ܸ�����ע����2������
            List<Double> list = entry.getValue();
            //����һ������ �����ֵ
            int index=0;
            int j=0;
            double sum = 0.0d;
            //����GSE7108ϵ��
            sum =list.get(0)+list.get(1)+list.get(4);
            arrayData[i][j++] = sum/3;
            sum =list.get(2)+list.get(3)+list.get(5);
            arrayData[i][j++] = sum/3;
            index = 6;

            //����GSE8432ϵ�е� ��ֵ��index����[6,32]
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

            //���� GSE29740��29741ϵ��, index����[33,158]
            index = 33;
            while(index <= 158){
                sum = list.get(index++) + list.get(index++) + list.get(index++);
                arrayData[i][j++] = sum/3;
            }
            //���� GSE33410ϵ��, index����[159,278]
            while(index <= 278){
                sum = list.get(index++) + list.get(index++) + list.get(index++);
                arrayData[i][j++] = sum/3;
            }
            //���� GSE41724ϵ��, index����[279,283]
            sum = list.get(index++) + list.get(index++);
            arrayData[i][j++] = sum/2;
            sum = list.get(index++) + list.get(index++) +list.get(index);
            arrayData[i][j++] = sum/3;

            i++;
        }
        return arrayData;
    }

    /**
     * �������ƶȵľ�ֵ
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
     *�������ƶ�ֵ��ɸѡ������ֵ�ıߣ���������Ҫ������ͼ����Ȩͼ or ��Ȩͼ��
     * @param similarityMatrix �������ƶȾ���
     * @param rowGeneIds ��������еĻ�����
     * @param isWeightedOrNot ѡ��Ҫ�������Ǽ�Ȩͼ������Ȩͼ��0��ʾ��Ȩͼ��1��ʾ��Ȩͼ
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
            MyPrint.print("�������������ƶȾ�ֵ������صľ�ֵ="+means[0]+", ����صľ�ֵ = "+means[1]);
            double threshold1Positive=means[0]*2;//��ֵȡ���ƶȵı���
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
                        //�������i��j֮��ľ������mean����i��j֮�䴴��һ����
                        String edge="";
                        if(isWeightedOrNot == 0) {
                            edge = rowGeneIds[i] + "\t" + rowGeneIds[j] + "\t"+ 1 + "\n";//�Զ�Ӧ�Ļ���ID��Ϊ �ڵ�
                        }else if(isWeightedOrNot == 1) {
                            edge = rowGeneIds[i] + "\t" + rowGeneIds[j] + "\t"+similarityMatrix[i][j] + "\n";//�Զ�Ӧ�Ļ���ID��Ϊ �ڵ�

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
                            edge = rowGeneIds[i] + "\t" + rowGeneIds[j] + "\t"+ 1 + "\n";//�Զ�Ӧ�Ļ���ID��Ϊ �ڵ㣬һ�У�from to ֵΪ1������Ȩͼ
                        }else if(isWeightedOrNot == 1) {
                            edge = rowGeneIds[i] + "\t" + rowGeneIds[j] + "\t"+similarityMatrix[i][j] + "\n";//�Զ�Ӧ�Ļ���ID��Ϊ �ڵ� ��һ�У�from to ��ӦȨֵ

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

            MyPrint.print("��ֵ > threshold1Positive:"+threshold1Positive+"�ı߹��У�",countTh1Positive+"��");
            MyPrint.print("��ֵ < threshold1Negative:"+threshold1Negative+"�ı߹��У�",countTh1Negative+"��");
            MyPrint.print("��ֵ > threshold2Positive:"+threshold2Positive+"�ı߹��У�",countTh2Positive+"��");
            MyPrint.print("��ֵ < threshold2Negative:"+threshold2Negative+"�ı߹��У�",countTh2Negative+"��");

            SimpleDateFormat df1 = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");//�������ڸ�ʽ
            System.out.println("��ʼ���������ֵ�ıߣ�"+df1.format(new Date()));// new Date()Ϊ��ȡ��ǰϵͳʱ��
            Iterator iter=edgesOfThreshold1.iterator();
            while(iter.hasNext()){
                String ed=(String)iter.next();
                br1.write(ed);
            }

            br1.close();
            ow1.close();
            MyPrint.print("��ֵ"+threshold1Positive+"�ı������ϣ�"+df1.format(new Date()));// new Date()Ϊ��ȡ��ǰϵͳʱ��
            Iterator iter2=edgesOfThreshold2.iterator();
            while(iter2.hasNext()){
                String ed=(String)iter2.next();
                br2.write(ed);
            }

            br2.close();
            ow2.close();
            MyPrint.print("��ֵ"+threshold2Positive+"�ı������ϣ�"+df1.format(new Date()));// new Date()Ϊ��ȡ��ǰϵͳʱ��
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

    //�������ƶȵ���ֵ���Ҿ��������ֵ�ı�
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
