package differexpression;

import serialprocess.OverlapPartition;
import serialprocess.Partition;
import utils.MyPrint;
import utils.MySerialization;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;

/**
 * Created by Liyanzhen on 2017/3/23.
 */
public class HeatmapData {
    public static void main(String[] args){
        //������һ�㣺ɸѡ���Ĳ��������7971���������ı�ɸѡ���������ķ�������allFoldChange�е����ݶ��ǲ���ġ�����������ѡ���������ͼG1�Ľ������ G2�Ľ����allFoldChange���ǲ����
        Map<String,Map<String,Double>> allFoldChange = (Map<String,Map<String,Double>>) MySerialization.antiSerializeObject("D:\\soybean project\\soybean-differ-expression\\allFoldChange\\allFoldChange.obj");
        MyPrint.print("��:"+allFoldChange.size()+" ��fold-change����");

        GetDifferExpGenesDataMatrix getDifferExpGenesDataMatrix = new GetDifferExpGenesDataMatrix();

        //�Ȱ��ҵ��Ĳ�������򼯺϶�ȡ����
        Set<String> differExpGenes = new HashSet<>();
        getDifferExpGenesDataMatrix.readDifferGenes(differExpGenes);

        //�ȷ�������ͼG1�Ļ��ֽ��������---����������ϵ
//        Partition bestNonOverlapPartition = (Partition) MySerialization.antiSerializeObject("D:\\paperdata\\soybean\\community detection\\��ʷ������\\2017.2.26\\bestNonOverlapPartition.obj");

        //����ͼG2���������ֽ��
        Partition bestNonOverlapPartition = (Partition) MySerialization.antiSerializeObject("D:\\paperdata\\soybean\\community detection\\��ʷ������\\2017.3.9����ͼG2\\bestNonOverlapPartition.obj");
        Set<String> genesInGraph = new HashSet<>();
        genesInGraph.addAll(bestNonOverlapPartition.getNodeCommunityMap().keySet());

        prepareDataForHeatmap(allFoldChange,genesInGraph,bestNonOverlapPartition);
        prepareDataForHeatmap(allFoldChange,differExpGenes);
    }

    /**
     *������Ҫ��fold-change����
     * @param allFoldChange ȫ������������fold-changeֵ
     * @param genes ��Ҫ����Ļ����б�����������ͼG2����G1�еĻ���
     * @param bestNonOverlapPartition
     */
    public static void prepareDataForHeatmap(Map<String,Map<String,Double>> allFoldChange ,Set<String>genes ,Partition bestNonOverlapPartition){
        try{
            FileWriter writer = new FileWriter("D:\\soybean project\\soybean-differ-expression\\allFoldChange\\allFoldChange.csv");
            BufferedWriter bw = new BufferedWriter(writer);

            FileWriter wr = new FileWriter("D:\\soybean project\\soybean-differ-expression\\allFoldChange\\geneMapComIndex.csv");
            BufferedWriter bw2 = new BufferedWriter(wr);

            //����geneMapComIndex.csv�ı�ͷ
            StringBuffer sb = new StringBuffer();
            sb.append("geneId");
            sb.append("\t");
            sb.append("communityIndex");
            sb.append("\t");
            sb.append("communityName");
            bw2.write(sb.toString());
            bw2.newLine();
            //׼��һ�ݡ�������---������š���ӳ���ϵ���ݣ���������size < 10������ͳһ�����0���������Ĵ�������1��ʼ�������
            Map<String,Integer> comNameMapIndex = new HashMap<>();
            int index=1;
            for(Map.Entry entry:bestNonOverlapPartition.getCommunities().entrySet()){
                String comName = (String) entry.getKey();
                List<String> community= (List<String>)entry.getValue();
                int size = community.size();

                if(size < 10){
                    comNameMapIndex.put(comName,0);//size < 10������ͳһ�������0������
                }else {
                    comNameMapIndex.put(comName,index++);
                }
            }
            MyPrint.print("size >= 10������ ���������index="+(index-1));


            //����allFoldChange.csv�ı�ͷ
            List<String> allExperments = new ArrayList<>();
            allExperments.addAll(allFoldChange.keySet());
            sb = new StringBuffer();
            sb.append("geneId");
            for(String expermentType :allExperments){
                sb.append("\t");
                sb.append(expermentType);
            }
            bw.write(sb.toString());//��ӡ��ͷ
            bw.newLine();
            int test = 0;
            for(String gene :genes){//�������в��������
                sb = new StringBuffer();
                sb.append(gene);

                for(String expName: allExperments){//��������ʵ���������е�GSE�µ������ҷֳ���2��ʵ��������fold-changeֵ
                    sb.append("\t");
                    sb.append(log2(allFoldChange.get(expName).get(gene)));
                }
                bw.write(sb.toString());
                bw.newLine();

                //���濪ʼ����geneMapComIndex.csv������
                if(bestNonOverlapPartition.getNodeCommunityMap().containsKey(gene)){
                    sb = new StringBuffer();
                    sb.append(gene);
                    sb.append("\t");
                    sb.append(comNameMapIndex.get(bestNonOverlapPartition.getNodeCommunityMap().get(gene)));//gene�����������
                    sb.append("\t");
                    sb.append(bestNonOverlapPartition.getNodeCommunityMap().get(gene));//gene����������
                    bw2.write(sb.toString());
                    bw2.newLine();
                    test++;
                }else {
                    sb = new StringBuffer();
                    sb.append(gene);
                    sb.append("\t");
                    sb.append(0);
                    sb.append("\t");
                    sb.append("none");
                    bw2.write(sb.toString());
                    bw2.newLine();
                    test++;
                }

            }

            bw.close();
            writer.close();

            bw2.close();
            wr.close();
        }catch (Exception e){
            MyPrint.print("���heatmap�����fold-changeֵ�쳣");
            e.printStackTrace();
        }
    }

    /**
     *�˷��������ȫ������������fold-change����
     * @param allFoldChange
     */
    public static void prepareDataForHeatmap(Map<String,Map<String,Double>> allFoldChange ,Set<String>differExpGenes ){
        try{
            FileWriter writer = new FileWriter("D:\\soybean project\\soybean-differ-expression\\allFoldChange\\allDifferGenesFoldChange.csv");
            BufferedWriter bw = new BufferedWriter(writer);

            //����allFoldChange.csv�ı�ͷ
            List<String> allExperments = new ArrayList<>();
            allExperments.addAll(allFoldChange.keySet());
            StringBuffer sb = new StringBuffer();
            sb.append("geneId");
            for(String expermentType :allExperments){
                sb.append("\t");
                sb.append(expermentType);
            }
            bw.write(sb.toString());//��ӡ��ͷ
            bw.newLine();
            int test = 0;
            for(String gene :differExpGenes){//�������в��������
                sb = new StringBuffer();
                sb.append(gene);

                for(String expName: allExperments){//��������ʵ���������е�GSE�µ������ҷֳ���2��ʵ��������fold-changeֵ
                    sb.append("\t");
                    sb.append(log2(allFoldChange.get(expName).get(gene)));
                }
                bw.write(sb.toString());
                bw.newLine();

            }

            bw.close();
            writer.close();
        }catch (Exception e){
            MyPrint.print("���heatmap�����fold-changeֵ�쳣");
            e.printStackTrace();
        }
    }

    private static double log2(double data){
        return Math.log(data)/Math.log(2);
    }
}
