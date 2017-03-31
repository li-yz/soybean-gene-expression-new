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
        //别忘了一点：筛选出的差异表达基因共7971个，若不改变筛选差异表达基因的方法，则allFoldChange中的数据都是不变的。。不管下面选择分析网络图G1的结果还是 G2的结果，allFoldChange都是不变的
        Map<String,Map<String,Double>> allFoldChange = (Map<String,Map<String,Double>>) MySerialization.antiSerializeObject("D:\\soybean project\\soybean-differ-expression\\allFoldChange\\allFoldChange.obj");
        MyPrint.print("共:"+allFoldChange.size()+" 组fold-change数据");

        GetDifferExpGenesDataMatrix getDifferExpGenesDataMatrix = new GetDifferExpGenesDataMatrix();

        //先把找到的差异表达基因集合读取进来
        Set<String> differExpGenes = new HashSet<>();
        getDifferExpGenesDataMatrix.readDifferGenes(differExpGenes);

        //先分析网络图G1的划分结果，基因---所属社区关系
//        Partition bestNonOverlapPartition = (Partition) MySerialization.antiSerializeObject("D:\\paperdata\\soybean\\community detection\\历史计算结果\\2017.2.26\\bestNonOverlapPartition.obj");

        //网络图G2的社区划分结果
        Partition bestNonOverlapPartition = (Partition) MySerialization.antiSerializeObject("D:\\paperdata\\soybean\\community detection\\历史计算结果\\2017.3.9网络图G2\\bestNonOverlapPartition.obj");
        Set<String> genesInGraph = new HashSet<>();
        genesInGraph.addAll(bestNonOverlapPartition.getNodeCommunityMap().keySet());

        prepareDataForHeatmap(allFoldChange,genesInGraph,bestNonOverlapPartition);
        prepareDataForHeatmap(allFoldChange,differExpGenes);
    }

    /**
     *构造需要的fold-change数据
     * @param allFoldChange 全部差异表达基因的fold-change值
     * @param genes 想要构造的基因列表，可以是网络图G2或者G1中的基因
     * @param bestNonOverlapPartition
     */
    public static void prepareDataForHeatmap(Map<String,Map<String,Double>> allFoldChange ,Set<String>genes ,Partition bestNonOverlapPartition){
        try{
            FileWriter writer = new FileWriter("D:\\soybean project\\soybean-differ-expression\\allFoldChange\\allFoldChange.csv");
            BufferedWriter bw = new BufferedWriter(writer);

            FileWriter wr = new FileWriter("D:\\soybean project\\soybean-differ-expression\\allFoldChange\\geneMapComIndex.csv");
            BufferedWriter bw2 = new BufferedWriter(wr);

            //制作geneMapComIndex.csv的表头
            StringBuffer sb = new StringBuffer();
            sb.append("geneId");
            sb.append("\t");
            sb.append("communityIndex");
            sb.append("\t");
            sb.append("communityName");
            bw2.write(sb.toString());
            bw2.newLine();
            //准备一份“社区名---社区序号”的映射关系数据，所有社区size < 10的社区统一用序号0代表，其他的大社区从1开始安排序号
            Map<String,Integer> comNameMapIndex = new HashMap<>();
            int index=1;
            for(Map.Entry entry:bestNonOverlapPartition.getCommunities().entrySet()){
                String comName = (String) entry.getKey();
                List<String> community= (List<String>)entry.getValue();
                int size = community.size();

                if(size < 10){
                    comNameMapIndex.put(comName,0);//size < 10的社区统一都用序号0来代表
                }else {
                    comNameMapIndex.put(comName,index++);
                }
            }
            MyPrint.print("size >= 10的社区 即社区序号index="+(index-1));


            //制作allFoldChange.csv的表头
            List<String> allExperments = new ArrayList<>();
            allExperments.addAll(allFoldChange.keySet());
            sb = new StringBuffer();
            sb.append("geneId");
            for(String expermentType :allExperments){
                sb.append("\t");
                sb.append(expermentType);
            }
            bw.write(sb.toString());//打印表头
            bw.newLine();
            int test = 0;
            for(String gene :genes){//遍历所有差异表达基因
                sb = new StringBuffer();
                sb.append(gene);

                for(String expName: allExperments){//遍历所有实验条件，有的GSE下的数据我分成了2组实验来计算fold-change值
                    sb.append("\t");
                    sb.append(log2(allFoldChange.get(expName).get(gene)));
                }
                bw.write(sb.toString());
                bw.newLine();

                //下面开始构造geneMapComIndex.csv的数据
                if(bestNonOverlapPartition.getNodeCommunityMap().containsKey(gene)){
                    sb = new StringBuffer();
                    sb.append(gene);
                    sb.append("\t");
                    sb.append(comNameMapIndex.get(bestNonOverlapPartition.getNodeCommunityMap().get(gene)));//gene所属社区序号
                    sb.append("\t");
                    sb.append(bestNonOverlapPartition.getNodeCommunityMap().get(gene));//gene所属社区名
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
            MyPrint.print("输出heatmap所需的fold-change值异常");
            e.printStackTrace();
        }
    }

    /**
     *此方法构造出全部差异表达基因的fold-change数据
     * @param allFoldChange
     */
    public static void prepareDataForHeatmap(Map<String,Map<String,Double>> allFoldChange ,Set<String>differExpGenes ){
        try{
            FileWriter writer = new FileWriter("D:\\soybean project\\soybean-differ-expression\\allFoldChange\\allDifferGenesFoldChange.csv");
            BufferedWriter bw = new BufferedWriter(writer);

            //制作allFoldChange.csv的表头
            List<String> allExperments = new ArrayList<>();
            allExperments.addAll(allFoldChange.keySet());
            StringBuffer sb = new StringBuffer();
            sb.append("geneId");
            for(String expermentType :allExperments){
                sb.append("\t");
                sb.append(expermentType);
            }
            bw.write(sb.toString());//打印表头
            bw.newLine();
            int test = 0;
            for(String gene :differExpGenes){//遍历所有差异表达基因
                sb = new StringBuffer();
                sb.append(gene);

                for(String expName: allExperments){//遍历所有实验条件，有的GSE下的数据我分成了2组实验来计算fold-change值
                    sb.append("\t");
                    sb.append(log2(allFoldChange.get(expName).get(gene)));
                }
                bw.write(sb.toString());
                bw.newLine();

            }

            bw.close();
            writer.close();
        }catch (Exception e){
            MyPrint.print("输出heatmap所需的fold-change值异常");
            e.printStackTrace();
        }
    }

    private static double log2(double data){
        return Math.log(data)/Math.log(2);
    }
}
