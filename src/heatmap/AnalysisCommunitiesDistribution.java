package heatmap;

import serialprocess.OverlapPartition;
import serialprocess.Partition;
import utils.MyPrint;
import utils.MySerialization;
import utils.MyUtils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;

/**
 * Created by Liyanzhen on 2017/2/15.
 */
public class AnalysisCommunitiesDistribution {
    public static void main(String[] args){
        //反序列化之前保存的网络图G、最优非重叠划分、重叠划分结果
//        Graph g = (Graph) MySerialization.antiSerializeObject("D:\\paperdata\\soybean\\community detection\\original graph structure\\graph.obj");
        Partition bestNonOverlapPartition = (Partition) MySerialization.antiSerializeObject("D:\\paperdata\\soybean\\community detection\\最终结果\\bestNonOverlapPartition.obj");
        OverlapPartition overlapPartition = (OverlapPartition) MySerialization.antiSerializeObject("D:\\paperdata\\soybean\\community detection\\最终结果\\overlapPartition.obj");

        MyPrint.print("非重叠社区个数： "+bestNonOverlapPartition.getCommunities().size());
        MyPrint.print("重叠社区个数： "+overlapPartition.getCommunities().size());

        compare1(bestNonOverlapPartition,overlapPartition);
    }

    private static void compare1(Partition bestNonOverlapPartition, OverlapPartition overlapPartition){
        Set<String> comNames = new HashSet<String>();
        int[]nonOverlapComSize = new int[bestNonOverlapPartition.getCommunities().size()];
        int[]overlapComSize = new int[overlapPartition.getCommunities().size()];
        Iterator iter = bestNonOverlapPartition.getCommunities().keySet().iterator();
        while(iter.hasNext()){
            String comName =(String) iter.next();
            comNames.add(comName);
        }
        Set<String> totalSeperateCommuNames = new HashSet<String>();
        Set<String> meaningfulSeparateCommuNames = new HashSet<String>();
        Set<String> meaningfullOverlapCommuNames = new HashSet<String>();
        Set<String> totalOverlapTagSet = new HashSet<String>();
        int i=0;
        for(String commuName :comNames){
            nonOverlapComSize[i] = bestNonOverlapPartition.getCommunities().get(commuName).size();
            overlapComSize[i] = overlapPartition.getCommunities().get(commuName).size();
            List<String> nonOverlapCommu = bestNonOverlapPartition.getCommunities().get(commuName);
            List<String> overlapCommu = overlapPartition.getCommunities().get(commuName);
            if(overlapCommu.size() <= nonOverlapCommu.size()){
                //该社区中没有划入重叠节点, 该社区可能是独立的社区，先统计了
                totalSeperateCommuNames.add(commuName);

                if(nonOverlapCommu.size() >= 10){ //如果独立的社区内的基因数大于7，则该“独立”社区也具有一定的生物学意义，值得分析
                    meaningfulSeparateCommuNames.add(commuName);
                }
            }else{
                //有重叠节点被划进了commuName的社区
                totalOverlapTagSet.add(commuName);
                Set<String> overlapCommunities = new HashSet<String>();
                overlapCommunities.add(commuName);
                MyPrint.print("社区标志"+commuName+" 的社区 吸入重叠节点："+(overlapCommu.size()-nonOverlapCommu.size())+"个");
                MyPrint.print("社区标志"+commuName+"非重叠时大小="+nonOverlapCommu.size()+"; 重叠下大小="+overlapCommu.size());
                MyPrint.print("---------------------------------------");
                if( overlapCommu.size() >= 10){//若吸入节点之前 commuName的社区大小 > 5，吸入节点之后名为commuName的重叠社区大小 > 7，则认为该commuName对应的社区是有意义的
                    meaningfullOverlapCommuNames.add(commuName);
                }


                //莫要忽略了一点：被“吸入”的节点v，其原来所属社区c也就成了重叠社区，若c没有吸收新的节点，则在上面的统计中体现不出来社区c
                //因此还要找出来这些 “被吸入”的重叠节点原来所属社区，看是否同样值得分析
                List<String> overlapNodes = MyUtils.getOverlapNodes(overlapCommu,nonOverlapCommu);
                for(String node: overlapNodes){
                    List<String> commuTags = overlapPartition.getNodeMapCommunities().get(node);
                    //遍历节点node所属的多个社区，筛选出不属于当前社区commuName的 社区大小 > 7的社区
                    for(String tag: commuTags){
                        if(!tag.equals(commuName)){//即该tag的社区不是当前社区，
                            if(bestNonOverlapPartition.getCommunities().get(tag).size() >= 10) {//且大小 > 7，该社区是有意义的重叠社区之一
                                meaningfullOverlapCommuNames.add(tag);
                            }
                            totalOverlapTagSet.add(tag);

                        }
                    }
                }
            }
            i++;
        }

        MyPrint.print("-----------------------------");
        //此时找到的totalSeperateCommuNames 与 totalOverlapTagSet可能存在交集，交集中的社区名应该是重叠的，因此要找到交集，这些交集不再是独立社区了，因此要从totalSeperateCommuNames与meaningfulSeperateCommuNum中移除交集元素
        Set<String> intersection = getInterSectionOf2Set(totalSeperateCommuNames,totalOverlapTagSet);
        for(String e:intersection){
            totalSeperateCommuNames.remove(e);
            meaningfulSeparateCommuNames.remove(e);
        }


        MyPrint.print("独立社区 总数 = "+totalSeperateCommuNames.size());
        MyPrint.print("有意义的独立社区：“独立的”且所包含的基因个数 >= 10的社区 "+meaningfulSeparateCommuNames.size()+" 个");

        MyPrint.print("-----------------------------");
        MyPrint.print("总的重叠社区个数 = "+totalOverlapTagSet.size()+"；其中 "+meaningfullOverlapCommuNames.size()+"个有意义");
        //判断有意义的重叠社区中，有没有小社区完全被大社区所包含的情况，若有，则小社区应该被删除
        List<String> list1 = new ArrayList<String>();
        Set<String> removeList = new HashSet<String>();
        list1.addAll(meaningfullOverlapCommuNames);
        for(int m=0;m < list1.size();m++){
            for(int n=m+1;n < list1.size();n++){
                String mName = list1.get(m);
                String nName = list1.get(n);
                if(overlapPartition.getCommunities().get(mName).containsAll(overlapPartition.getCommunities().get(nName))){
                    removeList.add(nName);
                }else if(overlapPartition.getCommunities().get(nName).containsAll(overlapPartition.getCommunities().get(mName))){
                    removeList.add(mName);
                }
            }
        }
        MyPrint.print("有"+removeList.size()+"个社区被完全包含！！！！！！把被包含的小社区应该被移除");
        for(String e:removeList){
            meaningfullOverlapCommuNames.remove(e);
        }
        MyPrint.print("移除"+removeList.size()+"个被完全包含的小社区后，还有"+meaningfullOverlapCommuNames.size()+"个有意义的重叠社区");

        MyPrint.print("有意义的重叠社区，社区名如下--------------------------");
        for(String e :meaningfullOverlapCommuNames){
            MyPrint.print(e);
        }
        MyPrint.print("-----------------------------");

        //序列化保存 nonOverlapComSize数组和 overlapComSize数组
        MySerialization.serializeObject(nonOverlapComSize,"D:\\paperdata\\soybean\\community detection\\community analysis\\nonOverlapComSizeArray.obj");
        MySerialization.serializeObject(overlapComSize,"D:\\paperdata\\soybean\\community detection\\community analysis\\overlapComSizeArray.obj");

        //将有意义的重叠社区、 有意义的独立社区按格式打印到文本文件中
        saveMeaningfulComuunities(overlapPartition,meaningfullOverlapCommuNames,meaningfulSeparateCommuNames);
        saveTotalOverlapCommunities(overlapPartition,totalOverlapTagSet);
        findIneractionBetweenTwoverlapCommunities(overlapPartition,meaningfullOverlapCommuNames);
    }

    private static void findIneractionBetweenTwoverlapCommunities(OverlapPartition overlapPartition, Set<String> overlapTagSet){
        try {
            FileWriter writer = new FileWriter("D:\\paperdata\\soybean\\community detection\\community analysis\\overlapNodesBetweenMeaningfulOverlapCommunities.txt");
            BufferedWriter bw = new BufferedWriter(writer);
            Set<String> allOverlapNodes = new HashSet<String>();
            List<String> totalOverlapCommunityNames = new ArrayList<String>();
            int oc = overlapTagSet.size();
            String[] overlapCommunityNames = new String[oc];
            int[] trueOverlapComSizeArray = new int[oc];//保存True重叠社区的社区大小，与上面的社区名数组顺序一致
            int[][] intersectionMatrix = new int[oc][oc];
            totalOverlapCommunityNames.addAll(overlapTagSet);
            for (int i = 0; i < totalOverlapCommunityNames.size(); i++) {
                overlapCommunityNames[i] = totalOverlapCommunityNames.get(i);
                List<String> nodesOfCommunityI = overlapPartition.getCommunities().get(totalOverlapCommunityNames.get(i));
                Set<String> setI = new HashSet<String>();
                setI.addAll(nodesOfCommunityI);
                trueOverlapComSizeArray[i] = nodesOfCommunityI.size();
                for (int j = i + 1; j < totalOverlapCommunityNames.size(); j++) {
                    List<String> nodesOfCommunityJ = overlapPartition.getCommunities().get(totalOverlapCommunityNames.get(j));
                    Set<String> setJ = new HashSet<String>();
                    setJ.addAll(nodesOfCommunityJ);

                    Set<String> interactionOfCommunityIAndJ = getInterSectionOf2Set(setI, setJ);
                    if (interactionOfCommunityIAndJ.size() > 0) {
                        StringBuffer sb = new StringBuffer();
                        intersectionMatrix[i][j] = interactionOfCommunityIAndJ.size();
                        intersectionMatrix[j][i] = interactionOfCommunityIAndJ.size();
                        MyPrint.print("社区序号 ("+i+") 社区标志" + totalOverlapCommunityNames.get(i) + "(" + overlapPartition.getCommunities().get(totalOverlapCommunityNames.get(i)).size() + ")" + "与社区序号 ("+j+") 社区标志" + totalOverlapCommunityNames.get(j) + "(" + overlapPartition.getCommunities().get(totalOverlapCommunityNames.get(j)).size() + ")" + " 之间的重叠节点个数=" + interactionOfCommunityIAndJ.size());
                        sb.append("社区序号("+i+") 社区标志" + totalOverlapCommunityNames.get(i) + "(" + overlapPartition.getCommunities().get(totalOverlapCommunityNames.get(i)).size() + "个)" + "与社区序号 ("+j+") 与社区标志" + totalOverlapCommunityNames.get(j) + "(" + overlapPartition.getCommunities().get(totalOverlapCommunityNames.get(j)).size() + "个)" + " 之间的重叠节点个数=" + interactionOfCommunityIAndJ.size()+" :");
                        sb.append("\n");
                        for(String e: interactionOfCommunityIAndJ){
                            sb.append(e);
                            sb.append(",");
                        }
                        sb.deleteCharAt(sb.lastIndexOf(","));
                        bw.write(sb.toString());
                        bw.newLine();
                        bw.write("--------------------------");
                        bw.newLine();
                    }

//                    if((double)interactionOfCommunityIAndJ.size()/setI.size() < 0.85 && (double)interactionOfCommunityIAndJ.size()/setJ.size() < 0.85){
                        allOverlapNodes.addAll(interactionOfCommunityIAndJ);
//                    }

                }
            }

            bw.close();
            writer.close();

            MySerialization.serializeObject(overlapCommunityNames,"D:\\paperdata\\soybean\\community detection\\community analysis\\trueOverlapCommunitiesNames.obj");
            MySerialization.serializeObject(trueOverlapComSizeArray,"D:\\paperdata\\soybean\\community detection\\community analysis\\trueOverlapComSizeArray.obj");
            MySerialization.serializeObject(intersectionMatrix,"D:\\paperdata\\soybean\\community detection\\community analysis\\trueOverlapCommunitiesIntersection.obj");
            printOverlapNodes(allOverlapNodes);
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    private static void printOverlapNodes(Set<String> allOverlapNodes){
        try {
            Set<String> set = new HashSet<String>();
            set.add("Gma.16735.2.S1_at");
            set.add("GmaAffx.21211.1.S1_at");
            set.add("GmaAffx.92386.1.S1_at");
            set.add("GmaAffx.91805.1.S1_s_at");
            set.add("Gma.6606.1.S1_at");
            set.add("GmaAffx.91805.1.S1_at");
            set.add("GmaAffx.80951.1.S1_at");
            set.add("Gma.7559.1.S1_s_at");
            set.add("GmaAffx.92383.1.S1_at");
            Set<String> intersection = getInterSectionOf2Set(set,allOverlapNodes);
            MyPrint.print("+++++++++++++有意义的重叠社区内的重叠节点包含几个“文献体现的重要基因”："+intersection.size());

            FileWriter writer = new FileWriter("D:\\paperdata\\soybean\\community detection\\community analysis\\allOverlapNodes.txt");
            BufferedWriter bw = new BufferedWriter(writer);

            for(String e: allOverlapNodes){
                bw.write(e);
                bw.newLine();
            }

            bw.close();
            writer.close();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * 保存overlapPartition 中有意义的重叠社区以及有意义的“独立的”社区
     * @param overlapPartition
     * @param meaningfullOverlapCommuNames
     * @param meaningfulSeperateCommuNames
     */
    private static void saveMeaningfulComuunities(OverlapPartition overlapPartition, Set<String> meaningfullOverlapCommuNames, Set<String> meaningfulSeperateCommuNames){
        try {
            FileWriter writer1 = new FileWriter("D:\\paperdata\\soybean\\community detection\\community analysis\\meaningfulSeperateResult.txt");
            BufferedWriter bw1 = new BufferedWriter(writer1);
            FileWriter writer2 = new FileWriter("D:\\paperdata\\soybean\\community detection\\community analysis\\meaningfulOverlapResult.txt");
            BufferedWriter bw2 = new BufferedWriter(writer2);

            StringBuffer sb = new StringBuffer();
            sb.append("本文件中的社区是有意义的独立社区，即“独立的无重叠节点的”且所包含的基因个数 >= 10的社区 \n");
            bw1.write(sb.toString());
            sb = new StringBuffer();
            sb.append("有意义的社区个数 ="+meaningfulSeperateCommuNames.size()+"\n");
            bw1.newLine();
            bw1.write(sb.toString());
            //输出有意义的“独立社区”
            for(String name :meaningfulSeperateCommuNames){
                bw1.newLine();
                bw1.write("-----------------------------------------------------------------------------");
                bw1.newLine();
                bw1.write("社区名： "+name+"包含"+overlapPartition.getCommunities().get(name).size()+" 个基因");
                for(String e:overlapPartition.getCommunities().get(name)){
                    bw1.newLine();
                    bw1.write(e);
                }
                bw1.newLine();
                bw1.write("-----------------------------------------------------------------------------");
            }

            sb = new StringBuffer();
            sb.append("本文件中是真正的重叠社区 \n");
            bw2.write(sb.toString());
            sb = new StringBuffer();
            sb.append("其中有意义的重叠社区个数 = "+meaningfullOverlapCommuNames.size()+"\n");
            bw2.newLine();
            bw2.write(sb.toString());
            sb = new StringBuffer();
            sb.append("若吸入节点之后名为commuName的重叠社区大小 >= 10，则认为该commuName对应的社区是有意义的");
            bw2.newLine();
            bw2.write(sb.toString());
            sb = new StringBuffer();
            sb.append("那些原本只有1个或很少节点的社区，经过识别重叠节点，其社区大小增长很大；如原来只有1个节点，最后增长到了1297的社区，原来仅有的1个节点基因可能有生物学意义");
            bw2.newLine();
            bw2.write(sb.toString());

            //输出有意义的重叠社区
            List<String> list = new ArrayList<String>();
            list.addAll(meaningfullOverlapCommuNames);
            for(int i=0;i < list.size();i++){
                String comName = list.get(i);
                bw2.newLine();
                bw2.write("-----------------------------------------------------------------------------");
                bw2.newLine();
                bw2.write("社区名： "+comName+"包含"+overlapPartition.getCommunities().get(comName).size()+" 个基因");
                for(String e:overlapPartition.getCommunities().get(comName)){
                    bw2.newLine();
                    bw2.write(e);
                }
                bw2.newLine();
                bw2.write("-----------------------------------------------------------------------------");
            }

            bw1.close();
            writer1.close();
            bw2.close();
            writer2.close();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * txt保存全部的重叠社区
     * @param overlapPartition
     * @param totalOverlapTagSet
     */
    private static void saveTotalOverlapCommunities(OverlapPartition overlapPartition, Set<String> totalOverlapTagSet){
        try {
            FileWriter writer = new FileWriter("D:\\paperdata\\soybean\\community detection\\community analysis\\totalOverlapCommunitiesResult.txt");
            BufferedWriter bw = new BufferedWriter(writer);

            StringBuffer sb = new StringBuffer();
            sb.append("本文件中是全部的重叠社区，个数 = "+totalOverlapTagSet.size());
            bw.write(sb.toString());
            sb = new StringBuffer();
            Iterator it = overlapPartition.getCommunities().entrySet().iterator();
            while (it.hasNext()) {
                Map.Entry<String,List<String>> entry = (Map.Entry<String,List<String>>)it.next();
                String commuName = entry.getKey();
                if(totalOverlapTagSet.contains(commuName)){
                    bw.newLine();
                    bw.write("-----------------------------------------------------------------------------");
                    bw.newLine();
                    bw.write("社区名： "+commuName+"包含"+overlapPartition.getCommunities().get(commuName).size()+" 个基因");
                    for(String e:overlapPartition.getCommunities().get(commuName)){
                        bw.newLine();
                        bw.write(e);
                    }
                    bw.newLine();
                    bw.write("-----------------------------------------------------------------------------");
                }

            }
            bw.close();
            writer.close();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * 求两个集合的交集
     * @param totalSeperateCommuNames
     * @param totalOverlapCommuNames
     * @return 交集
     */
    private static Set<String> getInterSectionOf2Set(Set<String> totalSeperateCommuNames, Set<String> totalOverlapCommuNames){
        Set<String> result = new HashSet<String>();

        result.clear();
        result.addAll(totalSeperateCommuNames);
        result.retainAll(totalOverlapCommuNames);
        return result;
    }

    /**
     * 返回非重叠社区划分的社区size分布
     * @return
     */
    public int[] getNonoverlapCommunitiesSizeArray(){
        return (int[]) MySerialization.antiSerializeObject("D:\\paperdata\\soybean\\community detection\\community analysis\\nonOverlapComSizeArray.obj");
    }


    /**
     * 返回重叠社区划分结果的社区size分布
     * @return
     */
    public int[] getOverlapCommunitiesSizeArray(){
        return (int[]) MySerialization.antiSerializeObject("D:\\paperdata\\soybean\\community detection\\community analysis\\overlapComSizeArray.obj");
    }

    /**
     * “真正的”重叠社区，两两之间的交集大小
     * @return
     */
    public int[][] getTrueOverlapCommunitiesIntersectionMatrix(){
        return (int[][]) MySerialization.antiSerializeObject("D:\\paperdata\\soybean\\community detection\\community analysis\\trueOverlapCommunitiesIntersection.obj");
    }

    public String[] getTrueOverlapCommunitiesNames(){
        return (String[]) MySerialization.antiSerializeObject("D:\\paperdata\\soybean\\community detection\\community analysis\\trueOverlapCommunitiesNames.obj");
    }

    public int[] getTrueOverlapComSizeArray(){
        return (int[]) MySerialization.antiSerializeObject("D:\\paperdata\\soybean\\community detection\\community analysis\\trueOverlapComSizeArray.obj");
    }

}
