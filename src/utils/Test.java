package utils;

import java.text.SimpleDateFormat;
import java.util.*;

/**
 * Created by Liyanzhen on 2016/12/18.
 */
public class Test {
    public static void main(String []args){
        int n=3,m=10,k=2;
        int[][] matrix={{4,6,8},{5,2,5}};
        int sum = find(m,k,n,matrix,0);
        MyPrint.print(sum+"");

        testDateTime();
        testClear();
    }
    public static int find(int m,int k,int n,int[][] matrix,int i){
        int sum = 0;
        if(i < k) {
            for (int j = 0; j < n; j++) {
                if (matrix[i][j] > m) {
                    sum += Math.pow(n, k - i - 1);
                } else {
                    sum += find(m - matrix[i][j], k, n, matrix, i + 1);
                }
            }
        }
        return sum;
    }

    /**
     * 获取时间差，测试方法
     */
    public static void testDateTime(){
        try {
            Date date = new Date();
            SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
            String startDate = sdf.format(date);

            Thread.sleep(6000);

            String endDate = sdf.format(new Date());
            int minutes = (int) (sdf.parse(endDate).getTime() - sdf.parse(startDate).getTime()) / (1000*60);
            MyPrint.print("当前时间差：" + minutes);
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    private static void testClear(){
        List<String> list = new ArrayList<>();
        Map<String,List<String>> map = new HashMap<>();
        for(int i=0;i < 3;i++){
            list.add(""+i);
        }
        map.put("first",list);
        for(String element: map.get("first")){
            MyPrint.print(element+"");
        }
        map.get("first").clear();
        MyPrint.print("clear后");
        for(String element: map.get("first")){
            MyPrint.print(element+"");
        }
    }
}
